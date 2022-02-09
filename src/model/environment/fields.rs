/*
Copyright 2021, 2022 Jakub Lewandowski

This file is part of Parcel Ascent Tracing System (PATS).

Parcel Ascent Tracing System (PATS) is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Parcel Ascent Tracing System (PATS) is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Parcel Ascent Tracing System (PATS). If not, see https://www.gnu.org/licenses/.
*/

//! Sub-module responsible for handling
//! pressure level data buffering.
use super::finite_difference;
use super::interpolation::{precompute_tricubic_coefficients, Point3D};
use super::projection::LambertConicConformal;
use crate::model::{configuration, LonLat};
use crate::{
    errors::{EnvironmentError, InputError},
    model::{configuration::Input, environment::DomainExtent},
    Float,
};
use eccodes::{CodesHandle, FallibleIterator, ProductKind::GRIB};
use eccodes::{
    KeyType::{self, FloatArray, Int, Str},
    KeyedMessage,
};
use floccus::constants::G;
use log::debug;
use ndarray::{concatenate, s, stack, Array, Array2, Array3, Axis, Zip};
use rustc_hash::FxHashSet;

#[derive(Debug)]
struct RawFields {
    height: Array3<Float>,
    pressure: Array3<Float>,
    temperature: Array3<Float>,
    specific_humidity: Array3<Float>,
    u_wind: Array3<Float>,
    v_wind: Array3<Float>,
    virtual_temperature: Array3<Float>,
}

/// Struct for storing environmental variables
/// from levels above ground (currently pressure levels).
///
/// To limit IO operations and reduce performance overhead
/// of the model boundary conditions data is stored in the
/// memory as 3D arrays.
#[derive(Debug)]
pub struct Fields {
    pub lons: Array2<Float>,
    pub lats: Array2<Float>,
    pub height: Array3<Float>,

    pub temperature_coeffs: Array3<[Float; 64]>,
    pub pressure_coeffs: Array3<[Float; 64]>,
    pub u_wind_coeffs: Array3<[Float; 64]>,
    pub v_wind_coeffs: Array3<[Float; 64]>,
    pub spec_humidity_coeffs: Array3<[Float; 64]>,
    pub virtual_temp_coeffs: Array3<[Float; 64]>,
}

impl Fields {
    pub(super) fn new(
        input: &Input,
        domain_edges: DomainExtent<usize>,
        proj: &LambertConicConformal,
    ) -> Result<Fields, EnvironmentError> {
        let data = collect_data(input)?;
        let surfaces = construct(input, &data, domain_edges, proj)?;

        Ok(surfaces)
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn collect_data(input: &configuration::Input) -> Result<Vec<KeyedMessage>, InputError> {
    let mut data_levels: Vec<KeyedMessage> = vec![];

    for file in &input.data_files {
        let handle = CodesHandle::new_from_file(file, GRIB)?;

        let mut data: Vec<KeyedMessage> = handle
            .filter(|msg| {
                Ok(
                    msg.read_key("typeOfLevel")?.value == Str(input.level_type.clone())
                        && (msg.read_key("shortName")?.value == Str("z".to_string())
                            || msg.read_key("shortName")?.value == Str("q".to_string())
                            || msg.read_key("shortName")?.value == Str("t".to_string())
                            || msg.read_key("shortName")?.value == Str("u".to_string())
                            || msg.read_key("shortName")?.value == Str("v".to_string())),
                )
            })
            .collect()?;

        data_levels.append(&mut data);
    }

    if data_levels.is_empty() {
        return Err(InputError::DataNotSufficient(
            "Not enough variables on isobaric levels, check your input data",
        ));
    }

    Ok(data_levels)
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn construct(
    input: &Input,
    data: &[KeyedMessage],
    domain_edges: DomainExtent<usize>,
    proj: &LambertConicConformal,
) -> Result<Fields, EnvironmentError> {
    debug!("Buffering fields");

    // need a margin of one for derivate and coefficients computation later on
    let domain_edges = DomainExtent::<usize> {
        north: domain_edges.north - 1,
        south: domain_edges.south + 1,
        east: domain_edges.east + 1,
        west: domain_edges.west - 1,
    };

    let (lons, lats) = obtain_lonlat_fields_coords(&input.distinct_lonlats, domain_edges);
    let raw_fields = obtain_raw_fields(input, domain_edges, data)?;
    let fields = compute_fields_data(&raw_fields, &(lons, lats), proj);

    Ok(fields)
}

/// Buffers longitudes and latitudes of pressure level data gridpoints.
fn obtain_lonlat_fields_coords(
    distinct_lonlats: &(Vec<Float>, Vec<Float>),
    domain_edges: DomainExtent<usize>,
) -> (Array2<Float>, Array2<Float>) {
    let lats = distinct_lonlats.1[domain_edges.north..=domain_edges.south].to_vec();
    let lons;

    if domain_edges.west < domain_edges.east {
        lons = distinct_lonlats.0[domain_edges.west..=domain_edges.east].to_vec();
    } else {
        let left_half = &distinct_lonlats.0[domain_edges.east..];
        let right_half = &distinct_lonlats.0[..=domain_edges.west];

        lons = [left_half, right_half].concat();
    }

    let lons = Array::from_vec(lons);
    let lats = Array::from_vec(lats);

    let lons_view = vec![lons.view(); lats.len()];
    let lats_view = vec![lats.view(); lons.len()];

    let lons = stack(Axis(1), lons_view.as_slice()).unwrap();
    let lats = stack(Axis(0), lats_view.as_slice()).unwrap();

    (lons, lats)
}

/// Reads variables on pressure levels from GRIB file
/// and buffers them in domain + margins extent.
fn obtain_raw_fields(
    input: &Input,
    domain_edges: DomainExtent<usize>,
    data: &[KeyedMessage],
) -> Result<RawFields, InputError> {
    let input_shape = input.shape;

    let pressure = read_truncated_pressure(data, domain_edges)?;

    let geopotential = read_raw_field("z", input_shape, data)?;
    let height = truncate_field_to_extent(&geopotential, domain_edges).mapv(|v| v / G);

    let temperature = read_raw_field("t", input_shape, data)?;
    let temperature = truncate_field_to_extent(&temperature, domain_edges);

    let u_wind = read_raw_field("u", input_shape, data)?;
    let u_wind = truncate_field_to_extent(&u_wind, domain_edges);

    let v_wind = read_raw_field("v", input_shape, data)?;
    let v_wind = truncate_field_to_extent(&v_wind, domain_edges);

    // specific humidity needs a check for negative values
    // which are replaced with the smallest positive value
    let specific_humidity = read_raw_field("q", input_shape, data)?;
    let specific_humidity = truncate_field_to_extent(&specific_humidity, domain_edges).mapv(|v| {
        if v <= 0.0 {
            1.0e-8
        } else {
            v
        }
    });

    let virtual_temperature = compute_vtemp_field(&temperature, &specific_humidity);

    Ok(RawFields {
        height,
        pressure,
        temperature,
        specific_humidity,
        u_wind,
        v_wind,
        virtual_temperature,
    })
}

/// Creates a 3d array of pressure data of shape
/// identical to other pressure level fields.
///
/// When data in GRIB file is provided on pressure levels
/// the information about pressure at each level is only
/// stored in message metadata. Thus, information what
/// pressure levels are available has to be extracted
/// and then casted to the 3d array expected by the
/// [`accesser`](super::super::accesser).
fn read_truncated_pressure(
    levels_data: &[KeyedMessage],
    domain_edges: DomainExtent<usize>,
) -> Result<Array3<Float>, InputError> {
    let xy_shape = (
        (domain_edges.east as isize - domain_edges.west as isize).abs() as usize + 1,
        (domain_edges.south as isize - domain_edges.north as isize).abs() as usize + 1,
    );

    let levels_list = list_levels(levels_data)?;
    let mut pressure_levels = vec![];

    for level in levels_list {
        let pressure_level = Array2::from_elem(xy_shape, level);
        let pressure_level = pressure_level.mapv(|v| (v as Float) * 100.0);
        pressure_levels.push(pressure_level);
    }

    let mut pressure_views = vec![];

    for level in &pressure_levels {
        pressure_views.push(level.view());
    }

    let pressure_levels = ndarray::stack(Axis(0), pressure_views.as_slice()).unwrap();

    Ok(pressure_levels)
}

/// Function to get the list of unique levels
/// of specified type in the provided GRIB files.
///
/// This function uses `FxHashSet` to get unique list of levels,
/// which in benchmarsk showed outstanding performance.
fn list_levels(data: &[KeyedMessage]) -> Result<Vec<i64>, InputError> {
    debug!("Getting levels list");

    let mut unique_levels: FxHashSet<i64> = FxHashSet::default();

    for msg in data {
        let level_id = msg.read_key("level")?;

        if let KeyType::Int(id) = level_id.value {
            unique_levels.insert(id);
        } else {
            return Err(InputError::IncorrectKeyType("level"));
        };
    }

    let mut unique_levels: Vec<i64> = unique_levels.into_iter().collect();
    unique_levels.sort_unstable();
    unique_levels.reverse();

    Ok(unique_levels)
}

/// Reads all values in GRIB file at specified level type
/// of variable with given `short_name` and collects them
/// into a 3d array.
///
/// In GRIB files data on each level is stored as a separate
/// message. Those need to be collected and only then can be
/// converted into a 3d array.
fn read_raw_field(
    short_name: &str,
    shape: (usize, usize),
    data: &[KeyedMessage],
) -> Result<Array3<Float>, InputError> {
    let data_levels = read_raw_messages(short_name, data)?;
    let result_data = messages_to_array(data_levels, shape)?;

    Ok(result_data)
}

/// Filters and read all GRIB messages that contain
/// variable with given `short_name` on specified level type.
fn read_raw_messages<'a>(
    short_name: &str,
    data: &'a [KeyedMessage],
) -> Result<Vec<&'a KeyedMessage>, InputError> {
    let mut data_levels: Vec<&KeyedMessage> = vec![];

    for msg in data {
        if msg.read_key("shortName")?.value == Str(short_name.to_string()) {
            data_levels.push(msg);
        }
    }

    if data_levels.is_empty() {
        return Err(InputError::DataNotSufficient(
            "Not enough data at isobaric levels",
        ));
    }

    Ok(data_levels)
}

/// Collects data from GRIB messages on specified level type
/// into a 3d array,
fn messages_to_array(
    data_levels: Vec<&KeyedMessage>,
    shape: (usize, usize),
) -> Result<Array3<Float>, InputError> {
    let mut sorted_data_levels = vec![];

    for msg in data_levels {
        let lvl_id = if let Int(id) = msg.read_key("level")?.value {
            id
        } else {
            return Err(InputError::IncorrectKeyType("level"));
        };

        let lvl_vals = if let FloatArray(vals) = msg.read_key("values")?.value {
            vals
        } else {
            return Err(InputError::IncorrectKeyType("values"));
        };

        // data values in GRIB are a vec of values row-by-row (x-axis is in WE direction)
        // we want a Array2 of provided `shape` with x-axis in WE direction
        // but from_shape_vec(final_shape, data) splits the data into final_shape.0 long chunks
        // and stacks them along Axis(0)
        // we then transpose that array to get axes along expected geographical directions
        // this solution is against contiguos memory convention, for historical reasons
        let lvl_vals = Array2::from_shape_vec((shape.1, shape.0), lvl_vals)?;
        let lvl_vals = lvl_vals.reversed_axes();

        sorted_data_levels.push((lvl_id, lvl_vals));
    }

    sorted_data_levels.sort_unstable_by_key(|k| k.0);
    sorted_data_levels.reverse();

    let sorted_data_levels: Vec<Array2<f64>> =
        sorted_data_levels.into_iter().map(|t| t.1).collect();

    let mut result_data = vec![];
    for lvl in &sorted_data_levels {
        result_data.push(lvl.view());
    }

    let result_data = ndarray::stack(Axis(0), result_data.as_slice())?;
    let result_data = result_data.mapv(|v| v as Float);

    Ok(result_data)
}

/// Truncates data on specified level type from GRIB file
/// to cover only the message + margins extent.
fn truncate_field_to_extent(
    raw_field: &Array3<Float>,
    domain_edges: DomainExtent<usize>,
) -> Array3<Float> {
    // truncate in NS axis
    let truncated_field = raw_field.slice(s![.., .., domain_edges.north..=domain_edges.south]);

    if domain_edges.west < domain_edges.east {
        let truncated_field =
            truncated_field.slice(s![.., domain_edges.west..=domain_edges.east, ..]);
        return truncated_field.to_owned();
    }

    let left_half = truncated_field.slice(s![.., domain_edges.west.., ..]);
    let right_half = truncated_field.slice(s![.., ..=domain_edges.east, ..]);
    let truncated_field = concatenate![Axis(1), left_half, right_half];
    truncated_field.to_owned()
}

/// Computes and buffers additional pressure level data from
/// values previously read from the GRIB file.
fn compute_vtemp_field(
    temperature: &Array3<Float>,
    specific_humidity: &Array3<Float>,
) -> Array3<Float> {
    let mut virtual_temperature: Array3<Float> = Array3::zeros(temperature.raw_dim());

    Zip::from(&mut virtual_temperature)
        .and(temperature)
        .and(specific_humidity)
        .for_each(|tv, &t, &q| {
            *tv = floccus::virtual_temperature::general3(t, q).expect(
                "Error while computing virtual temperature: variable out of reasonable bounds",
            );
        });

    virtual_temperature
}

fn compute_fields_data(
    raw_fields: &RawFields,
    coords: &LonLat<Array2<Float>>,
    proj: &LambertConicConformal,
) -> Fields {
    let coords_xy = project_lonlats(&coords, proj, raw_fields.height.shape()[0]);

    // compute derivatives
    let pressure_data_points = finite_difference::compute_3d_points(
        &raw_fields.pressure,
        &coords_xy.0,
        &coords_xy.1,
        &raw_fields.height,
    );

    let temperature_data_points = finite_difference::compute_3d_points(
        &raw_fields.temperature,
        &coords_xy.0,
        &coords_xy.1,
        &raw_fields.height,
    );

    let spec_humidity_data_points = finite_difference::compute_3d_points(
        &raw_fields.specific_humidity,
        &coords_xy.0,
        &coords_xy.1,
        &raw_fields.height,
    );

    let uwind_data_points = finite_difference::compute_3d_points(
        &raw_fields.u_wind,
        &coords_xy.0,
        &coords_xy.1,
        &raw_fields.height,
    );

    let vwind_data_points = finite_difference::compute_3d_points(
        &raw_fields.v_wind,
        &coords_xy.0,
        &coords_xy.1,
        &raw_fields.height,
    );

    let virtual_temp_data_points = finite_difference::compute_3d_points(
        &raw_fields.virtual_temperature,
        &coords_xy.0,
        &coords_xy.1,
        &raw_fields.height,
    );

    // compute coefficients
    let pressure_coeffs = compute_field_coeffs(&pressure_data_points);
    let temperature_coeffs = compute_field_coeffs(&temperature_data_points);
    let spec_humidity_coeffs = compute_field_coeffs(&spec_humidity_data_points);
    let u_wind_coeffs = compute_field_coeffs(&uwind_data_points);
    let v_wind_coeffs = compute_field_coeffs(&vwind_data_points);
    let virtual_temp_coeffs = compute_field_coeffs(&virtual_temp_data_points);

    Fields {
        lons: coords.0.slice(s![1..-1, 1..-1]).to_owned(),
        lats: coords.1.slice(s![1..-1, 1..-1]).to_owned(),
        height: raw_fields.height.slice(s![.., 1..-1, 1..-1]).to_owned(),
        temperature_coeffs,
        pressure_coeffs,
        u_wind_coeffs,
        v_wind_coeffs,
        spec_humidity_coeffs,
        virtual_temp_coeffs,
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn project_lonlats(
    lonlats: &LonLat<Array2<Float>>,
    proj: &LambertConicConformal,
    levels_count: usize,
) -> (Array3<Float>, Array3<Float>) {
    let mut x = Array2::default(lonlats.0.raw_dim());
    let mut y = Array2::default(lonlats.1.raw_dim());

    Zip::from(&mut x)
        .and(&mut y)
        .and(&lonlats.0)
        .and(&lonlats.1)
        .for_each(|x, y, &lon, &lat| {
            let projected_xy = proj.project(lon, lat);
            *x = projected_xy.0;
            *y = projected_xy.1;
        });

    let x = ndarray::stack(Axis(0), vec![x.view(); levels_count].as_slice()).unwrap();
    let y = ndarray::stack(Axis(0), vec![y.view(); levels_count].as_slice()).unwrap();

    (x, y)
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn compute_field_coeffs(points: &Array3<Point3D>) -> Array3<[Float; 64]> {
    let coeffs_shape = (points.dim().0 - 1, points.dim().1 - 1, points.dim().2 - 1);

    let mut points_lower: Array3<[Point3D; 4]> = Array3::default(coeffs_shape);
    let mut points_upper: Array3<[Point3D; 4]> = Array3::default(coeffs_shape);
    let mut coeffs = Array3::from_elem(coeffs_shape, [0.0; 64]);

    Zip::from(&mut points_lower)
        .and(&points.slice(s![0..-1, 0..-1, 0..-1]))
        .and(&points.slice(s![0..-1, 1.., 0..-1]))
        .and(&points.slice(s![0..-1, 0..-1, 1..]))
        .and(&points.slice(s![0..-1, 1.., 1..]))
        .for_each(|points, &p1, &p2, &p3, &p4| {
            *points = [p1, p2, p3, p4];
        });

    Zip::from(&mut points_upper)
        .and(&points.slice(s![1.., 0..-1, 0..-1]))
        .and(&points.slice(s![1.., 1.., 0..-1]))
        .and(&points.slice(s![1.., 0..-1, 1..]))
        .and(&points.slice(s![1.., 1.., 1..]))
        .for_each(|points, &p1, &p2, &p3, &p4| {
            *points = [p1, p2, p3, p4];
        });

    Zip::from(&mut coeffs)
        .and(&points_lower)
        .and(&points_upper)
        .for_each(|coeffs, &pl, &pu| {
            *coeffs = precompute_tricubic_coefficients(&[
                pl[0], pl[1], pl[2], pl[3], pu[0], pu[1], pu[2], pu[3],
            ]);
        });

    coeffs
}
