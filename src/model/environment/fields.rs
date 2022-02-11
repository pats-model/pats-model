/*
Copyright 2021 - 2022 Jakub Lewandowski

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

    pub temperature: Array3<Float>,
    pub pressure: Array3<Float>,
    pub u_wind: Array3<Float>,
    pub v_wind: Array3<Float>,
    pub spec_humidity: Array3<Float>,
    pub virtual_temp: Array3<Float>,
    pub vertical_vel: Array3<Float>,
}

impl Fields {
    pub(super) fn new(
        input: &Input,
        domain_edges: DomainExtent<usize>,
    ) -> Result<Self, EnvironmentError> {
        let data = collect(input)?;
        let fields = construct_fields(input, &data, domain_edges)?;

        Ok(fields)
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub(super) fn collect(input: &configuration::Input) -> Result<Vec<KeyedMessage>, InputError> {
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
                            || msg.read_key("shortName")?.value == Str("v".to_string())
                            || msg.read_key("shortName")?.value == Str("w".to_string())),
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

/// Function to read pressure level data from GRIB input
/// in extent covering domain and margins and buffer it.
///
/// Data from GRIB files can only be read in a whole
/// GRIB domain, therefore it needs to be truncated.
///
/// Some useful variables are not provided by
/// most NWP models, therefore they need to be computed.
fn construct_fields(
    input: &Input,
    data: &[KeyedMessage],
    domain_edges: DomainExtent<usize>,
) -> Result<Fields, EnvironmentError> {
    debug!("Buffering fields");

    let coords = cast_lonlat_fields_coords(&input.distinct_lonlats, domain_edges);
    let fields = assign_fields(input, domain_edges, data, coords)?;

    Ok(fields)
}

/// Buffers longitudes and latitudes of pressure level data gridpoints.
fn cast_lonlat_fields_coords(
    distinct_lonlats: &(Vec<Float>, Vec<Float>),
    domain_edges: DomainExtent<usize>,
) -> LonLat<Array2<Float>> {
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
fn assign_fields(
    input: &Input,
    domain_edges: DomainExtent<usize>,
    data: &[KeyedMessage],
    coords: LonLat<Array2<Float>>,
) -> Result<Fields, InputError> {
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

    let spec_humidity = read_raw_field("q", input_shape, data)?;
    // check for negative values of specific humidity and replace them with the smallest positive value
    let spec_humidity = truncate_field_to_extent(&spec_humidity, domain_edges).mapv(|v| {
        if v <= 0.0 {
            1.0e-8
        } else {
            v
        }
    });

    let virtual_temp = compute_virtual_temperature(&temperature, &spec_humidity);

    let vertical_motion = read_raw_field("w", input_shape, data)?;
    let vertical_motion = truncate_field_to_extent(&vertical_motion, domain_edges);
    let vertical_vel = compute_vertical_velocity(&pressure, &height, &vertical_motion);

    Ok(Fields {
        lons: coords.0,
        lats: coords.1,
        height,
        temperature,
        pressure,
        u_wind,
        v_wind,
        spec_humidity,
        virtual_temp,
        vertical_vel,
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

        // a bit of magic
        // data values in GRIB are a vec of values row-by-row (x-axis is in WE direction)
        // we want a Array2 of provided `shape` with x-axis in WE direction
        // but from_shape_vec(final_shape, data) splits the data into final_shape.1 long chunks
        // and puts them in columns
        // so we need to correctly split the data in GRIB vector into Array2 and then transpose
        // that array to get axes along expected geographical directions
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
    //truncate in NS axis
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
fn compute_virtual_temperature(
    temperature: &Array3<Float>,
    spec_humidity: &Array3<Float>,
) -> Array3<Float> {
    let mut virtual_temperature: Array3<Float> = Array3::zeros(temperature.raw_dim());

    Zip::from(&mut virtual_temperature)
        .and(temperature)
        .and(spec_humidity)
        .for_each(|tv, &t, &q| {
            *tv = floccus::virtual_temperature::general3(t, q).expect(
                "Error while computing virtual temperature: variable out of reasonable bounds",
            );
        });

    virtual_temperature
}

/// What it is?
fn compute_vertical_velocity(
    pressure: &Array3<Float>,
    height: &Array3<Float>,
    vertical_motion: &Array3<Float>,
) -> Array3<Float> {
    // compute thickness in negative m Pa^-1
    let mut thickness = (&height.slice(s![1.., .., ..]) - &height.slice(s![0..-1, .., ..]))
        / (&pressure.slice(s![1.., .., ..]) - &pressure.slice(s![0..-1, .., ..]));

    // thickness array doesn't have the top level, so we will copy it
    let thickness_top = thickness.slice(s![-1, .., ..]).to_owned();
    let thickness_top = Array2::from(thickness_top);
    thickness.push(Axis(0), thickness_top.view()).unwrap();

    // multiply vertical motion and thickness to get velocity
    vertical_motion * thickness
}
