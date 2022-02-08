/*
Copyright 2021 Jakub Lewandowski

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
//! surface data buffering.

use super::{
    finite_difference,
    interpolation::{self, Point2D},
    projection::LambertConicConformal,
};
use crate::{
    errors::{EnvironmentError, InputError},
    model::{configuration::Input, environment::DomainExtent, LonLat},
    Float,
};
use eccodes::{CodesHandle, FallibleIterator, ProductKind::GRIB};
use eccodes::{
    KeyType::{FloatArray, Str},
    KeyedMessage,
};
use floccus::constants::G;
use log::debug;
use ndarray::{concatenate, s, stack, Array, Array2, Axis, Zip};

/// Struct for storing environmental variables at/near surface.
///
/// To limit IO operations and reduce performance overhead
/// of the model surface data is stored in the
/// memory as 2D arrays.
#[derive(Debug)]
pub struct Surfaces {
    pub lons: Array2<Float>,
    pub lats: Array2<Float>,
    pub height_coeffs: Array2<[Float; 16]>,
    pub pressure_coeffs: Array2<[Float; 16]>,
    pub temperature_coeffs: Array2<[Float; 16]>,
    pub dewpoint_coeffs: Array2<[Float; 16]>,
    pub u_wind_coeffs: Array2<[Float; 16]>,
    pub v_wind_coeffs: Array2<[Float; 16]>,
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Debug)]
struct RawSurfaces {
    height: Array2<Float>,
    pressure: Array2<Float>,
    temperature: Array2<Float>,
    dewpoint: Array2<Float>,
    u_wind: Array2<Float>,
    v_wind: Array2<Float>,
}

impl Surfaces {
    pub(super) fn new(
        input: &Input,
        domain_edges: DomainExtent<usize>,
        proj: &LambertConicConformal,
    ) -> Result<Surfaces, EnvironmentError> {
        let data = collect_data(input)?;
        let surfaces = construct(input, &data, domain_edges, proj)?;

        Ok(surfaces)
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn collect_data(input: &Input) -> Result<Vec<KeyedMessage>, InputError> {
    let mut data_levels: Vec<KeyedMessage> = vec![];

    for file in &input.data_files {
        let handle = CodesHandle::new_from_file(file, GRIB)?;

        let mut data: Vec<KeyedMessage> = handle
            .filter(|msg| {
                Ok(
                    msg.read_key("typeOfLevel")?.value == Str("surface".to_string())
                        && (msg.read_key("shortName")?.value == Str("10u".to_string())
                            || msg.read_key("shortName")?.value == Str("10v".to_string())
                            || msg.read_key("shortName")?.value == Str("2t".to_string())
                            || msg.read_key("shortName")?.value == Str("2d".to_string())
                            || msg.read_key("shortName")?.value == Str("sp".to_string())
                            || msg.read_key("shortName")?.value == Str("z".to_string())),
                )
            })
            .collect()?;

        data_levels.append(&mut data);
    }

    if data_levels.is_empty() {
        return Err(InputError::DataNotSufficient(
            "Not enough variables on surface levels, check your input data",
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
) -> Result<Surfaces, EnvironmentError> {
    debug!("Buffering surfaces");

    // need a margin of one for derivate and coefficients computation later on
    let domain_edges = DomainExtent::<usize> {
        north: domain_edges.north - 1,
        south: domain_edges.south + 1,
        east: domain_edges.east + 1,
        west: domain_edges.west - 1,
    };

    let (lons, lats) = obtain_lonlat_surface_coords(&input.distinct_lonlats, domain_edges);
    let raw_surfaces = obtain_raw_surfaces(input, data, domain_edges)?;
    let surfaces = compute_surfaces_data(raw_surfaces, (lons, lats), proj);

    Ok(surfaces)
}

/// Obtains longitudes and latitudes of surface data gridpoints.
fn obtain_lonlat_surface_coords(
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

    let lons = stack(Axis(1), vec![lons.view(); lats.len()].as_slice()).unwrap();
    let lats = stack(Axis(0), vec![lats.view(); lons.len()].as_slice()).unwrap();

    (lons, lats)
}

/// Reads variables on surface level from GRIB file
/// and buffers them in domain + margins extent.
fn obtain_raw_surfaces(
    input: &Input,
    data: &[KeyedMessage],
    domain_edges: DomainExtent<usize>,
) -> Result<RawSurfaces, InputError> {
    let input_shape = input.shape;

    let geopotential = read_raw_surface("z", input_shape, data)?;
    let height = truncate_surface(&geopotential, domain_edges).mapv(|v| v / G);

    let pressure = read_raw_surface("sp", input_shape, data)?;
    let pressure = truncate_surface(&pressure, domain_edges);

    let temperature = read_raw_surface("2t", input_shape, data)?;
    let temperature = truncate_surface(&temperature, domain_edges);

    let dewpoint = read_raw_surface("2d", input_shape, data)?;
    let dewpoint = truncate_surface(&dewpoint, domain_edges);

    let u_wind = read_raw_surface("10u", input_shape, data)?;
    let u_wind = truncate_surface(&u_wind, domain_edges);

    let v_wind = read_raw_surface("10v", input_shape, data)?;
    let v_wind = truncate_surface(&v_wind, domain_edges);

    Ok(RawSurfaces {
        height,
        pressure,
        temperature,
        dewpoint,
        u_wind,
        v_wind,
    })
}

/// Reads all values in GRIB file at surface level
/// of variable with given `short_name`.
fn read_raw_surface(
    short_name: &str,
    shape: (usize, usize),
    data: &[KeyedMessage],
) -> Result<Array2<Float>, InputError> {
    let mut data_level = None;

    for msg in data {
        if msg.read_key("shortName")?.value == Str(short_name.to_string()) {
            data_level = Some(msg);
            break;
        }
    }

    if data_level.is_none() {
        return Err(InputError::DataNotSufficient(
            "Not enough data on surface levels, check your input data",
        ));
    }

    let data_level = data_level.unwrap();
    let data_level = data_level.read_key("values")?.value;
    let data_level = if let FloatArray(v) = data_level {
        v
    } else {
        return Err(InputError::IncorrectKeyType("values"));
    };

    // data values in GRIB are a vec of values row-by-row (x-axis is in WE direction)
    // we want a Array2 of provided `shape` with x-axis in WE direction
    // but from_shape_vec(final_shape, data) splits the data into final_shape.0 long chunks
    // and stacks them along Axis(0)
    // we then transpose that array to get axes along expected geographical directions
    // this solution is against contiguos memory convention, for historical reasons
    let result_data = Array2::from_shape_vec((shape.1, shape.0), data_level)?;
    let result_data = result_data.reversed_axes();
    let result_data = result_data.mapv(|v| v as Float);

    Ok(result_data)
}

/// Truncates surface data array from GRIB file to
/// cover only the domain + margins extent.
fn truncate_surface(raw_field: &Array2<Float>, domain_edges: DomainExtent<usize>) -> Array2<Float> {
    // truncate in NS axis
    let truncated_field = raw_field.slice(s![.., domain_edges.north..=domain_edges.south]);

    if domain_edges.west < domain_edges.east {
        let truncated_field = truncated_field.slice(s![domain_edges.west..=domain_edges.east, ..]);
        return truncated_field.to_owned();
    }

    let left_half = truncated_field.slice(s![domain_edges.west.., ..]);
    let right_half = truncated_field.slice(s![..=domain_edges.east, ..]);
    let truncated_field = concatenate![Axis(0), left_half, right_half];
    truncated_field.to_owned()
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn compute_surfaces_data(
    raw_surfaces: RawSurfaces,
    coords: LonLat<Array2<Float>>,
    proj: &LambertConicConformal,
) -> Surfaces {
    let coords_xy = project_lonlats(&coords, proj);

    // compute derivatives
    let height_data_points =
        finite_difference::compute_2d_points(&raw_surfaces.height, &coords_xy.0, &coords_xy.1);

    let pressure_data_points =
        finite_difference::compute_2d_points(&raw_surfaces.pressure, &coords_xy.0, &coords_xy.1);

    let temperature_data_points =
        finite_difference::compute_2d_points(&raw_surfaces.temperature, &coords_xy.0, &coords_xy.1);

    let dewpoint_data_points =
        finite_difference::compute_2d_points(&raw_surfaces.dewpoint, &coords_xy.0, &coords_xy.1);

    let u_wind_data_points =
        finite_difference::compute_2d_points(&raw_surfaces.u_wind, &coords_xy.0, &coords_xy.1);

    let v_wind_data_points =
        finite_difference::compute_2d_points(&raw_surfaces.v_wind, &coords_xy.0, &coords_xy.1);

    // compute coefficients
    let height_coeffs = compute_surface_coeffs(height_data_points);
    let pressure_coeffs = compute_surface_coeffs(pressure_data_points);
    let temperature_coeffs = compute_surface_coeffs(temperature_data_points);
    let dewpoint_coeffs = compute_surface_coeffs(dewpoint_data_points);
    let u_wind_coeffs = compute_surface_coeffs(u_wind_data_points);
    let v_wind_coeffs = compute_surface_coeffs(v_wind_data_points);

    // save coefficients to struct

    Surfaces {
        lons: coords.0.slice(s![1..-1, 1..-1]).to_owned(),
        lats: coords.1.slice(s![1..-1, 1..-1]).to_owned(),
        height_coeffs,
        pressure_coeffs,
        temperature_coeffs,
        dewpoint_coeffs,
        u_wind_coeffs,
        v_wind_coeffs,
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn project_lonlats(
    lonlats: &LonLat<Array2<Float>>,
    proj: &LambertConicConformal,
) -> (Array2<Float>, Array2<Float>) {
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

    (x, y)
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn compute_surface_coeffs(points: Array2<Point2D>) -> Array2<[Float; 16]> {
    let mut coeffs = Array2::default((points.dim().0 - 1, points.dim().1 - 1));

    Zip::from(&mut coeffs)
        .and(&points.slice(s![0..-1, 0..-1]))
        .and(&points.slice(s![1.., 0..-1]))
        .and(&points.slice(s![0..-1, 1..]))
        .and(&points.slice(s![1.., 1..]))
        .for_each(|coeffs, &p1, &p2, &p3, &p4| {
            *coeffs = interpolation::precompute_bicubic_coefficients([p1, p2, p3, p4]);
        });

    coeffs
}
