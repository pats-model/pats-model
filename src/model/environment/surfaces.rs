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
//! surface data buffering.

use crate::model::{configuration, LonLat};
use crate::{
    errors::{EnvironmentError, InputError},
    model::{configuration::Input, environment::DomainExtent},
    Float,
};
use eccodes::{CodesHandle, FallibleIterator, ProductKind::GRIB};
use eccodes::{
    KeyType::{FloatArray, Str},
    KeyedMessage,
};
use floccus::constants::G;
use log::debug;
use ndarray::{concatenate, s, stack, Array, Array2, Axis};

/// Struct for storing environmental variables at/near surface.
///
/// To limit IO operations and reduce performance overhead
/// of the model surface data is stored in the
/// memory as 2D arrays.
#[derive(Debug)]
pub struct Surfaces {
    pub lons: Array2<Float>,
    pub lats: Array2<Float>,

    pub temperature: Array2<Float>,
    pub dewpoint: Array2<Float>,
    pub pressure: Array2<Float>,
    pub height: Array2<Float>,
    pub u_wind: Array2<Float>,
    pub v_wind: Array2<Float>,
}

impl Surfaces {
    pub(super) fn new(
        input: &Input,
        domain_edges: DomainExtent<usize>,
    ) -> Result<Self, EnvironmentError> {
        let data = collect(input)?;
        let surfs = construct_surfaces(input, &data, domain_edges)?;

        Ok(surfs)
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn collect(input: &configuration::Input) -> Result<Vec<KeyedMessage>, InputError> {
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

/// Function to read surface data from GRIB input
/// in extent covering domain and margins and buffer it.
///
/// Data from GRIB files can only be read in a whole
/// GRIB domain, therefore it needs to be truncated.
///
/// Some useful surface variables are not provided by
/// most NWP models, therefore they need to be computed.
fn construct_surfaces(
    input: &Input,
    data: &[KeyedMessage],
    domain_edges: DomainExtent<usize>,
) -> Result<Surfaces, EnvironmentError> {
    debug!("Buffering surfaces");

    let coords = cast_lonlat_surface_coords(&input.distinct_lonlats, domain_edges);
    let surfaces = assign_surfaces(input, data, domain_edges, coords)?;

    Ok(surfaces)
}

/// Buffers longitudes and latitudes of surface data gridpoints.
fn cast_lonlat_surface_coords(
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

/// Reads variables on surface level from GRIB file
/// and buffers them in domain + margins extent.
fn assign_surfaces(
    input: &Input,
    data: &[KeyedMessage],
    domain_edges: DomainExtent<usize>,
    coords: LonLat<Array2<Float>>,
) -> Result<Surfaces, InputError> {
    let input_shape = input.shape;

    let geopotential = read_raw_surface("z", input_shape, data)?;
    let height = truncate_surface_to_extent(&geopotential, domain_edges).mapv(|v| v / G);

    let pressure = read_raw_surface("sp", input_shape, data)?;
    let pressure = truncate_surface_to_extent(&pressure, domain_edges);

    let temperature = read_raw_surface("2t", input_shape, data)?;
    let temperature = truncate_surface_to_extent(&temperature, domain_edges);

    let dewpoint = read_raw_surface("2d", input_shape, data)?;
    let dewpoint = truncate_surface_to_extent(&dewpoint, domain_edges);

    let u_wind = read_raw_surface("10u", input_shape, data)?;
    let u_wind = truncate_surface_to_extent(&u_wind, domain_edges);

    let v_wind = read_raw_surface("10v", input_shape, data)?;
    let v_wind = truncate_surface_to_extent(&v_wind, domain_edges);

    Ok(Surfaces {
        lons: coords.0,
        lats: coords.1,
        temperature,
        dewpoint,
        pressure,
        height,
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

    // a bit of magic
    // data values in GRIB are a vec of values row-by-row (x-axis is in WE direction)
    // we want a Array2 of provided `shape` with x-axis in WE direction
    // but from_shape_vec(final_shape, data) splits the data into final_shape.1 long chunks
    // and puts them in columns
    // so we need to correctly split the data in GRIB vector into Array2 and then transpose
    // that array to get axes along expected geographical directions
    let result_data = Array2::from_shape_vec((shape.1, shape.0), data_level)?;
    let result_data = result_data.reversed_axes();
    let result_data = result_data.mapv(|v| v as Float);

    Ok(result_data)
}

/// Truncates surface data array from GRIB file to
/// cover only the domain + margins extent.
fn truncate_surface_to_extent(
    raw_field: &Array2<Float>,
    domain_edges: DomainExtent<usize>,
) -> Array2<Float> {
    //truncate in NS axis
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
