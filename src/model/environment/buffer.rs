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

//! Module containing methods responsible for
//! buffering data from GRIB input.

mod fields;
mod surface;

use crate::{
    errors::InputError,
    model::environment::{bisection, Environment, LonLat},
    Float,
};
use eccodes::{
    codes_handle::{
        CodesHandle,
        KeyType::{FloatArray, Int},
        ProductKind::GRIB,
    },
    FallibleIterator,
};

impl Environment {
    /// Function to read distinct longitudes and latitudes
    /// in input GRIB files.
    fn read_distinct_latlons(&self) -> Result<LonLat<Vec<Float>>, InputError> {
        // We can read any message from any file as we assume that lat-lons
        // are aligned in all GRIB messages
        let any_file = &self.input.data_files[0];
        let mut any_file = CodesHandle::new_from_file(any_file, GRIB)?;

        let any_message = any_file.next()?.ok_or(InputError::DataNotSufficient(
            "One or more input files does not contain any valid GRIB message",
        ))?;

        let mut distinct_latitudes: Vec<Float>;
        let mut distinct_longitudes: Vec<Float>;

        if let FloatArray(lats) = any_message.read_key("distinctLatitudes")?.value {
            distinct_latitudes = lats.into_iter().map(|v| v as Float).collect();
        } else {
            return Err(InputError::IncorrectKeyType("distinctLatitudes"));
        }

        if let FloatArray(lons) = any_message.read_key("distinctLongitudes")?.value {
            distinct_longitudes = lons.into_iter().map(|v| v as Float).collect();
        } else {
            return Err(InputError::IncorrectKeyType("distinctLongitudes"));
        }

        distinct_latitudes
            .sort_by(|a, b| a.partial_cmp(b).expect("Sorting distinct latitudes failed"));
        distinct_longitudes.sort_by(|a, b| {
            a.partial_cmp(b)
                .expect("Sorting distinct longitudes failed")
        });

        Ok((distinct_longitudes, distinct_latitudes))
    }

    /// Function to read a grid shape of input GRIB files.
    fn read_input_shape(&self) -> Result<(usize, usize), InputError> {
        let any_file = &self.input.data_files[0];
        let mut any_file = CodesHandle::new_from_file(any_file, GRIB)?;

        let any_message = any_file.next()?.ok_or(InputError::DataNotSufficient(
            "One or more input files does not contain any valid GRIB message",
        ))?;

        let ni;
        let nj;

        if let Int(val) = any_message.read_key("Ni")?.value {
            ni = val as usize;
        } else {
            return Err(InputError::IncorrectKeyType("Ni"));
        }

        if let Int(val) = any_message.read_key("Nj")?.value {
            nj = val as usize;
        } else {
            return Err(InputError::IncorrectKeyType("Nj"));
        }

        Ok((ni, nj))
    }
}

/// Finds closests indices in the GRIB input files
/// grid that fully cover domain with margins (it is
/// with some excess).
fn find_extent_edge_indices(
    distinct_lonlats: &(Vec<Float>, Vec<Float>),
    west_south: (Float, Float),
    east_north: (Float, Float),
) -> ((usize, usize), (usize, usize)) {
    let edge_lats = (
        bisection::find_left_closest(&distinct_lonlats.1, &west_south.1).unwrap(),
        bisection::find_right_closest(&distinct_lonlats.1, &east_north.1).unwrap(),
    );
    let edge_lons = (
        bisection::find_left_closest(
            &distinct_lonlats.0,
            &convert_to_grib_longitudes(west_south.0),
        )
        .unwrap(),
        bisection::find_right_closest(
            &distinct_lonlats.0,
            &convert_to_grib_longitudes(east_north.0),
        )
        .unwrap(),
    );
    (edge_lats, edge_lons)
}

/// Converts the longitude in convention used by model
/// (longitude between -180 and 180) to longitude 
/// in GRIB convention (any positive integer).
fn convert_to_grib_longitudes(longitude: Float) -> Float {
    if longitude < 0.0 {
        return 360.0 + longitude;
    }

    longitude
}
