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

//! Module with methods for accessing the
//! environment and surface boundary
//! conditions data.

use super::{bisection, Environment, FieldTypes, SurfaceTypes};
use crate::{
    errors::{EnvironmentError, SearchError},
    model::environment::interpolation::{interpolate_bicubic, interpolate_tricubic},
    Float,
};
use ndarray::s;

impl Environment {
    /// Function to get interpolated value of given
    /// surface field at given (cartographic) coordinates.
    pub fn get_surface_value(
        &self,
        x: Float,
        y: Float,
        field: SurfaceTypes,
    ) -> Result<Float, EnvironmentError> {
        let (lon, lat) = self.projection.inverse_project(x, y);

        let west_lon_index = bisection::find_left_closest(
            self.surfaces.lons.slice(s![.., 0]).as_slice().unwrap(),
            &lon,
        )?;

        let south_lat_index = bisection::find_left_closest(
            self.surfaces
                .lats
                .slice(s![west_lon_index, ..])
                .as_slice()
                .unwrap(),
            &lat,
        )?;

        // check if indexes are not on the rightmost edge

        if west_lon_index == self.fields.lons.slice(s![.., 0]).len() - 1
            || south_lat_index == self.fields.lats.slice(s![west_lon_index, ..]).len() - 1
        {
            return Err(SearchError::OutOfBounds.into());
        }

        let field = match field {
            SurfaceTypes::Temperature => self.surfaces.temperature_coeffs.view(),
            SurfaceTypes::Dewpoint => self.surfaces.dewpoint_coeffs.view(),
            SurfaceTypes::Pressure => self.surfaces.pressure_coeffs.view(),
            SurfaceTypes::Height => self.surfaces.height_coeffs.view(),
            #[cfg(feature = "3d")]
            SurfaceTypes::UWind => self.surface.u_wind_coeffs.view(),
            #[cfg(feature = "3d")]
            SurfaceTypes::VWind => self.surface.v_wind_coeffs.view(),
        };

        let result_val = interpolate_bicubic(x, y, &field[[west_lon_index, south_lat_index]]);

        Ok(result_val)
    }

    /// Function to get interpolated value of given
    /// environment field at given (cartographic) coordinates.
    pub fn get_field_value(
        &self,
        x: Float,
        y: Float,
        z: Float,
        field: FieldTypes,
    ) -> Result<Float, EnvironmentError> {
        let (lon, lat) = self.projection.inverse_project(x, y);

        let west_lon_index = bisection::find_left_closest(
            self.fields.lons.slice(s![.., 0]).as_slice().unwrap(),
            &lon,
        )?;

        let south_lat_index = bisection::find_left_closest(
            self.fields
                .lats
                .slice(s![west_lon_index, ..])
                .as_slice()
                .unwrap(),
            &lat,
        )?;

        let z_index_search_array = self
            .fields
            .height
            .slice(s![.., west_lon_index, south_lat_index])
            .to_vec();

        let z_index = bisection::find_left_closest(&z_index_search_array, &z).or_else(|err| {
            // when searched height is below the lowest level
            // we set lowest point to 0-level for extrapolation
            // in all other cases error is returned

            match err {
                SearchError::OutOfBounds => {
                    if z < self.fields.height[[0, west_lon_index, south_lat_index]] {
                        Ok(0)
                    } else {
                        Err(err)
                    }
                }
                SearchError::EmptyArray => Err(err),
            }
        })?;

        // check if indexes are not on the rightmost edge

        if west_lon_index == self.fields.lons.slice(s![.., 0]).len() - 1
            || south_lat_index == self.fields.lats.slice(s![west_lon_index, ..]).len() - 1
            || z_index == z_index_search_array.len() - 1
        {
            return Err(SearchError::OutOfBounds.into());
        }

        let field = match field {
            FieldTypes::Pressure => self.fields.pressure_coeffs.view(),
            FieldTypes::VirtualTemperature => self.fields.virtual_temp_coeffs.view(),
            FieldTypes::UWind => self.fields.u_wind_coeffs.view(),
            FieldTypes::VWind => self.fields.v_wind_coeffs.view(),
        };

        let result_val =
            interpolate_tricubic(x, y, z, &field[[z_index, west_lon_index, south_lat_index]]);

        Ok(result_val)
    }
}
