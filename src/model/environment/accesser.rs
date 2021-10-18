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

//! Module with methods for accessing the
//! environment and surface boundary
//! conditions data.

use ndarray::s;

use crate::{
    errors::{EnvironmentError, SearchError},
    model::environment::interpolation::{
        interpolate_bilinear, interpolate_tilinear, Point2D, Point3D,
    },
    Float,
};

use super::{bisection, EnvFields, Environment, SurfaceFields};

impl Environment {
    /// Function to get interpolated value of given
    /// surface field at given (cartographic) coordinates.
    pub fn get_surface_value(
        &self,
        x: Float,
        y: Float,
        field: SurfaceFields,
    ) -> Result<Float, EnvironmentError> {
        let (lon, lat) = self.projection.inverse_project(x, y);

        let west_lon_index = bisection::find_left_closest(
            self.surface.lons.slice(s![.., 0]).as_slice().unwrap(),
            &lon,
        )?;

        let south_lat_index = bisection::find_left_closest(
            self.surface
                .lats
                .slice(s![west_lon_index, ..])
                .as_slice()
                .unwrap(),
            &lat,
        )?;

        let field = match field {
            SurfaceFields::Temperature => self.surface.temperature.view(),
            SurfaceFields::Dewpoint => self.surface.dewpoint.view(),
            SurfaceFields::Pressure => self.surface.pressure.view(),
            SurfaceFields::Height => self.surface.height.view(),
            SurfaceFields::UWind => self.surface.u_wind.view(),
            SurfaceFields::VWind => self.surface.v_wind.view(),
        };

        let horizontal_points = [
            (west_lon_index, south_lat_index),
            (west_lon_index, south_lat_index + 1),
            (west_lon_index + 1, south_lat_index),
            (west_lon_index + 1, south_lat_index + 1),
        ];

        let mut ref_points = [Point2D::default(); 4];

        for (i, (x_index, y_index)) in horizontal_points.iter().enumerate() {
            let (lon, lat) = (
                self.fields.lons[[*x_index, *y_index]],
                self.fields.lats[[*x_index, *y_index]],
            );
            let (x, y) = self.projection.project(lon, lat);

            ref_points[i] = Point2D {
                x,
                y,
                value: field[[*x_index, *y_index]],
            };
        }

        let result_val = interpolate_bilinear(x, y, ref_points);

        Ok(result_val)
    }

    /// Function to get interpolated value of given
    /// environment field at given (cartographic) coordinates.
    pub fn get_field_value(
        &self,
        x: Float,
        y: Float,
        z: Float,
        field: EnvFields,
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

        let field = match field {
            EnvFields::Temperature => self.fields.temperature.view(),
            EnvFields::Pressure => self.fields.pressure.view(),
            EnvFields::VirtualTemperature => self.fields.virtual_temp.view(),
            EnvFields::SpecificHumidity => self.fields.spec_humidity.view(),
            EnvFields::UWind => self.fields.u_wind.view(),
            EnvFields::VWind => self.fields.v_wind.view(),
        };

        let horizontal_points = [
            (west_lon_index, south_lat_index),
            (west_lon_index, south_lat_index + 1),
            (west_lon_index + 1, south_lat_index),
            (west_lon_index + 1, south_lat_index + 1),
        ];

        let mut ref_points = [Point3D::default(); 8];

        for (i, (x_index, y_index)) in horizontal_points.iter().enumerate() {
            let z_index_search_array = self
                .fields
                .height
                .slice(s![.., *x_index, *y_index])
                .to_vec();

            let z_index =
                bisection::find_left_closest(&z_index_search_array, &z).or_else(|err| {
                    // when searched height is below the lowest level
                    // we set lowest point to 0-level for extrapolation
                    // in all other cases error is returned

                    match err {
                        SearchError::OutOfBounds => {}
                        _ => {
                            return Err(err);
                        }
                    }

                    // height decreases with height
                    if &z > &self.fields.height[[0, *x_index, *y_index]] {
                        return Ok(0);
                    }

                    return Err(err);
                })?;

            let (lon, lat) = (
                self.fields.lons[[*x_index, *y_index]],
                self.fields.lats[[*x_index, *y_index]],
            );
            let (x, y) = self.projection.project(lon, lat);

            // bottom point
            ref_points[i] = Point3D {
                x,
                y,
                z: self.fields.height[[z_index, *x_index, *y_index]],
                value: field[[z_index, *x_index, *y_index]],
            };

            // upper point
            ref_points[i + 4] = Point3D {
                x,
                y,
                z: self.fields.height[[z_index + 1, *x_index, *y_index]],
                value: field[[z_index + 1, *x_index, *y_index]],
            };
        }

        let result_val = interpolate_tilinear(x, y, z, ref_points);

        Ok(result_val)
    }
}
