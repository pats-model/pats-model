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

//! Module responsible for reading and storing boundary conditions
//! (environment) data, and providing that data to parcels.

mod accesser;
mod bisection;
mod buffer;
mod interpolation;
mod projection;

use super::configuration::{self, Config, Domain};
use crate::constants::{NS_C_EARTH, WE_C_EARTH};
use crate::errors::InputError;
use crate::model::environment::projection::LambertConicConformal;
use crate::{errors::EnvironmentError, Float};
use eccodes::codes_handle::{CodesHandle, KeyType, ProductKind::GRIB};
use eccodes::FallibleIterator;
use log::debug;
use ndarray::{Array2, Array3};
use rustc_hash::FxHashSet;

/// Convenience type to store lat-lon coordinates.
type LonLat<T> = (T, T);

/// Enum containing fields on pressure
/// levels that can be requested.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum EnvFields {
    Pressure,
    VirtualTemperature,
    UWind,
    VWind,
}

/// Struct for storing environmental variables
/// from levels above ground (currently pressure levels).
///
/// To limit IO operations and reduce performance overhead
/// of the model boundary conditions data is stored in the
/// memory as 3D arrays.
#[derive(Debug)]
pub struct Fields {
    temperature: Array3<Float>,
    pressure: Array3<Float>,
    height: Array3<Float>,
    u_wind: Array3<Float>,
    v_wind: Array3<Float>,
    spec_humidity: Array3<Float>,
    virtual_temp: Array3<Float>,
    lons: Array2<Float>,
    lats: Array2<Float>,
}

impl Fields {
    fn new_empty() -> Self {
        Fields {
            temperature: Array3::default((0, 0, 0)),
            pressure: Array3::default((0, 0, 0)),
            height: Array3::default((0, 0, 0)),
            u_wind: Array3::default((0, 0, 0)),
            v_wind: Array3::default((0, 0, 0)),
            spec_humidity: Array3::default((0, 0, 0)),
            virtual_temp: Array3::default((0, 0, 0)),
            lons: Array2::default((0, 0)),
            lats: Array2::default((0, 0)),
        }
    }
}

/// Enum containing surface fields
/// that can be requested.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum SurfaceFields {
    Temperature,
    Dewpoint,
    Pressure,
    Height,
    #[cfg(feature="3d")]
    UWind,
    #[cfg(feature="3d")]
    VWind,
}

/// Struct for storing environmental variables at/near surface.
///
/// To limit IO operations and reduce performance overhead
/// of the model surface data is stored in the
/// memory as 2D arrays.
#[derive(Debug)]
pub struct Surface {
    temperature: Array2<Float>,
    dewpoint: Array2<Float>,
    pressure: Array2<Float>,
    height: Array2<Float>,
    u_wind: Array2<Float>,
    v_wind: Array2<Float>,
    lons: Array2<Float>,
    lats: Array2<Float>,
}

impl Surface {
    fn new_empty() -> Self {
        Surface {
            temperature: Array2::default((0, 0)),
            dewpoint: Array2::default((0, 0)),
            pressure: Array2::default((0, 0)),
            height: Array2::default((0, 0)),
            u_wind: Array2::default((0, 0)),
            v_wind: Array2::default((0, 0)),
            lons: Array2::default((0, 0)),
            lats: Array2::default((0, 0)),
        }
    }
}

/// Environment main struct storing and providing
/// boundary condition (environment) data.
///
/// Use of a separate struct for handling boundary
/// conditions allows to have a clean API for
/// other model parts to use.
#[derive(Debug)]
pub struct Environment {
    fields: Fields,
    surface: Surface,
    input: configuration::Input,
    pub projection: LambertConicConformal,
    levels: Vec<i64>,
}

impl Environment {
    /// Environment struct constructor
    /// responsible for reading GRIB files
    /// and buffering data in domain extent.
    pub fn new(config: &Config) -> Result<Self, EnvironmentError> {
        debug!("Creating new enviroment");

        let fields = Fields::new_empty();
        let surface = Surface::new_empty();

        let projection = generate_domain_projection(&config.domain)?;

        let (west_south, east_north) = compute_domain_extent(&config.domain, &projection);

        let levels = list_levels(&config.input)?;

        let mut new_env = Environment {
            fields,
            surface,
            input: config.input.clone(),
            projection,
            levels,
        };

        new_env.buffer_fields(west_south, east_north)?;
        new_env.buffer_surface(west_south, east_north)?;

        Ok(new_env)
    }
}

/// Function to create a geographic projection struct
/// with parameters that allow for lowest distorion
/// for a given domain.
fn generate_domain_projection(domain: &Domain) -> Result<LambertConicConformal, EnvironmentError> {
    let sides = measure_domain_sides(domain);

    // if there's only one parcel to release in some direction
    // we use the special case for projection
    let lon_0;
    let lat_1;
    let lat_2;

    if sides.0 < 0.1 {
        lon_0 = domain.ref_lon;
    } else {
        lon_0 = approx_central_lon(domain.ref_lon, domain.ref_lat, sides.0);
    }

    if sides.1 < 0.1 {
        lat_1 = domain.ref_lat - 1.0;
        lat_2 = domain.ref_lat + 1.0;
    } else {
        lat_1 = domain.ref_lat;
        lat_2 = compute_top_lat(domain.ref_lat, sides.1);
    }

    let projection = LambertConicConformal::new(lon_0, lat_1, lat_2)?;

    Ok(projection)
}

/// Function to get domain sides length
/// in meters.
fn measure_domain_sides(domain: &Domain) -> (Float, Float) {
    let x_side = (domain.shape.0 - 1) as Float * domain.spacing;
    let y_side = (domain.shape.1 - 1) as Float * domain.spacing;

    (x_side, y_side)
}

/// Function to compute the latitude of domain top
/// on the WGS84 ellipsoid.
fn compute_top_lat(lat: Float, distance: Float) -> Float {
    let degree_length = NS_C_EARTH / 360.0;
    let arc_distance = distance / degree_length;

    lat + arc_distance
}

/// Function to approximate the longitude of domain centre
/// on the WGS84 ellipsoid.
fn approx_central_lon(lon_0: Float, lat_0: Float, distance: Float) -> Float {
    let degree_length = lat_0.to_radians().cos() * (WE_C_EARTH / 360.0);
    let half_arc_length = (distance / 2.0) / degree_length;

    lon_0 + half_arc_length
}

/// Function to get a lat-lon extent of domain with margins.
fn compute_domain_extent(
    domain: &Domain,
    projection: &LambertConicConformal,
) -> (LonLat<Float>, LonLat<Float>) {
    let sw_xy = projection.project(domain.ref_lon, domain.ref_lat);

    let ne_xy = (
        sw_xy.0 + ((domain.shape.0 - 1) as Float * domain.spacing),
        sw_xy.1 + ((domain.shape.1 - 1) as Float * domain.spacing),
    );

    let ne_lonlat = projection.inverse_project(ne_xy.0, ne_xy.1);

    (
        (
            domain.ref_lon - domain.margins.0,
            domain.ref_lat - domain.margins.1,
        ),
        (
            ne_lonlat.0 + domain.margins.0,
            ne_lonlat.1 + domain.margins.1,
        ),
    )
}

/// Function to get the list of unique levels
/// of specified type in the provided GRIB files.
fn list_levels(input: &configuration::Input) -> Result<Vec<i64>, InputError> {
    debug!("Getting levels list");

    let mut unique_levels: FxHashSet<i64> = FxHashSet::default();

    for file_name in &input.data_files {
        let mut handle = CodesHandle::new_from_file(file_name, GRIB)?;

        while let Some(level) = &mut handle.next()? {
            if level.read_key("typeOfLevel")?.value == KeyType::Str(input.level_type.clone()) {
                let level_id = level.read_key("level")?;

                if let KeyType::Int(id) = level_id.value {
                    unique_levels.insert(id);
                } else {
                    return Err(InputError::IncorrectKeyType("level"));
                };
            }
        }
    }

    let mut unique_levels: Vec<i64> = unique_levels.into_iter().collect();
    unique_levels.sort_unstable();
    unique_levels.reverse();

    Ok(unique_levels)
}
