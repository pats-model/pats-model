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

//! Module responsible for reading and storing boundary conditions
//! (environment) data, and providing that data to parcels.

mod accesser;
mod bisection;
mod fields;
mod interpolation;
mod projection;
mod surfaces;

use self::fields::Fields;
use self::surfaces::Surfaces;
use super::configuration::{Config, Domain};
use crate::constants::{NS_C_EARTH, WE_C_EARTH};
use crate::model::environment::projection::LambertConicConformal;
use crate::{errors::EnvironmentError, Float};
use log::debug;

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
struct DomainExtent<T> {
    north: T,
    south: T,
    west: T,
    east: T,
}

/// Enum containing fields on pressure
/// levels that can be requested.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum EnvFields {
    Pressure,
    Temperature,
    VirtualTemperature,
    UWind,
    VWind,
    VerticalVel,
}

/// Enum containing surface fields
/// that can be requested.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum SurfaceFields {
    Temperature,
    Dewpoint,
    Pressure,
    Height,
    #[cfg(feature = "3d")]
    UWind,
    #[cfg(feature = "3d")]
    VWind,
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
    surfaces: Surfaces,
    pub projection: LambertConicConformal,
}

impl Environment {
    /// Environment struct constructor
    /// responsible for reading GRIB files
    /// and buffering data in domain extent.
    pub fn new(config: &Config) -> Result<Self, EnvironmentError> {
        debug!("Creating new enviroment");

        let projection = generate_domain_projection(&config.domain)?;
        let domain_edges = compute_domain_edges(config, &projection);

        let fields = Fields::new(&config.input, domain_edges)?;
        let surfaces = Surfaces::new(&config.input, domain_edges)?;

        Ok(Environment {
            fields,
            surfaces,
            projection,
        })
    }
}

/// Function to create a geographic projection struct
/// with parameters that allow for lowest distorion
/// for a given domain.
fn generate_domain_projection(domain: &Domain) -> Result<LambertConicConformal, EnvironmentError> {
    let sides = measure_domain_sides(domain);

    // if there's only one parcel to release in some direction
    // we use the special case for projection
    let lat_1;
    let lat_2;

    let lon_0 = if sides.0 < 0.1 {
        domain.ref_lon
    } else {
        approx_central_lon(domain.ref_lon, domain.ref_lat, sides.0)
    };

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
    let x_side = Float::from(domain.shape.0 - 1) * domain.spacing;
    let y_side = Float::from(domain.shape.1 - 1) * domain.spacing;

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
fn compute_domain_edges(
    config: &Config,
    projection: &LambertConicConformal,
) -> DomainExtent<usize> {
    let sw_xy = projection.project(config.domain.ref_lon, config.domain.ref_lat);

    let ne_xy = (
        sw_xy.0 + (Float::from(config.domain.shape.0 - 1) * config.domain.spacing),
        sw_xy.1 + (Float::from(config.domain.shape.1 - 1) * config.domain.spacing),
    );

    let ne_lonlat = projection.inverse_project(ne_xy.0, ne_xy.1);

    let domain_extent = DomainExtent {
        west: config.domain.ref_lon - config.domain.margins.0,
        south: config.domain.ref_lat - config.domain.margins.1,
        east: ne_lonlat.0 + config.domain.margins.0,
        north: ne_lonlat.1 + config.domain.margins.1,
    };

    debug!(
        "Computed buffering extent: N{:.2}<->{:.2} E{:.2}<->{:.2}",
        domain_extent.south, domain_extent.north, domain_extent.west, domain_extent.east
    );

    let distinct_lonlats = &config.input.distinct_lonlats;
    find_extent_edge_indices(distinct_lonlats, domain_extent)
}

/// Finds closests indices in the GRIB input files
/// grid that fully cover domain with margins (it is
/// with some excess).
fn find_extent_edge_indices(
    distinct_lonlats: &(Vec<Float>, Vec<Float>),
    domain_extent: DomainExtent<Float>,
) -> DomainExtent<usize> {
    let edge_lats = (
        bisection::find_left_closest(&distinct_lonlats.1, &domain_extent.north).unwrap(),
        bisection::find_right_closest(&distinct_lonlats.1, &domain_extent.south).unwrap(),
    );
    let edge_lons = (
        bisection::find_left_closest(
            &distinct_lonlats.0,
            &convert_to_grib_longitudes(domain_extent.west),
        )
        .unwrap(),
        bisection::find_right_closest(
            &distinct_lonlats.0,
            &convert_to_grib_longitudes(domain_extent.east),
        )
        .unwrap(),
    );

    DomainExtent {
        north: edge_lats.0,
        south: edge_lats.1,
        west: edge_lons.0,
        east: edge_lons.1,
    }
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
