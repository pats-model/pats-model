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

//! (TODO: What it is)
//!
//! (Why it is neccessary)

pub(super) mod conv_params;
mod logger;
mod runge_kutta;

use self::conv_params::ConvectiveParams;
use super::{
    configuration::Config,
    environment::{
        Environment,
        SurfaceFields::{Dewpoint, Height, Pressure, Temperature},
    },
    vec3::Vec3,
};
use crate::{errors::ParcelError, model::parcel::conv_params::compute_conv_params, Float};
use chrono::NaiveDateTime;
use floccus::{mixing_ratio, virtual_temperature};
use log::debug;
use runge_kutta::RungeKuttaDynamics;
use std::sync::Arc;

#[cfg(feature = "3d")]
use super::environment::SurfaceFields::{UWind, VWind};

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
struct ParcelState {
    datetime: NaiveDateTime,
    position: Vec3,
    velocity: Vec3,
    pres: Float,
    temp: Float,
    mxng_rto: Float,
    satr_mxng_rto: Float,
    vrt_temp: Float,
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub fn deploy(
    start_coords: (Float, Float),
    config: &Arc<Config>,
    environment: &Arc<Environment>,
) -> Result<ConvectiveParams, ParcelError> {
    let initial_state = prepare_parcel(start_coords, config, environment)?;

    let mut dynamic_scheme =
        RungeKuttaDynamics::new(initial_state, config.datetime.timestep, environment);

    let parcel_result = dynamic_scheme.run_simulation();

    // if the parcel simulation stops with error
    // we report compute parcel's initial geographic
    // coords and return the error with that additional info
    if let Err(err) = parcel_result {
        let (lon, lat) = environment
            .projection
            .inverse_project(start_coords.0, start_coords.1);

        return Err(ParcelError::AscentStopped(lat, lon, err));
    }

    if cfg!(feature = "raw_output") {
        logger::save_parcel_log(&dynamic_scheme.parcel_log, environment)?;
    }

    let parcel_params = compute_conv_params(&dynamic_scheme.parcel_log, environment)?;

    Ok(parcel_params)
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn prepare_parcel(
    start_coords: (Float, Float),
    config: &Arc<Config>,
    environment: &Arc<Environment>,
) -> Result<ParcelState, ParcelError> {
    debug!("Preparing parcel at: {:?}", start_coords);
    // currently, parcel deployed directly from surface
    // but then (configurable) mixed parcel
    let initial_time = config.datetime.start;

    let x_pos = start_coords.0;
    let y_pos = start_coords.1;
    let z_pos = environment.get_surface_value(x_pos, y_pos, Height)?;

    #[cfg(feature = "3d")]
    let x_vel = environment.get_surface_value(x_pos, y_pos, UWind)?;
    #[cfg(feature = "3d")]
    let y_vel = environment.get_surface_value(x_pos, y_pos, VWind)?;

    #[cfg(not(feature = "3d"))]
    let x_vel = 0.0;
    #[cfg(not(feature = "3d"))]
    let y_vel = 0.0;

    // currently, constant initial vertical velocity (0.2 m/s)
    // but then lifiting can be taken into account
    // also as initial acceleration
    let z_vel = 0.2;

    let pres = environment.get_surface_value(x_pos, y_pos, Pressure)?;
    let temp = environment.get_surface_value(x_pos, y_pos, Temperature)?;
    let dwpt = environment.get_surface_value(x_pos, y_pos, Dewpoint)?;

    let mxng_rto = mixing_ratio::accuracy1(dwpt, pres)?;
    let satr_mxng_rto = mixing_ratio::accuracy1(temp, pres)?;
    let vrt_temp = virtual_temperature::general1(temp, mxng_rto)?;

    Ok(ParcelState {
        datetime: initial_time,
        position: Vec3 {
            x: x_pos,
            y: y_pos,
            z: z_pos,
        },
        velocity: Vec3 {
            x: x_vel,
            y: y_vel,
            z: z_vel,
        },
        pres,
        temp,
        mxng_rto,
        satr_mxng_rto,
        vrt_temp,
    })
}
