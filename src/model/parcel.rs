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

//! (TODO: What it is)
//!
//! (Why it is neccessary)

use std::sync::Arc;

use chrono::NaiveDateTime;
use floccus::{mixing_ratio, virtual_temperature};

use crate::{
    errors::{ModelError, ParcelError},
    Float,
};

use super::{
    configuration::Config,
    environment::{
        Environment,
        SurfaceFields::{Dewpoint, Height, Pressure, Temperature, UWind, VWind},
    },
};

type Vec3<T> = (T, T, T);

struct ParcelState {
    datetime: NaiveDateTime,
    position: Vec3<Float>,
    velocity: Vec3<Float>,
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
) -> Result<(), ModelError> {
    // currently, parcel deployed directly from surface
    // but then (configurable) mixed parcel
    let initial_time = config.datetime.start;

    let x_pos = start_coords.0;
    let y_pos = start_coords.1;
    let z_pos = environment.get_surface_value(x_pos, y_pos, Height)?;

    // currently, constant initial vertical velocity (0.2 m/s)
    // but then lifiting can be taken into account
    // also as initial acceleration
    let x_vel = environment.get_surface_value(x_pos, y_pos, UWind)?;
    let y_vel = environment.get_surface_value(x_pos, y_pos, VWind)?;
    let z_vel = 0.2;

    let pres = environment.get_surface_value(x_pos, y_pos, Pressure)?;
    let temp = environment.get_surface_value(x_pos, y_pos, Temperature)?;
    let dwpt = environment.get_surface_value(x_pos, y_pos, Dewpoint)?;

    let mxng_rto = mixing_ratio::accuracy1(dwpt, pres)?;
    let satr_mxng_rto = mixing_ratio::accuracy1(temp, pres)?;
    let vrt_temp = virtual_temperature::general1(temp, mxng_rto)?;

    let initial_state = ParcelState {
        datetime: initial_time,
        position: (x_pos, y_pos, z_pos),
        velocity: (x_vel, y_vel, z_vel),
        pres,
        temp,
        mxng_rto,
        satr_mxng_rto,
        vrt_temp,
    };

    let parcel_result = run_simulation(initial_state, config.datetime.timestep);

    Ok(())
}

fn run_simulation(
    initial_state: ParcelState,
    timestep: Float,
) -> (Vec<ParcelState>, Result<(), ParcelError>) {
    (vec![], Ok(()))
}
