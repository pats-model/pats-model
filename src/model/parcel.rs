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

mod runge_kutta;

use log::error;
use runge_kutta::RungeKuttaDynamics;
use std::{
    ops::{Add, AddAssign, Mul},
    sync::Arc,
};

use chrono::NaiveDateTime;
use floccus::{mixing_ratio, virtual_temperature};

use crate::{errors::ModelError, Float};

use super::{
    configuration::Config,
    environment::{
        Environment,
        SurfaceFields::{Dewpoint, Height, Pressure, Temperature, UWind, VWind},
    },
};

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
struct Vec3 {
    x: Float,
    y: Float,
    z: Float,
}

impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Self::Output {
        Vec3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        *self = Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        };
    }
}

impl Mul<Vec3> for Float {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Self::Output {
        Vec3 {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

impl Mul<Float> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: Float) -> Self::Output {
        Vec3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

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
) -> Result<(), ModelError> {
    let initial_state = prepare_parcel(start_coords, config, environment)?;

    let mut dynamic_scheme =
        RungeKuttaDynamics::new(initial_state, config.datetime.timestep, environment);

    let parcel_result = dynamic_scheme.run_simulation();

    // if the parcel simulation stops with error
    // we report that in log, but do not return the error
    // as parcel has been deployed
    if let Err(err) = parcel_result {
        error!("Parcel released from x: {:.2} y: {:.2} has stopped its ascent with error: {}. Check your configuration.", 
        start_coords.0, start_coords.1, err);
        return Ok(());
    }

    

    Ok(())
}

fn prepare_parcel(
    start_coords: (Float, Float),
    config: &Arc<Config>,
    environment: &Arc<Environment>,
) -> Result<ParcelState, ModelError> {
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

fn save_parcel_log(parcel_log: Vec<ParcelState>) {

}
