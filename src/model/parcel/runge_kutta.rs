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

mod schemes;

use super::{ParcelState, Vec3};
use floccus::constants::G;
use schemes::{AdiabaticScheme, PseudoAdiabaticScheme};
use crate::errors::ParcelSimulationError;
use crate::model::environment::EnvFields::{UWind, VWind, VirtualTemperature};
use crate::{model::environment::Environment, Float};
use chrono::Duration;
use log::debug;
use std::sync::Arc;

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Clone, Debug)]
pub(super) struct RungeKuttaDynamics<'a> {
    timestep: Float,
    env: &'a Arc<Environment>,
    pub parcel_log: Vec<ParcelState>,
}

impl<'a> RungeKuttaDynamics<'a> {
    pub fn new(
        initial_state: ParcelState,
        timestep: Float,
        environment: &'a Arc<Environment>,
    ) -> Self {
        let parcel_log = vec![initial_state];

        RungeKuttaDynamics {
            timestep,
            env: environment,
            parcel_log,
        }
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    pub fn run_simulation(&mut self) -> Result<(), ParcelSimulationError> {
        // from parcel theory: ascent adiabatic until saturation
        self.ascent_adiabatically()?;

        // from parcel theory: ascent pseudoadiabatic after saturation
        self.ascent_pseudoadiabatically()?;

        // for dry parcel pseudoadiabatic process is effectively adiabatic
        // so changing ascent for performance and accuracy
        self.ascent_adiabatically()?;

        Ok(())
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn ascent_adiabatically(&mut self) -> Result<(), ParcelSimulationError> {
        let initial_state = self.parcel_log.last().unwrap();

        if initial_state.velocity.z <= 0.0 {
            return Ok(());
        }

        debug!("Starting adiabatic ascent");
        debug!("Init state: {:?}", initial_state);

        let adiabatic_scheme = AdiabaticScheme::new(initial_state, self.env);

        loop {
            let ref_parcel = *self.parcel_log.last().unwrap();

            // holographic parcel is a virtual parcel that is moved
            // around for RK4 computations but doesn't change its
            // thermodynamic properties in reference to the prestep state
            let holo_parcel = ref_parcel;
            let c_0 = ref_parcel.velocity;
            let k_0 =
                self.calculate_bouyancy_force(&adiabatic_scheme.state_at_position(&holo_parcel)?)?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += 0.5 * self.timestep * c_0;
            let c_1 = ref_parcel.velocity + 0.5 * self.timestep * k_0;
            let k_1 =
                self.calculate_bouyancy_force(&adiabatic_scheme.state_at_position(&holo_parcel)?)?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += 0.5 * self.timestep * c_1;
            let c_2 = ref_parcel.velocity + 0.5 * self.timestep * k_1;
            let k_2 =
                self.calculate_bouyancy_force(&adiabatic_scheme.state_at_position(&holo_parcel)?)?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += self.timestep * c_2;
            let c_3 = ref_parcel.velocity + self.timestep * k_2;
            let k_3 =
                self.calculate_bouyancy_force(&adiabatic_scheme.state_at_position(&holo_parcel)?)?;

            let delta_pos = (self.timestep / 6.0) * (c_0 + 2.0 * c_1 + 2.0 * c_2 + c_3);
            let delta_vel = (self.timestep / 6.0) * (k_0 + 2.0 * k_1 + 2.0 * k_2 + k_3);

            let mut result_parcel = ref_parcel;
            result_parcel.datetime += Duration::milliseconds((self.timestep * 1000.0) as i64);
            result_parcel.position += delta_pos;
            result_parcel.velocity += delta_vel;

            if cfg!(feature = "3d") {
                result_parcel.velocity.x = self.env.get_field_value(
                    result_parcel.position.x,
                    result_parcel.position.y,
                    result_parcel.position.z,
                    UWind,
                )?;

                result_parcel.velocity.y = self.env.get_field_value(
                    result_parcel.position.x,
                    result_parcel.position.y,
                    result_parcel.position.z,
                    VWind,
                )?;
            }

            result_parcel = adiabatic_scheme.state_at_position(&result_parcel)?;

            if result_parcel.velocity.z <= 0.0
                || result_parcel.mxng_rto > result_parcel.satr_mxng_rto
            {
                break;
            }

            self.parcel_log.push(result_parcel);
        }

        Ok(())
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn ascent_pseudoadiabatically(&mut self) -> Result<(), ParcelSimulationError> {
        let initial_state = self.parcel_log.last().unwrap();

        if initial_state.velocity.z <= 0.0 || initial_state.mxng_rto < 0.000_001 {
            return Ok(());
        }

        debug!("Starting pseudoadiabatic ascent");
        debug!("Init state: {:?}", initial_state);

        let mut pseudoadiabatic_scheme = PseudoAdiabaticScheme::new(initial_state, self.env);

        loop {
            let ref_parcel = *self.parcel_log.last().unwrap();

            // holographic parcel is a virtual parcel that is moved
            // around for RK4 computations but doesn't change its
            // thermodynamic properties in reference to the prestep state
            let holo_parcel = ref_parcel;
            let c_0 = ref_parcel.velocity;
            let k_0 = self.calculate_bouyancy_force(
                &pseudoadiabatic_scheme.state_at_position(&holo_parcel)?,
            )?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += 0.5 * self.timestep * c_0;
            let c_1 = ref_parcel.velocity + 0.5 * self.timestep * k_0;
            let k_1 = self.calculate_bouyancy_force(
                &pseudoadiabatic_scheme.state_at_position(&holo_parcel)?,
            )?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += 0.5 * self.timestep * c_1;
            let c_2 = ref_parcel.velocity + 0.5 * self.timestep * k_1;
            let k_2 = self.calculate_bouyancy_force(
                &pseudoadiabatic_scheme.state_at_position(&holo_parcel)?,
            )?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += self.timestep * c_2;
            let c_3 = ref_parcel.velocity + self.timestep * k_2;
            let k_3 = self.calculate_bouyancy_force(
                &pseudoadiabatic_scheme.state_at_position(&holo_parcel)?,
            )?;

            let delta_pos = (self.timestep / 6.0) * (c_0 + 2.0 * c_1 + 2.0 * c_2 + c_3);
            let delta_vel = (self.timestep / 6.0) * (k_0 + 2.0 * k_1 + 2.0 * k_2 + k_3);

            let mut result_parcel = ref_parcel;
            result_parcel.datetime += Duration::milliseconds((self.timestep * 1000.0) as i64);
            result_parcel.position += delta_pos;
            result_parcel.velocity += delta_vel;

            if cfg!(feature = "3d") {
                result_parcel.velocity.x = self.env.get_field_value(
                    result_parcel.position.x,
                    result_parcel.position.y,
                    result_parcel.position.z,
                    UWind,
                )?;

                result_parcel.velocity.y = self.env.get_field_value(
                    result_parcel.position.x,
                    result_parcel.position.y,
                    result_parcel.position.z,
                    VWind,
                )?;
            }

            result_parcel = pseudoadiabatic_scheme.state_at_position(&result_parcel)?;

            if result_parcel.velocity.z <= 0.0 || result_parcel.mxng_rto < 0.000_001 {
                break;
            }

            pseudoadiabatic_scheme.update_ref_state(&result_parcel);
            self.parcel_log.push(result_parcel);
        }

        Ok(())
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn calculate_bouyancy_force(
        &self,
        parcel: &ParcelState,
    ) -> Result<Vec3, ParcelSimulationError> {
        let tv_env = self.env.get_field_value(
            parcel.position.x,
            parcel.position.y,
            parcel.position.z,
            VirtualTemperature,
        )?;
        let bouyancy_force = G * ((parcel.vrt_temp - tv_env) / tv_env);

        Ok(Vec3 {
            x: 0.0,
            y: 0.0,
            z: bouyancy_force,
        })
    }
}
