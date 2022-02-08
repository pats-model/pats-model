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

//! (TODO: What it is)
//!
//! (Why it is neccessary)

use super::ParcelState;
use crate::{
    errors::ParcelError,
    model::environment::{Environment, FieldTypes::VirtualTemperature},
    Float,
};
use float_cmp::approx_eq;
use floccus::constants::G;
use serde::Serialize;
use std::sync::Arc;

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default, Serialize)]
pub struct ConvectiveParams {
    start_lon: Float,
    start_lat: Float,

    // Parcel Top Height
    parcel_top: Float,

    // Parcel displacement from initial point
    x_displac: Float,
    y_displac: Float,

    // Parcel Maximum Vertical Velocity
    max_vert_vel: Float,

    // Condensation Level
    // (similar to Convective Condensation Level)
    condens_lvl: Option<Float>,

    // Level of Free Convection
    lfc: Option<Float>,

    // Equilibrium Level
    el: Option<Float>,

    // Convective Available Potential Energy
    cape: Option<Float>,

    // Convective Inhibition
    cin: Option<Float>,
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub(super) fn compute_conv_params(
    parcel_log: &[ParcelState],
    environment: &Arc<Environment>,
) -> Result<ConvectiveParams, ParcelError> {
    let mut result_params = ConvectiveParams::default();

    // add parcel identification
    let parcel_start_coords = environment.projection.inverse_project(
        parcel_log.first().unwrap().position.x,
        parcel_log.first().unwrap().position.y,
    );

    result_params.start_lon = parcel_start_coords.0;
    result_params.start_lat = parcel_start_coords.1;

    // get environmental virtual temperature along parcel trace
    // to avoid calls to Environment
    let env_vrt_tmp = get_env_vtemp(parcel_log, environment)?;

    result_params.update_displacements(parcel_log);
    result_params.update_levels(parcel_log, &env_vrt_tmp);
    result_params.update_thermodynamic_vars(parcel_log, &env_vrt_tmp);

    Ok(result_params)
}

impl ConvectiveParams {
    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn update_displacements(&mut self, parcel_log: &[ParcelState]) {
        self.parcel_top = parcel_log.last().unwrap().position.z;

        self.x_displac =
            parcel_log.last().unwrap().position.x - parcel_log.first().unwrap().position.x;
        self.y_displac =
            parcel_log.last().unwrap().position.y - parcel_log.first().unwrap().position.y;

        self.max_vert_vel = parcel_log
            .iter()
            .max_by(|x, y| {
                x.velocity
                    .z
                    .partial_cmp(&y.velocity.z)
                    .expect("Float comparison failed")
            })
            .expect("Parcel log is empty")
            .velocity
            .z;
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn update_levels(&mut self, parcel_log: &[ParcelState], env_vrt_tmp: &[Float]) {
        // searched levels are subsequent and interdependent, so we look for them in loops
        // iterating from log beginning so from ascent bottom
        let mut ccl_index = 0;

        for (i, point) in parcel_log.iter().enumerate() {
            // first time this is true is condensation level
            if point.mxng_rto >= point.satr_mxng_rto {
                self.condens_lvl = Some(point.position.z);
                ccl_index = i;
                break;
            }
        }

        let mut lfc_index = 0;

        if self.condens_lvl.is_some() {
            // we check the condensation level as it might be a level of free convection
            for i in ccl_index..parcel_log.len() {
                let point = parcel_log[i];

                // first time this is true is LFC
                if point.vrt_temp > env_vrt_tmp[i] {
                    self.lfc = Some(point.position.z);
                    lfc_index = i;
                    break;
                }
            }
        }

        if self.lfc.is_some() {
            // start checking from level after LFC for rare case when virtual temperatures are equal
            for i in (lfc_index + 1)..parcel_log.len() {
                let point = parcel_log[i];

                // first time this is true is LFC
                if point.vrt_temp <= env_vrt_tmp[i] {
                    self.el = Some(point.position.z);
                    break;
                }
            }
        }
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn update_thermodynamic_vars(&mut self, parcel_log: &[ParcelState], env_vrt_tmp: &[Float]) {
        let mut lfc_id = 0;

        // compute CIN if LFC is present
        let mut cin: Float = 0.0;
        if self.lfc.is_some() {
            // we start from the 2nd point of parcel log to not go out of bounds
            for i in 1..parcel_log.len() {
                let point = parcel_log[i];

                let y_1 = (point.vrt_temp - env_vrt_tmp[i]) / env_vrt_tmp[i];
                let y_0 = (parcel_log[i - 1].vrt_temp - env_vrt_tmp[i - 1]) / env_vrt_tmp[i - 1];

                let delta_z = point.position.z - parcel_log[i - 1].position.z;

                cin += ((y_0 + y_1) / 2.0) * delta_z;

                if approx_eq!(Float, point.position.z, self.lfc.unwrap()) {
                    lfc_id = i;
                    break;
                }
            }
        }

        self.cin = Some(-G * cin);

        // compute CAPE if LFC and EL is present
        let mut cape: Float = 0.0;
        if self.lfc.is_some() && self.el.is_some() {
            // we start integration from LFC
            for i in (lfc_id + 1)..parcel_log.len() {
                let point = parcel_log[i];

                let y_1 = (point.vrt_temp - env_vrt_tmp[i]) / env_vrt_tmp[i];
                let y_0 = (parcel_log[i - 1].vrt_temp - env_vrt_tmp[i - 1]) / env_vrt_tmp[i - 1];

                let delta_z = point.position.z - parcel_log[i - 1].position.z;

                cape += ((y_0 + y_1) / 2.0) * delta_z;

                if approx_eq!(Float, point.position.z, self.el.unwrap()) {
                    break;
                }
            }
        }

        self.cape = Some(G * cape);
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn get_env_vtemp(
    parcel_log: &[ParcelState],
    environment: &Arc<Environment>,
) -> Result<Vec<Float>, ParcelError> {
    let env_vtemp: Result<Vec<_>, _> = parcel_log
        .iter()
        .map(|pst| {
            environment.get_field_value(
                pst.position.x,
                pst.position.y,
                pst.position.z,
                VirtualTemperature,
            )
        })
        .collect();

    Ok(env_vtemp?)
}
