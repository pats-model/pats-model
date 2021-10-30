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

use super::ParcelState;
use crate::{
    errors::ParcelError,
    model::environment::{EnvFields::VirtualTemperature, Environment},
    Float,
};
use floccus::constants::G;
use serde::Serialize;
use std::sync::Arc;

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default, Serialize)]
pub struct ConvectiveParams {
    // Parcel Top Height
    parcel_top: Float,

    // Parcel displacement from initial point
    horztl_displac: (Float, Float),

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

pub(super) fn compute_conv_params(
    parcel_log: &Vec<ParcelState>,
    environment: &Arc<Environment>,
) -> Result<ConvectiveParams, ParcelError> {
    let mut result_params = ConvectiveParams::default();

    let env_vrt_tmp = get_env_vtemp(parcel_log, environment)?;

    result_params.update_displacements(parcel_log);
    result_params.update_levels(parcel_log, &env_vrt_tmp)?;
    result_params.update_thermodynamic_vars(parcel_log, &env_vrt_tmp);

    Ok(result_params)
}

impl ConvectiveParams {
    fn update_displacements(&mut self, parcel_log: &Vec<ParcelState>) {
        self.parcel_top = parcel_log.last().unwrap().position.z;

        self.horztl_displac = (
            parcel_log.last().unwrap().position.x - parcel_log.first().unwrap().position.x,
            parcel_log.last().unwrap().position.y - parcel_log.first().unwrap().position.y,
        );

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

    fn update_levels(
        &mut self,
        parcel_log: &Vec<ParcelState>,
        env_vrt_tmp: &Vec<Float>,
    ) -> Result<(), ParcelError> {
        // searched levels are subsequent and interdependent, so we look for them in loops
        // iterating from log beginning so from ascent bottom
        let mut current_point = 0;

        for (i, point) in parcel_log.iter().enumerate() {
            // first time this is true is condensation level
            if point.mxng_rto >= point.satr_mxng_rto {
                self.condens_lvl = Some(point.position.z);
                current_point = i;
                break;
            }
        }

        if let Some(_) = self.condens_lvl {
            // we check the condensation level as it might be a level of free convection
            for (i, point) in parcel_log.iter().skip(current_point).enumerate() {
                // first time this is true is LFC
                if point.vrt_temp > env_vrt_tmp[i] {
                    self.lfc = Some(point.position.z);
                    current_point = i;
                    break;
                }
            }
        }

        if let Some(_) = self.lfc {
            // start checking from level after LFC for rare case when virtual temperatures are equal
            for (i, point) in parcel_log.iter().skip(current_point + 1).enumerate() {
                // first time this is true is LFC
                if point.vrt_temp <= env_vrt_tmp[i] {
                    self.el = Some(point.position.z);
                    break;
                }
            }
        }

        Ok(())
    }

    fn update_thermodynamic_vars(
        &mut self,
        parcel_log: &Vec<ParcelState>,
        env_vrt_tmp: &Vec<Float>,
    ) {
        let mut lfc_id = 0;

        // compute CIN if LFC is present
        let mut cin: Float = 0.0;
        if let Some(_) = self.lfc {
            //we start from the 2nd point of parcel log to not go out of bounds
            for (i, point) in parcel_log.iter().skip(1).enumerate() {
                let y_1 = (point.vrt_temp - env_vrt_tmp[i]) / env_vrt_tmp[i];
                let y_0 = (parcel_log[i - 1].vrt_temp - env_vrt_tmp[i - 1]) / env_vrt_tmp[i - 1];

                let delta_z = point.position.z - parcel_log[i - 1].position.z;

                cin += ((y_0 + y_1) / 2.0) * delta_z;

                if point.position.z == self.lfc.unwrap() {
                    lfc_id = i;
                    break;
                }
            }
        }

        self.cin = Some(-G * cin);

        // compute CAPE if LFC and EL is present
        let mut cape: Float = 0.0;
        if let Some(_) = self.lfc {
            if let Some(_) = self.el {
                // we start integration from LFC
                for (i, point) in parcel_log.iter().skip(lfc_id + 1).enumerate() {
                    let y_1 = (point.vrt_temp - env_vrt_tmp[i]) / env_vrt_tmp[i];
                    let y_0 =
                        (parcel_log[i - 1].vrt_temp - env_vrt_tmp[i - 1]) / env_vrt_tmp[i - 1];

                    let delta_z = point.position.z - parcel_log[i - 1].position.z;

                    cape += ((y_0 + y_1) / 2.0) * delta_z;

                    if point.position.z == self.el.unwrap() {
                        break;
                    }
                }
            }
        }

        self.cape = Some(G * cape);
    }
}

fn get_env_vtemp(
    parcel_log: &Vec<ParcelState>,
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

    let env_temp = env_vtemp?;

    Ok(env_temp)
}
