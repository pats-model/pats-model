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
use crate::errors::ParcelSimulationError;
use crate::model::environment::EnvFields::Pressure;
use crate::{model::environment::Environment, Float};
use floccus::{
    constants::{C_P, C_PV, C_V, C_VV, EPSILON, L_V, R_D},
    mixing_ratio, vapour_pressure, virtual_temperature,
};
use std::sync::Arc;

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Clone, Debug)]
pub(super) struct AdiabaticScheme<'a> {
    lambda: Float,
    gamma: Float,
    env: &'a Arc<Environment>,
}

impl<'a> AdiabaticScheme<'a> {
    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    pub fn new(refrence: &ParcelState, environment: &'a Arc<Environment>) -> Self {
        let gamma = (C_P * ((1.0 + refrence.mxng_rto * (C_PV / C_P)) / (1.0 + refrence.mxng_rto)))
            / (C_V * ((1.0 + refrence.mxng_rto * (C_VV / C_V)) / (1.0 + refrence.mxng_rto)));

        let lambda = refrence.pres.powf(1.0 - gamma) * refrence.temp.powf(gamma);

        Self {
            lambda,
            gamma,
            env: environment,
        }
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    pub fn state_at_position(
        &self,
        ref_state: &ParcelState,
    ) -> Result<ParcelState, ParcelSimulationError> {
        let mut updated_state = *ref_state;

        updated_state.pres = self.env.get_field_value(
            ref_state.position.x,
            ref_state.position.y,
            ref_state.position.z,
            Pressure,
        )?;

        updated_state.temp =
            (self.lambda / updated_state.pres.powf(1.0 - self.gamma)).powf(1.0 / self.gamma);

        let satr_vap_pres;
        if updated_state.temp > 273.15 {
            // for most ranges use usual buck formula over water
            satr_vap_pres = vapour_pressure::buck1(updated_state.temp, updated_state.pres)?;
        } else if updated_state.temp > 193.0 {
            // if the temperature is very low use dedicated formula
            satr_vap_pres = vapour_pressure::buck2(updated_state.temp, updated_state.pres)?;
        } else {
            // as last resort if the temperature is very very low use more expensive dedicated formula
            satr_vap_pres = vapour_pressure::wexler2(updated_state.temp)?;
        }

        updated_state.satr_mxng_rto = mixing_ratio::general1(updated_state.pres, satr_vap_pres)?;
        updated_state.vrt_temp =
            virtual_temperature::general1(updated_state.temp, updated_state.mxng_rto)?;

        Ok(updated_state)
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Clone, Debug)]
pub(super) struct PseudoAdiabaticScheme<'a> {
    ref_temp: Float,
    ref_pres: Float,
    ref_mxng_rto: Float,
    ref_satr_mxng_rto: Float,
    env: &'a Arc<Environment>,
}

impl<'a> PseudoAdiabaticScheme<'a> {
    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    pub fn new(refrence: &ParcelState, environment: &'a Arc<Environment>) -> Self {
        PseudoAdiabaticScheme {
            ref_temp: refrence.temp,
            ref_pres: refrence.pres,
            env: environment,
            ref_mxng_rto: refrence.mxng_rto,
            ref_satr_mxng_rto: refrence.satr_mxng_rto,
        }
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    pub fn update_ref_state(&mut self, ref_state: &ParcelState) {
        self.ref_temp = ref_state.temp;
        self.ref_pres = ref_state.pres;
        self.ref_mxng_rto = ref_state.mxng_rto;
        self.ref_satr_mxng_rto = ref_state.satr_mxng_rto;
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    pub fn state_at_position(
        &self,
        ref_state: &ParcelState,
    ) -> Result<ParcelState, ParcelSimulationError> {
        let mut updated_state = *ref_state;

        updated_state.pres = self.env.get_field_value(
            ref_state.position.x,
            ref_state.position.y,
            ref_state.position.z,
            Pressure,
        )?;

        updated_state.temp = self.iterate_to_temperature(updated_state.pres);

        let satr_vap_pres;
        if updated_state.temp > 273.15 {
            // for most ranges use usual buck formula over water
            satr_vap_pres = vapour_pressure::buck1(updated_state.temp, updated_state.pres)?;
        } else if updated_state.temp > 193.0 {
            // if the temperature is very low use dedicated formula
            satr_vap_pres = vapour_pressure::buck2(updated_state.temp, updated_state.pres)?;
        } else {
            // as last resort if the temperature is very very low use more expensive dedicated formula
            satr_vap_pres = vapour_pressure::wexler2(updated_state.temp)?;
        }

        updated_state.satr_mxng_rto = mixing_ratio::general1(updated_state.pres, satr_vap_pres)?;

        // if saturation mixing ratio dropped we bring the parcel back to
        // 100% saturation
        if updated_state.satr_mxng_rto < updated_state.mxng_rto {
            updated_state.mxng_rto = updated_state.satr_mxng_rto;
        }

        updated_state.vrt_temp =
            virtual_temperature::general1(updated_state.temp, updated_state.mxng_rto)?;

        Ok(updated_state)
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn iterate_to_temperature(&self, target_pressure: Float) -> Float {
        let step_count = ((self.ref_pres - target_pressure).abs() / 1.0).ceil() as usize;
        let step = (target_pressure - self.ref_pres) / step_count as Float;

        let mut temp_n = self.ref_temp;
        let mut pres_n = self.ref_pres;

        // throughout the derivation we're keeping mixing ratios constant
        // as the derivative is a partial derivative of the pressure and temperature
        for _ in 0..step_count {
            let k_0 = pseudoadiabatic_derivative(
                temp_n,
                pres_n,
                self.ref_mxng_rto,
                self.ref_satr_mxng_rto,
            );
            let k_1 = pseudoadiabatic_derivative(
                temp_n + 0.5 * step * k_0,
                pres_n + 0.5 * step,
                self.ref_mxng_rto,
                self.ref_satr_mxng_rto,
            );
            let k_2 = pseudoadiabatic_derivative(
                temp_n + 0.5 * step * k_1,
                pres_n + 0.5 * step,
                self.ref_mxng_rto,
                self.ref_satr_mxng_rto,
            );
            let k_3 = pseudoadiabatic_derivative(
                temp_n + step * k_2,
                pres_n + step,
                self.ref_mxng_rto,
                self.ref_satr_mxng_rto,
            );

            pres_n += step;
            temp_n += (step / 6.0) * (k_0 + 2.0 * k_1 + 2.0 * k_2 + k_3);
        }

        temp_n
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn pseudoadiabatic_derivative(
    temp: Float,
    pres: Float,
    mxng_rto: Float,
    satr_mxng_rto: Float,
) -> Float {
    let b = (1.0 + (mxng_rto / EPSILON)) / (1.0 + (mxng_rto / (C_P / C_PV)));

    (b / pres)
        * ((R_D * temp + L_V * satr_mxng_rto)
            / (C_P + ((L_V * L_V * satr_mxng_rto * EPSILON * b) / (R_D * temp * temp))))
}
