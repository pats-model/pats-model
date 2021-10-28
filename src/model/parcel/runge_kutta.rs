use super::{ParcelState, Vec3};
use crate::errors::ParcelSimulationError;
use crate::model::environment::EnvFields::{Pressure, UWind, VWind, VirtualTemperature};
use crate::{model::environment::Environment, Float};
use chrono::Duration;
use floccus::constants::{C_P, C_PV, C_V, C_VV, EPSILON, G, L_V, R_D};
use floccus::{mixing_ratio, vapour_pressure, virtual_temperature};
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
        let mut parcel_log = vec![];
        parcel_log.push(initial_state);

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
            let ref_parcel = self.parcel_log.last().unwrap().clone();

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

        if initial_state.velocity.z <= 0.0 || initial_state.mxng_rto < 0.000001 {
            return Ok(());
        }

        debug!("Starting pseudoadiabatic ascent");
        debug!("Init state: {:?}", initial_state);

        let pseudoadiabatic_scheme = PseudoAdiabaticScheme::new(initial_state, self.env);

        loop {
            let ref_parcel = self.parcel_log.last().unwrap().clone();

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

            if result_parcel.velocity.z <= 0.0 || result_parcel.mxng_rto < 0.000001 {
                break;
            }

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

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Clone, Debug)]
struct AdiabaticScheme<'a> {
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
        let mut updated_state = ref_state.clone();

        updated_state.pres = self.env.get_field_value(
            ref_state.position.x,
            ref_state.position.y,
            ref_state.position.z,
            Pressure,
        )?;

        updated_state.temp =
            (self.lambda / updated_state.pres.powf(1.0 - self.gamma)).powf(1.0 / self.gamma);

        let satr_vap_pres;
        if updated_state.temp < 253.0 {
            // if the temperature is very low use dedicated formula
            satr_vap_pres = vapour_pressure::buck2(updated_state.temp, updated_state.pres)?;
        } else {
            // else use usual buck formula over water
            satr_vap_pres = vapour_pressure::buck1(updated_state.temp, updated_state.pres)?;
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
struct PseudoAdiabaticScheme<'a> {
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
    pub fn state_at_position(
        &self,
        ref_state: &ParcelState,
    ) -> Result<ParcelState, ParcelSimulationError> {
        let mut updated_state = ref_state.clone();

        updated_state.pres = self.env.get_field_value(
            ref_state.position.x,
            ref_state.position.y,
            ref_state.position.z,
            Pressure,
        )?;

        updated_state.temp = self.iterate_to_temperature(updated_state.pres);

        let satr_vap_pres;
        if updated_state.temp < 253.0 {
            // if the temperature is very low use dedicated formula
            satr_vap_pres = vapour_pressure::buck2(updated_state.temp, updated_state.pres)?;
        } else {
            // else use usual buck formula over water
            satr_vap_pres = vapour_pressure::buck1(updated_state.temp, updated_state.pres)?;
        }

        updated_state.satr_mxng_rto = mixing_ratio::general1(updated_state.pres, satr_vap_pres)?;

        // if saturation mixing ratio dropped we bring the parcel back to
        // 100% saturation
        if updated_state.satr_mxng_rto > updated_state.mxng_rto {
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
