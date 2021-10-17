use super::{ParcelState, Vec3};
use crate::model::environment::EnvFields::{Pressure, VirtualTemperature};
use crate::{errors::ParcelError, model::environment::Environment, Float};
use floccus::constants::{C_P, C_PV, C_V, C_VV, G};
use floccus::{mixing_ratio, virtual_temperature};
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
    pub fn run_simulation(&mut self) -> Result<(), ParcelError> {
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
    fn ascent_adiabatically(&mut self) -> Result<(), ParcelError> {
        let initial_state = self.parcel_log.last().unwrap();
        let adiabatic_scheme = AdiabaticScheme::new(initial_state, self.env);

        loop {
            let ref_parcel = self.parcel_log.last().unwrap().clone();

            if ref_parcel.velocity.z <= 0.0 && ref_parcel.mxng_rto <= ref_parcel.satr_mxng_rto {
                break;
            }

            // holographic parcel is a virtual parcel that is moved
            // around for RK4 computations but doesn't change its
            // thermodynamic properties in reference to the prestep state
            let holo_parcel = ref_parcel;
            let c_0 = holo_parcel.velocity;
            let k_0 =
                self.calculate_bouyancy_force(&adiabatic_scheme.state_at_position(&holo_parcel)?)?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += 0.5 * self.timestep * c_0;
            let c_1 = holo_parcel.velocity + 0.5 * self.timestep * k_0;
            let k_1 =
                self.calculate_bouyancy_force(&adiabatic_scheme.state_at_position(&holo_parcel)?)?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += 0.5 * self.timestep * c_1;
            let c_2 = holo_parcel.velocity + 0.5 * self.timestep * k_1;
            let k_2 =
                self.calculate_bouyancy_force(&adiabatic_scheme.state_at_position(&holo_parcel)?)?;

            let mut holo_parcel = ref_parcel;
            holo_parcel.position += self.timestep * c_2;
            let c_3 = holo_parcel.velocity + self.timestep * k_2;
            let k_3 =
                self.calculate_bouyancy_force(&adiabatic_scheme.state_at_position(&holo_parcel)?)?;

            let delta_pos = (self.timestep / 6.0) * (c_0 + 2.0 * c_1 + 2.0 * c_2 + c_3);
            let delta_vel = (self.timestep / 6.0) * (k_0 + 2.0 * k_1 + 2.0 * k_2 + k_3);

            let mut result_parcel = ref_parcel;
            result_parcel.position += delta_pos;
            result_parcel.velocity += delta_vel;
            adiabatic_scheme.state_at_position(&result_parcel)?;

            self.parcel_log.push(result_parcel);
        }

        Ok(())
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn ascent_pseudoadiabatically(&mut self) -> Result<(), ParcelError> {
        Ok(())
    }

    /// (TODO: What it is)
    ///
    /// (Why it is neccessary)
    fn calculate_bouyancy_force(&self, parcel: &ParcelState) -> Result<Vec3, ParcelError> {
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
    pub fn state_at_position(&self, ref_state: &ParcelState) -> Result<ParcelState, ParcelError> {
        let mut updated_state = ref_state.clone();

        updated_state.pres = self.env.get_field_value(
            ref_state.position.x,
            ref_state.position.y,
            ref_state.position.z,
            Pressure,
        )?;

        updated_state.temp =
            (self.lambda / updated_state.pres.powf(1.0 - self.gamma)).powf(1.0 / self.gamma);

        updated_state.satr_mxng_rto =
            mixing_ratio::accuracy1(updated_state.temp, updated_state.pres)?;
        updated_state.vrt_temp =
            virtual_temperature::general1(updated_state.temp, updated_state.mxng_rto)?;

        Ok(updated_state)
    }
}
