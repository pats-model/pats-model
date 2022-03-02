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

use super::ParcelState;
use crate::{
    errors::{EnvironmentError, ParcelError},
    model::{
        environment::{
            EnvFields::{Temperature, VirtualTemperature},
            Environment,
        },
        vec3::Vec3,
    },
    Float,
};
use chrono::NaiveDateTime;
use std::{path::Path, sync::Arc};

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
struct AnnotatedParcelState {
    datetime: NaiveDateTime,
    lon: Float,
    lat: Float,
    height: Float,
    velocity: Vec3,
    pres: Float,
    temp: Float,
    mxng_rto: Float,
    satr_mxng_rto: Float,
    vrt_temp: Float,
    env_temp: Float,
    env_vrt_temp: Float,
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub(super) fn save_parcel_log(
    parcel_log: &[ParcelState],
    environment: &Arc<Environment>,
) -> Result<(), ParcelError> {
    let parcel_id = construct_parcel_id(parcel_log.first().unwrap(), environment);

    let parcel_log = annotate_parcel_log(parcel_log, environment)?;

    let out_path = format!("./output/{}.csv", parcel_id);
    let out_path = Path::new(&out_path);

    let mut out_file = csv::Writer::from_path(out_path)?;

    out_file.write_record(&[
        "dateTime",
        "longitude",
        "latitude",
        "height",
        "velocityX",
        "velocityY",
        "velocityZ",
        "pressure",
        "temperature",
        "mixingRatio",
        "saturationMixingRatio",
        "virtualTemperature",
        "envTemperature",
        "envVirtualTemperature",
    ])?;

    for parcel in parcel_log {
        out_file.write_record(&[
            parcel.datetime.to_string(),
            parcel.lon.to_string(),
            parcel.lat.to_string(),
            parcel.height.to_string(),
            parcel.velocity.x.to_string(),
            parcel.velocity.y.to_string(),
            parcel.velocity.z.to_string(),
            parcel.pres.to_string(),
            parcel.temp.to_string(),
            parcel.mxng_rto.to_string(),
            parcel.satr_mxng_rto.to_string(),
            parcel.vrt_temp.to_string(),
            parcel.env_temp.to_string(),
            parcel.env_vrt_temp.to_string(),
        ])?;
    }

    out_file.flush()?;

    Ok(())
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn annotate_parcel_log(
    parcel_log: &[ParcelState],
    environment: &Arc<Environment>,
) -> Result<Vec<AnnotatedParcelState>, EnvironmentError> {
    let mut result_log = Vec::<AnnotatedParcelState>::with_capacity(parcel_log.len());

    for parcel in parcel_log {
        let (lon, lat) = environment
            .projection
            .inverse_project(parcel.position.x, parcel.position.y);

        let env_temp = environment.get_field_value(
            parcel.position.x,
            parcel.position.y,
            parcel.position.z,
            Temperature,
        )?;

        let env_vrt_temp = environment.get_field_value(
            parcel.position.x,
            parcel.position.y,
            parcel.position.z,
            VirtualTemperature,
        )?;

        result_log.push(AnnotatedParcelState {
            datetime: parcel.datetime,
            lon,
            lat,
            height: parcel.position.z,
            velocity: parcel.velocity,
            pres: parcel.pres,
            temp: parcel.temp,
            mxng_rto: parcel.mxng_rto,
            satr_mxng_rto: parcel.satr_mxng_rto,
            vrt_temp: parcel.vrt_temp,
            env_temp,
            env_vrt_temp,
        });
    }

    Ok(result_log)
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn construct_parcel_id(initial_state: &ParcelState, environment: &Arc<Environment>) -> String {
    let time_stamp = initial_state.datetime.format("%Y-%m-%dT%H%M%S").to_string();
    let (lon, lat) = environment
        .projection
        .inverse_project(initial_state.position.x, initial_state.position.y);

    let position_stamp = format!("N{:.4}_E{:.4}", lon, lat);

    format!("parcel_{}_{}", position_stamp, time_stamp)
}
