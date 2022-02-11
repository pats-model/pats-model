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
use crate::model::environment::Environment;
use std::{io::Error, path::Path, sync::Arc};

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub(super) fn save_parcel_log(
    parcel_log: &[ParcelState],
    environment: &Arc<Environment>,
) -> Result<(), Error> {
    let parcel_id = construct_parcel_id(parcel_log.first().unwrap());

    let out_path = format!("./output/{}.csv", parcel_id);
    let out_path = Path::new(&out_path);

    let mut out_file = csv::Writer::from_path(out_path)?;

    out_file.write_record(&[
        "dateTime",
        "longitude",
        "latitude",
        "positionZ",
        "velocityX",
        "velocityY",
        "velocityZ",
        "pressure",
        "temperature",
        "mixingRatio",
        "saturationMixingRatio",
        "virtualTemperature",
    ])?;

    for parcel in parcel_log {
        let (lon, lat) = environment
            .projection
            .inverse_project(parcel.position.x, parcel.position.y);

        out_file.write_record(&[
            parcel.datetime.to_string(),
            lon.to_string(),
            lat.to_string(),
            parcel.position.z.to_string(),
            parcel.velocity.x.to_string(),
            parcel.velocity.y.to_string(),
            parcel.velocity.z.to_string(),
            parcel.pres.to_string(),
            parcel.temp.to_string(),
            parcel.mxng_rto.to_string(),
            parcel.satr_mxng_rto.to_string(),
            parcel.vrt_temp.to_string(),
        ])?;
    }

    out_file.flush()?;

    Ok(())
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn construct_parcel_id(initial_state: &ParcelState) -> String {
    let time_stamp = initial_state.datetime.format("%Y-%m-%dT%H%M%S").to_string();
    let position_stamp = format!(
        "x-{}_y-{}",
        initial_state.position.x as i64, initial_state.position.y as i64
    );

    format!("parcel_{}_{}", position_stamp, time_stamp)
}
