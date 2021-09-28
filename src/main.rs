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

mod errors;
mod model;

use env_logger::Env;
use log::{error, info};

fn main() {
    let logger_env = Env::new().filter_or("PATS_LOG_LEVEL", "info");
    env_logger::Builder::from_env(logger_env)
        .format_timestamp_millis()
        .init();

    match model::main() {
        Ok(out_name) => info!(
            "Model finished successfully. Check the output file {}",
            out_name
        ),
        Err(err) => error!("Model failed with error: {}", err),
    }
}
