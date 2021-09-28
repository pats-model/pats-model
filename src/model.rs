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

mod configuration;

use crate::{errors::ModelError, model::configuration::Config};
use log::{debug, info};
use std::path::Path;

#[cfg(feature = "double_precision")]
type Float = f64;

#[cfg(not(feature = "double_precision"))]
type Float = f32;

pub fn main() -> Result<String, ModelError> {
    info!("Model computation started");

    debug!("Reading configuration from config.yaml");
    let global_config = Config::new_from_file(Path::new("config.yaml"))?;

    debug!("Setting up ThreadPool");
    rayon::ThreadPoolBuilder::new()
        .num_threads(global_config.threads as usize)
        .build_global()?;

    Ok(" ".to_string())
}
