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

//! Parcel Ascent Tracing System (PATS) is the numerical model
//! for convective parcel ascent simulation in three-dimension
//! developed as a Master's project at the School of Environmental
//! Sciences at the University of East Anglia.
//!
//! The model predicts how convective air parcel behaves in conditions
//! much more realistic then those assumed by classical parcel theory.
//!
//! **The model is currently under development.**

mod constants;
mod errors;
mod model;

use cap::Cap;
use env_logger::Env;
use log::{error, info};
use std::alloc;

#[cfg(feature = "double_precision")]
type Float = f64;

/// Type used by the model floating point variables.
///
/// ECMWF IFS developers [have shown](https://www.ecmwf.int/en/newsletter/148/meteorology/single-precision-ifs)
/// that single-precision (32 bit) floating point variables are sufficient in most cases.
/// However, the increased precision (at a cost of performance) might me necessary sometimes.
/// Therefore all floating point variables in the model are of this type.
/// It defaults to [`f32`] and changes to [`f64`] when `double_precision` feature is used.
#[cfg(not(feature = "double_precision"))]
type Float = f32;

/// Global allocator used by the model.
///
/// Use of static global allocator allows for capping the memory to the limit set by user
/// in configuration file and in effect provide better [OOM error](https://en.wikipedia.org/wiki/Out_of_memory) handling.
#[global_allocator]
static ALLOCATOR: Cap<alloc::System> = Cap::new(alloc::System, usize::MAX);

/// The main program function.
/// Prepares the runtime environment and calls the [`model::main`].
///
/// To provide meaningful and high-quality error messages the `env_logger`
/// needs to be initiated before any log messages are possible to occur.
/// Furthermore, errors can occur also during model shutdown and they also
/// can be handled.
fn main() {
    let logger_env = Env::new().filter_or("PATS_LOG_LEVEL", "info");
    env_logger::Builder::from_env(logger_env)
        .format_timestamp_millis()
        .init();

    match model::main() {
        Ok(_) => info!("Model execution finished. Check the output directory and log."),
        Err(err) => error!("Model execution failed with error: {}", err),
    }
}
