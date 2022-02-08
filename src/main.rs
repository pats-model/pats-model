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

//! Technical documentation of Parcel Ascent Tracing System (PATS) -
//! the numerical model for convective parcel ascent simulation in three-dimensions.
//! 
//! This documentation provides a description of functions and structures
//! used in the model. Its main purpose is to make it easier to maintain
//! and contribute to the project codebase. However, it can be also useful
//! for users who want to understand the model in more detail.

mod constants;
mod errors;
mod model;

use cap::Cap;
use env_logger::Env;
use log::{error, info};
use std::alloc;


type Float = f64;

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
    #[cfg(not(feature="debug"))]
    let logger_env = Env::new().filter_or("PATS_LOG_LEVEL", "info");

    #[cfg(feature="debug")]
    let logger_env = Env::new().filter_or("PATS_LOG_LEVEL", "debug");
    
    env_logger::Builder::from_env(logger_env)
        .format_timestamp_millis()
        .init();

    match model::main() {
        Ok(_) => info!("Model execution finished. Check the output directory and log."),
        Err(err) => error!("Model execution failed with error: {}", err),
    }
}
