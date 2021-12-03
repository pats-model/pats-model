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

//! Module containing the actual model code.
//! Whole documentation of how the model works is provided here.

mod configuration;
mod environment;
mod parcel;
mod vec3;

#[cfg(test)]
mod super_tests;

use crate::model::parcel::conv_params::ConvectiveParams;
use crate::{
    errors::ModelError,
    model::{configuration::Config, environment::Environment},
    Float, ALLOCATOR,
};
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, error, info};
use ndarray::Array1;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::{
    fs,
    io::Error,
    path::Path,
    sync::{mpsc, Arc},
};

/// Convenience type to store lon-lat coordinates.
type LonLat<T> = (T, T);

/// Main model function, responsible for all simulation steps.
///
/// It reads the provided configuration and input data
/// and then deploys parcels within the domain onto the threadpool
/// and checks for errors.
pub fn main() -> Result<(), ModelError> {
    info!("Preparing the model core");

    // prepare all prerequisites for running the model
    prepare_output_dir()?;

    let model_core = Core::new()?;

    let parcels = prepare_parcels_list(&model_core);
    let parcels_count = parcels.len();

    let mut parcels_params: Vec<ConvectiveParams> = Vec::with_capacity(parcels_count);

    let config = Arc::new(model_core.config);
    let environment = Arc::new(model_core.environ);

    info!("Deploying parcels");

    // set progress bar for simulated parcels
    let parcels_bar = ProgressBar::new(parcels_count as u64);
    parcels_bar.set_style(
        ProgressStyle::default_bar()
            .template("{prefix} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} {msg}")
            .progress_chars("#>-"),
    );
    parcels_bar.set_prefix("Simulated parcels");

    // deploy parcels on to the threadpool
    let (tx, rx) = mpsc::channel();

    for parcel_coords in parcels {
        let tx = tx.clone();
        let config = Arc::clone(&config);
        let environment = Arc::clone(&environment);

        model_core.threadpool.spawn(move || {
            tx.send(parcel::deploy(parcel_coords, &config, &environment))
                .unwrap();
        });
    }

    // receive parcels status and computed convective parameters
    for _ in 0..parcels_count {
        let parcel_result = rx.recv().expect("Receiving parcel result failed");

        match parcel_result {
            Ok(params) => {
                parcels_params.push(params);
            }
            Err(err) => {
                error!("Parcel simulation handling failed due to an error, check the details and rerun the model: {}", err);
                // this is neccessary to make sure that all error messages
                // are fully written to stdout before the progress bar updates
                println!("");
            }
        }
        parcels_bar.inc(1);
    }

    parcels_bar.finish_with_message("All parcels finished");
    info!("Writing output");

    //write convective parameters to file
    save_conv_params(parcels_params)?;

    Ok(())
}

/// Structure containing model parameters.
///
/// To run the simulation model needs to load and compute some initial
/// datam which is then stored in this structure.
#[derive(Debug)]
pub struct Core {
    pub config: Config,
    pub threadpool: ThreadPool,
    pub environ: Environment,
}

impl Core {
    /// Model [`Core`] constructor.
    ///
    /// Before the simulation can start (and to run it safely),
    /// configuration and input data provided by the user must be
    /// loaded and checked.
    pub fn new() -> Result<Self, ModelError> {
        debug!("Reading configuration from config.yaml");
        let mut config = Config::new_from_file(Path::new("config.yaml"))?;

        debug!("Setting memory limit");
        ALLOCATOR
            .set_limit(config.resources.memory * 1024 * 1024)
            .unwrap();

        debug!("Setting up ThreadPool");
        let threadpool = ThreadPoolBuilder::new()
            .num_threads(config.resources.threads as usize)
            .stack_size(2 * 1024 * 1024)
            .build()?;

        debug!("Reading environmental boundary conditions from GRIB");
        let environ = Environment::new(&mut config)?;

        Ok(Core {
            config,
            threadpool,
            environ,
        })
    }
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn prepare_output_dir() -> Result<(), ModelError> {
    debug!("Checking and setting output directory");

    let out_path = Path::new("./output/");

    if out_path.is_dir() {
        if out_path.read_dir()?.next().is_none() {
            debug!("Output directory exists but is empty so continuing");
        } else {
            return Err(ModelError::FaultyOutput(
                "Output directory exists and is not empty",
            ));
        }
    } else {
        debug!("Output directory does not exist so creating a new one");
        fs::create_dir(out_path)?;
    }

    Ok(())
}

/// Function calculating initial parcels positions from configuration
/// and gathering it into a list.
///
/// In configuration only south-west corner of the domain is provided.
/// Thus it is neccessary to compute the starting position of each parcel.
fn prepare_parcels_list(model_core: &Core) -> Vec<(Float, Float)> {
    let domain_anchor = model_core.environ.projection.project(
        model_core.config.domain.ref_lon,
        model_core.config.domain.ref_lat,
    );

    let x_coords = Array1::linspace(
        domain_anchor.0,
        domain_anchor.0
            + (Float::from(model_core.config.domain.shape.0 - 1)
                * model_core.config.domain.spacing),
        model_core.config.domain.shape.0 as usize,
    )
    .to_vec();

    let y_coords = Array1::linspace(
        domain_anchor.1,
        domain_anchor.1
            + (Float::from(model_core.config.domain.shape.1 - 1)
                * model_core.config.domain.spacing),
        model_core.config.domain.shape.1 as usize,
    )
    .to_vec();

    let mut xy_coords = vec![];

    for x in &x_coords {
        for y in &y_coords {
            xy_coords.push((*x, *y));
        }
    }

    xy_coords
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
fn save_conv_params(convective_params_list: Vec<ConvectiveParams>) -> Result<(), Error> {
    let out_path = Path::new("./output/model_convective_params.csv");

    let mut out_file = csv::Writer::from_path(out_path)?;

    out_file.write_record(&[
        "Lon",
        "Lat",
        "TopHeight",
        "xHorizontalDisplacement",
        "yHorizontalDisplacement",
        "maxVerticalVelocity",
        "condensationLevel",
        "levelOfFreeConvection",
        "equilibriumLevel",
        "CAPE",
        "CIN",
    ])?;

    for conv_params in convective_params_list {
        out_file.serialize(conv_params)?;
    }

    out_file.flush()?;

    Ok(())
}
