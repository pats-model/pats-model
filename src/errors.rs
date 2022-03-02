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

//! Module with error definitions for all
//! struct and function in the model.

use crate::Float;
use thiserror::Error;

/// General errors gathering all errors that can be
/// returned by the model.
#[derive(Error, Debug)]
pub enum ModelError {
    #[error("Error while reading config.yaml: {0}")]
    Config(#[from] ConfigError),

    #[error("Error while creating ThreadPool: {0}")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError),

    #[error("Error occured in Environment struct: {0}")]
    Environment(#[from] EnvironmentError),

    #[error("Error while doing thermodynamic computation, check your input data: {0}")]
    UnreasonableVariable(#[from] floccus::errors::InputError),

    #[error("Error while handling the file: {0}")]
    FileHandling(#[from] std::io::Error),

    #[error("Error with output directory: {0}")]
    FaultyOutput(&'static str),
}

/// Errors related to reading and handling the model configuration.
#[derive(Error, Debug)]
pub enum ConfigError {
    #[error("Cannot open config.yaml: {0}")]
    CantOpenFile(#[from] std::io::Error),

    #[error("Cannot deserialize config.yaml: {0}")]
    CantDeserialize(#[from] serde_yaml::Error),

    #[error("Configuration component is out of bounds: {0}")]
    OutOfBounds(&'static str),

    #[error("Error while reading GRIB input: {0}")]
    CannotReadInput(#[from] InputError),
}

/// Errors related to reading and handling
/// boundary conditions (environment data).
#[derive(Error, Debug)]
pub enum EnvironmentError {
    #[error(
        "Error while handling projection ({0}), please check your domain or report it on Github"
    )]
    ProjectionError(#[from] ProjectionError),

    #[error("Error of GRIB input handling: {0}")]
    GRIBInput(#[from] InputError),

    #[error("Could not find the value using bisection: {0}")]
    SearchUnable(#[from] SearchError),
}

/// Errors related to reading input GRIB files.
#[derive(Error, Debug)]
pub enum InputError {
    #[error("Error while reading the GRIB file: {0}")]
    CannotReadGrib(#[from] eccodes::errors::CodesError),

    #[error("Error while parsing string to datetime: {0}")]
    CannotParseDatetime(#[from] chrono::format::ParseError),

    #[error("The type of key {0} is incorrect")]
    IncorrectKeyType(&'static str),

    #[error("Provided input data is not sufficient to run the model, please check the documentation: {0}")]
    DataNotSufficient(&'static str),

    #[error("Values shape mismatch in GRIB, please check your input data: {0}")]
    IncorrectShape(#[from] ndarray::ShapeError),
}

/// Errors related to searching datasets with bisection.
#[derive(Error, Debug)]
pub enum SearchError {
    #[error("Provided array is empty")]
    EmptyArray,

    #[error("Provided target is out of array bounds")]
    OutOfBounds,
}

/// Errors related to parcel handling.
#[derive(Error, Debug)]
pub enum ParcelError {
    #[error("Error while doing thermodynamic computation, check your input data: {0}")]
    UnreasonableVariable(#[from] floccus::errors::InputError),

    #[error("Error while accessing environmental variable: {0}")]
    EnvironmentAccess(#[from] EnvironmentError),

    #[error("Error while handling the file: {0}")]
    FileHandling(#[from] std::io::Error),

    #[error("Error while handling the csv file: {0}")]
    CSVHandling(#[from] csv::Error),

    #[error("Parcel released from N{0:.3} E{1:.3} has stopped its ascent with error: {2} Check your configuration.")]
    AscentStopped(Float, Float, ParcelSimulationError),
}

/// Errors related to parcel simulation.
#[derive(Error, Debug)]
pub enum ParcelSimulationError {
    #[error("Error while doing thermodynamic computation, check your input data: {0}")]
    UnreasonableVariable(#[from] floccus::errors::InputError),

    #[error("Error while accessing environmental variable: {0}")]
    EnvironmentAccess(#[from] EnvironmentError),
}

/// Errors realted to geographic projection.
#[derive(Error, Debug)]
pub enum ProjectionError {
    #[error("Incorrect projection parameters: {0}")]
    IncorrectParams(&'static str),
}
