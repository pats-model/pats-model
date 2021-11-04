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

//! Module responsible for parsing and checking the configuration file.
//!
//! To provide meaningful error messages. The configuration file uses
//! [YAML](https://en.wikipedia.org/wiki/YAML) and `serde` to enforce
//! strong typing and automatic type checking.
//!
//! The structures and their fields in this module directly correspond to
//! the fields inside `config.yaml` so you can check this documentation
//! for more details how to set the config file.

use crate::errors::ConfigError;
use chrono::NaiveDateTime;
use serde::Deserialize;
use std::{
    fs,
    path::{Path, PathBuf},
};

use crate::Float;

/// Fields with model domain information.
///
/// Model domain is defined as the area from which parcels
/// start their plus margins for parcels released near the domain edge.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Deserialize)]
pub struct Domain {
    /// Longitude (in degrees) of south-west domain corner.
    ///
    /// Must meet the condition: `-180 < ref_lon < 180`
    pub ref_lon: Float,

    /// Latitude (in degrees) of south-west domain corner.
    ///
    /// Must meet the condition: `-90 < ref_lon < 90`
    pub ref_lat: Float,

    /// Domain spacing in meters. Represents the distance between parcels
    /// in x and y directions.
    /// 
    /// Cannot be smaller than `1`.
    pub spacing: Float,

    /// Domain shape (in model gridpoints/parcels). Represents
    /// how much parcels will be released along each axis.
    /// 
    /// Total number of released parcels cannot be smaller than `1`.
    pub shape: (u16, u16),

    /// _(Optional)_ Domain margins (in degrees) for lon and lat
    /// axis respectively. Parcels will not be released in the margins
    /// area, but the input data will be read there so that parcels can use it.
    /// 
    /// Defaults to `1.0`. Cannot be less than `0.1`.
    #[serde(default = "Domain::default_margins")]
    pub margins: (Float, Float),
}

impl Domain {
    /// Checks if domain specification follows conventions
    /// and limits.
    pub fn check_bounds(&self) -> Result<(), ConfigError> {
        if !(-90.0..90.0).contains(&self.ref_lat) {
            return Err(ConfigError::OutOfBounds(
                "Reference latitude is too low or too high",
            ));
        }

        if !(-180.0..180.0).contains(&self.ref_lon) {
            return Err(ConfigError::OutOfBounds(
                "Reference longitude is too low or too high",
            ));
        }

        if (u64::from(self.shape.0) *  u64::from(self.shape.1)) < 1 {
            return Err(ConfigError::OutOfBounds(
                "Total number of gridpoints cannot be less than 1",
            ));
        }

        if self.spacing < 1.0 {
            return Err(ConfigError::OutOfBounds(
                "Grid spacing cannot be smaller than 1 m",
            ));
        }

        if self.margins.0 < 0.1 || self.margins.1 < 0.1 {
            return Err(ConfigError::OutOfBounds(
                "Margins cannot be smaller than 0.1 degree",
            ));
        }

        Ok(())
    }

    fn default_margins() -> (Float, Float) {
        (1.0, 1.0)
    }
}

/// Fields with information about time used by model.
#[derive(Clone, PartialEq, PartialOrd, Debug, Deserialize)]
pub struct DateTime {
    /// Timestep (in seconds) used by the model.
    ///
    /// This setting will directly affect the execution time,
    /// memory usage, output size of the model, as well as its
    /// accuracy and stability.
    pub timestep: Float,

    /// Start datetime for the model. Currently, it is used
    /// only as a reference to provide more helpful output
    /// and does not affect background conditions.
    pub start: NaiveDateTime,
}

/// Fields with information about model input data
/// for providing boundary conditions.
#[derive(Clone, PartialEq, PartialOrd, Debug, Deserialize)]
pub struct Input {
    /// Level type of GRIB messages in input files from
    /// which 3D boundary conditions data should be read.
    ///
    /// Currently the model will only work when this field
    /// is set to "isobaricInhPa".
    pub level_type: String,

    /// List of input GRIB files to read boundary coonditions.
    ///
    /// Currently those files must meet following criteria:
    ///
    /// - Data inside files must cover at least whole with margins.
    /// - Required variables for surface levels are: temperature, dewpoint,
    /// u and v wind components, pressure and geopotential.
    /// - Required variables for pressure levels are: temperature, geopotential,
    /// specific humidity and u and v wind components.
    /// - For each variable all levels must be unique.
    /// - Files must contain data only for one datetime.
    /// - None of the files can be empty.
    /// - Ideally, there should be only data actually used by model in files.
    pub data_files: Vec<PathBuf>,
}

/// _(Optional)_ Fields with information about
/// resources available for model.
#[derive(Clone, PartialEq, PartialOrd, Debug, Deserialize)]
pub struct Resources {
    /// _(Optional)_ Thread count used by the model.
    /// The thread pool initiated by this model will use
    /// up to this number of workers.
    ///
    /// Cannot be less than `1`. Defaults to `1`.
    #[serde(default = "Resources::default_threads")]
    pub threads: u16,

    /// _(Optional)_ Heap memory limit for the model in MB.
    /// Useful for enabling meaningful Out-of-memory error messages.
    ///
    /// Cannot be less than `128`. Defaults to whole addressable-space
    /// (`2^32` or `2^64` bytes).
    ///
    /// Currently, Rust memory allocator aborts the program
    /// when OOM occurs. By default, Rust memory allocator
    /// doesn't know (read: doesn't care) about memory available
    /// in your system, and defaults its memory limit to a whole
    /// addresable space which most likely (read: for sure) exceeds
    /// your system memory.
    ///
    /// This model allocate substantial amount of memory in small chunks,
    /// which makes OOM unlikely to occur. Instead, your system will slow down
    /// and eventually kill the process without any additional information.
    /// However, when the allocator has a capped memory amount available
    /// it will abort the process with (somehow useful) OOM error message.
    ///
    /// If the model gets killed for no apparent reason try setting the
    /// memory limit lower than your avilable system memory and check if
    /// OOM error occurs. Be generous when setting the limit but leave some
    /// space for other processes.
    #[serde(default = "Resources::default_memory")]
    pub memory: usize,
}

impl Resources {
    fn default_threads() -> u16 {
        1
    }

    fn default_memory() -> usize {
        usize::MAX / (1024 * 1024)
    }

    /// Checks if thread count and memory limit are
    /// above limits.
    pub fn check_bounds(&self) -> Result<(), ConfigError> {
        if self.threads < 1 {
            return Err(ConfigError::OutOfBounds(
                "Available threads cannot be less than 1",
            ));
        }

        if self.memory < 128 {
            return Err(ConfigError::OutOfBounds(
                "Available memory cannot be less than 128 MB",
            ));
        }

        Ok(())
    }
}

impl Default for Resources {
    fn default() -> Self {
        Resources {
            threads: Resources::default_threads(),
            memory: Resources::default_memory(),
        }
    }
}

/// Main config structure representing the fields in
/// configuration file.
#[derive(Clone, PartialEq, PartialOrd, Debug, Deserialize)]
pub struct Config {
    pub domain: Domain,

    pub datetime: DateTime,

    pub input: Input,

    #[serde(default)]
    pub resources: Resources,
}

impl Config {
    /// Config structure constructor, responsible for
    /// deserializing configuration and checking it.
    pub fn new_from_file(file_path: &Path) -> Result<Config, ConfigError> {
        let data = fs::read(file_path)?;
        let config: Config = serde_yaml::from_slice(data.as_slice())?;

        config.domain.check_bounds()?;
        config.resources.check_bounds()?;

        Ok(config)
    }
}
