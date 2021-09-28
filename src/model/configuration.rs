use crate::errors::ConfigError;
use serde::Deserialize;
use std::{
    fs,
    path::{Path, PathBuf},
};

use super::Float;

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Deserialize)]
pub struct Domain {
    pub lat_0: Float,
    pub lon_0: Float,
    pub spacing: Float,
    pub ni: u16,
    pub nj: u16,
}

impl Domain {
    pub fn check_bounds(domain: &Domain) -> Result<(), ConfigError> {
        if !(-90.0..90.0).contains(&domain.lat_0) {
            return Err(ConfigError::OutOfBounds(
                "Reference latitude is too low or too high",
            ));
        }

        if !(-180.0..180.0).contains(&domain.lon_0) {
            return Err(ConfigError::OutOfBounds(
                "Reference longitude is too low or too high",
            ));
        }

        Ok(())
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Deserialize)]
pub struct Config {
    pub domain: Domain,
    pub level_type: String,
    pub timestep: Float,
    pub data_files: Vec<PathBuf>,
    pub threads: u16,
    pub memory: u32,
}

impl Config {
    pub fn new_from_file(file_path: &Path) -> Result<Config, ConfigError> {
        let data = fs::read(file_path)?;
        let config: Config = serde_yaml::from_slice(data.as_slice())?;

        Domain::check_bounds(&config.domain)?;

        Ok(config)
    }
}
