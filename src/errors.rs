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

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ModelError {
    #[error("Error while reading config.yaml: {0}")]
    Config(#[from] ConfigError),

    #[error("Error while creating ThreadPool: {0}")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError)
}   

#[derive(Error, Debug)]
pub enum ConfigError {
    #[error("Cannot open config.yaml: {0}")]
    CantOpenFile(#[from] std::io::Error),

    #[error("Cannot deserialize config.yaml: {0}")]
    CantDeserialize(#[from] serde_yaml::Error),

    #[error("Configuration component is out of bounds {0}")]
    OutOfBounds(&'static str),
}