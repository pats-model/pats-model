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

//! (TODO: What it is)
//! 
//! (Why it is neccessary)

use std::sync::Arc;

use crate::{errors::ParcelError, Float};

use super::{configuration::Config, environment::Environment};

/// (TODO: What it is)
/// 
/// (Why it is neccessary)
pub fn deploy(
    start_coords: (Float, Float),
    config: Arc<Config>,
    environment: Arc<Environment>,
) -> Result<(), ParcelError> {
    Ok(())
}
