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
use crate::{errors::ParcelError, model::environment::Environment};
use super::ParcelState;

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct ConvectiveParams {}

pub(super) fn compute_conv_params(
    parcel_log: &Vec<ParcelState>,
    environment: &Arc<Environment>,
) -> Result<ConvectiveParams, ParcelError> {
    Ok(ConvectiveParams {})
}
