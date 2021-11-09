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

//! Module containing constants used by the model.

use crate::Float;
use std::f64::consts::PI;

///WGS84 ellipsoid semi-major axis
pub const WGS84_A: Float = 6_378_137.0;

///WGS84 ellipsoid semi-minor axis
#[allow(clippy::excessive_precision)]
pub const WGS84_B: Float = 6_356_752.314_245;

///WGS84 ellipsoid eccentricity
#[allow(clippy::excessive_precision)]
pub const WGS84_E: Float =
    0.081_819_190_842_965_558_441_157_725_155_790_103_599_429_130_554_199_218_75;

///WGS84 ellipsoid Ramanujan's $h$ parameter
pub const WGS84_H: Float =
    ((WGS84_A - WGS84_B) * (WGS84_A - WGS84_B)) / ((WGS84_A + WGS84_B) * (WGS84_A - WGS84_B));

///WGS84 ellipsoid circumference along meridian
///
///Computed with first 6 terms of infinite series:
///`C = \pi(a+b)\sum_{n=0}^{+\infty}\binom{0.5}{n}h^n`
pub const NS_C_EARTH: Float = PI
    * (WGS84_A + WGS84_B)
    * (1.0
        + (1.0 / 4.0) * (WGS84_H)
        + (1.0 / 64.0) * (WGS84_H * WGS84_H)
        + (1.0 / 256.0) * (WGS84_H * WGS84_H * WGS84_H)
        + (25.0 / 16384.0) * (WGS84_H * WGS84_H * WGS84_H * WGS84_H)
        + (49.0 / 65536.0) * (WGS84_H * WGS84_H * WGS84_H * WGS84_H * WGS84_H));

///WGS84 ellipsoid circumference along equator
pub const WE_C_EARTH: Float = 2.0 * PI * WGS84_A;
