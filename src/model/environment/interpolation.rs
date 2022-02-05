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

//! Module containing interpolation methods.

use crate::Float;
use nalgebra::{SMatrix, SVector};

pub type Vector16 = SVector<Float, 16>;
type Matrix16 = SMatrix<Float, 16, 16>;

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
pub struct Point2D {
    pub value: Float,

    pub x: Float,
    pub y: Float,
    pub dx: Float,
    pub dy: Float,
    pub dxdy: Float,
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
pub struct Point3D {
    pub x: Float,
    pub y: Float,
    pub z: Float,
    pub value: Float,
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub fn precompute_bicubic_coefficients(points: [Point2D; 4]) -> [Float; 16] {
    let mut lhs = Vec::<Float>::with_capacity(256);
    let mut rhs = Vec::<Float>::with_capacity(16);

    for point in points {
        let x = point.x;
        let y = point.y;

        lhs.append(&mut vec![
            1.0,
            y,
            y * y,
            y * y * y,
            x,
            x * y,
            x * y * y,
            x * y * y * y,
            x * x,
            x * x * y,
            x * x * y * y,
            x * x * y * y * y,
            x * x * x,
            x * x * x * y,
            x * x * x * y * y,
            x * x * x * y * y * y,
        ]);
        lhs.append(&mut vec![
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            y,
            y * y,
            y * y * y,
            2.0 * x,
            2.0 * x * y,
            2.0 * x * y * y,
            2.0 * x * y * y * y,
            3.0 * x * x,
            3.0 * x * x * y,
            3.0 * x * x * y * y,
            3.0 * x * x * y * y * y,
        ]);
        lhs.append(&mut vec![
            0.0,
            1.0,
            2.0 * y,
            3.0 * y * y,
            0.0,
            x,
            2.0 * x * y,
            3.0 * x * y * y,
            0.0,
            x * x,
            2.0 * x * x * y,
            3.0 * x * x * y * y,
            0.0,
            x * x * x,
            2.0 * x * x * x * y,
            3.0 * x * x * x * y * y,
        ]);
        lhs.append(&mut vec![
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            2.0 * y,
            3.0 * y * y,
            0.0,
            2.0 * x,
            4.0 * x * y,
            6.0 * x * y * y,
            0.0,
            3.0 * x * x,
            6.0 * x * x * y,
            9.0 * x * x * y * y,
        ]);

        rhs.push(point.value);
        rhs.push(point.dx);
        rhs.push(point.dy);
        rhs.push(point.dxdy);
    }

    let lhs = Matrix16::from_row_slice(&lhs);

    let rhs = Vector16::from_vec(rhs);

    let lhs = lhs.try_inverse().unwrap();
    let coeffs = lhs * rhs;

    coeffs.try_into().unwrap()
}
