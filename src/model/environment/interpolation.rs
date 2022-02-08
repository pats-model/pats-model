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

//! Module containing interpolation methods.

use crate::Float;
use nalgebra::{SMatrix, SVector};

type Vector16 = SVector<Float, 16>;
type Matrix16 = SMatrix<Float, 16, 16>;
type Vector64 = SVector<Float, 64>;
type Matrix64 = SMatrix<Float, 64, 64>;

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
    pub value: Float,

    pub x: Float,
    pub y: Float,
    pub z: Float,

    pub dx: Float,
    pub dy: Float,
    pub dz: Float,

    pub dxdy: Float,
    pub dxdz: Float,
    pub dydz: Float,

    pub dxdydz: Float,
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

    coeffs.into()
}

pub fn precompute_tricubic_coefficients(points: [Point3D; 8]) -> [Float; 64] {
    let mut lhs = Vec::<Float>::with_capacity(4096);
    let mut rhs = Vec::<Float>::with_capacity(64);

    for point in points {
        let x = point.x;
        let y = point.y;
        let z = point.z;

        lhs.append(&mut vec![
            1.0,
            z,
            z * z,
            z * z * z,
            y,
            y * z,
            y * z * z,
            y * z * z * z,
            y * y,
            y * y * z,
            y * y * z * z,
            y * y * z * z * z,
            y * y * y,
            y * y * y * z,
            y * y * y * z * z,
            y * y * y * z * z * z,
            x,
            x * z,
            x * z * z,
            x * z * z * z,
            x * y,
            x * y * z,
            x * y * z * z,
            x * y * z * z * z,
            x * y * y,
            x * y * y * z,
            x * y * y * z * z,
            x * y * y * z * z * z,
            x * y * y * y,
            x * y * y * y * z,
            x * y * y * y * z * z,
            x * y * y * y * z * z * z,
            x * x,
            x * x * z,
            x * x * z * z,
            x * x * z * z * z,
            x * x * y,
            x * x * y * z,
            x * x * y * z * z,
            x * x * y * z * z * z,
            x * x * y * y,
            x * x * y * y * z,
            x * x * y * y * z * z,
            x * x * y * y * z * z * z,
            x * x * y * y * y,
            x * x * y * y * y * z,
            x * x * y * y * y * z * z,
            x * x * y * y * y * z * z * z,
            x * x * x,
            x * x * x * z,
            x * x * x * z * z,
            x * x * x * z * z * z,
            x * x * x * y,
            x * x * x * y * z,
            x * x * x * y * z * z,
            x * x * x * y * z * z * z,
            x * x * x * y * y,
            x * x * x * y * y * z,
            x * x * x * y * y * z * z,
            x * x * x * y * y * z * z * z,
            x * x * x * y * y * y,
            x * x * x * y * y * y * z,
            x * x * x * y * y * y * z * z,
            x * x * x * y * y * y * z * z * z,
        ]);
        lhs.append(&mut vec![
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            z,
            z * z,
            z * z * z,
            y,
            y * z,
            y * z * z,
            y * z * z * z,
            y * y,
            y * y * z,
            y * y * z * z,
            y * y * z * z * z,
            y * y * y,
            y * y * y * z,
            y * y * y * z * z,
            y * y * y * z * z * z,
            2.0 * x,
            2.0 * x * z,
            2.0 * x * z * z,
            2.0 * x * z * z * z,
            2.0 * x * y,
            2.0 * x * y * z,
            2.0 * x * y * z * z,
            2.0 * x * y * z * z * z,
            2.0 * x * y * y,
            2.0 * x * y * y * z,
            2.0 * x * y * y * z * z,
            2.0 * x * y * y * z * z * z,
            2.0 * x * y * y * y,
            2.0 * x * y * y * y * z,
            2.0 * x * y * y * y * z * z,
            2.0 * x * y * y * y * z * z * z,
            3.0 * x * x,
            3.0 * x * x * z,
            3.0 * x * x * z * z,
            3.0 * x * x * z * z * z,
            3.0 * x * x * y,
            3.0 * x * x * y * z,
            3.0 * x * x * y * z * z,
            3.0 * x * x * y * z * z * z,
            3.0 * x * x * y * y,
            3.0 * x * x * y * y * z,
            3.0 * x * x * y * y * z * z,
            3.0 * x * x * y * y * z * z * z,
            3.0 * x * x * y * y * y,
            3.0 * x * x * y * y * y * z,
            3.0 * x * x * y * y * y * z * z,
            3.0 * x * x * y * y * y * z * z * z,
        ]);
        lhs.append(&mut vec![
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            z,
            z * z,
            z * z * z,
            2.0 * y,
            2.0 * y * z,
            2.0 * y * z * z,
            2.0 * y * z * z * z,
            3.0 * y * y,
            3.0 * y * y * z,
            3.0 * y * y * z * z,
            3.0 * y * y * z * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            x,
            x * z,
            x * z * z,
            x * z * z * z,
            2.0 * x * y,
            2.0 * x * y * z,
            2.0 * x * y * z * z,
            2.0 * x * y * z * z * z,
            3.0 * x * y * y,
            3.0 * x * y * y * z,
            3.0 * x * y * y * z * z,
            3.0 * x * y * y * z * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            x * x,
            x * x * z,
            x * x * z * z,
            x * x * z * z * z,
            2.0 * x * x * y,
            2.0 * x * x * y * z,
            2.0 * x * x * y * z * z,
            2.0 * x * x * y * z * z * z,
            3.0 * x * x * y * y,
            3.0 * x * x * y * y * z,
            3.0 * x * x * y * y * z * z,
            3.0 * x * x * y * y * z * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            x * x * x,
            x * x * x * z,
            x * x * x * z * z,
            x * x * x * z * z * z,
            2.0 * x * x * x * y,
            2.0 * x * x * x * y * z,
            2.0 * x * x * x * y * z * z,
            2.0 * x * x * x * y * z * z * z,
            3.0 * x * x * x * y * y,
            3.0 * x * x * x * y * y * z,
            3.0 * x * x * x * y * y * z * z,
            3.0 * x * x * x * y * y * z * z * z,
        ]);
        lhs.append(&mut vec![
            0.0,
            1.0,
            2.0 * z,
            3.0 * z * z,
            0.0,
            y,
            2.0 * y * z,
            3.0 * y * z * z,
            0.0,
            y * y,
            2.0 * y * y * z,
            3.0 * y * y * z * z,
            0.0,
            y * y * y,
            2.0 * y * y * y * z,
            3.0 * y * y * y * z * z,
            0.0,
            x,
            2.0 * x * z,
            3.0 * x * z * z,
            0.0,
            x * y,
            2.0 * x * y * z,
            3.0 * x * y * z * z,
            0.0,
            x * y * y,
            2.0 * x * y * y * z,
            3.0 * x * y * y * z * z,
            0.0,
            x * y * y * y,
            2.0 * x * y * y * y * z,
            3.0 * x * y * y * y * z * z,
            0.0,
            x * x,
            2.0 * x * x * z,
            3.0 * x * x * z * z,
            0.0,
            x * x * y,
            2.0 * x * x * y * z,
            3.0 * x * x * y * z * z,
            0.0,
            x * x * y * y,
            2.0 * x * x * y * y * z,
            3.0 * x * x * y * y * z * z,
            0.0,
            x * x * y * y * y,
            2.0 * x * x * y * y * y * z,
            3.0 * x * x * y * y * y * z * z,
            0.0,
            x * x * x,
            2.0 * x * x * x * z,
            3.0 * x * x * x * z * z,
            0.0,
            x * x * x * y,
            2.0 * x * x * x * y * z,
            3.0 * x * x * x * y * z * z,
            0.0,
            x * x * x * y * y,
            2.0 * x * x * x * y * y * z,
            3.0 * x * x * x * y * y * z * z,
            0.0,
            x * x * x * y * y * y,
            2.0 * x * x * x * y * y * y * z,
            3.0 * x * x * x * y * y * y * z * z,
        ]);
        lhs.append(&mut vec![
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            z,
            z * z,
            z * z * z,
            2.0 * y,
            2.0 * y * z,
            2.0 * y * z * z,
            2.0 * y * z * z * z,
            3.0 * y * y,
            3.0 * y * y * z,
            3.0 * y * y * z * z,
            3.0 * y * y * z * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            2.0 * x,
            2.0 * x * z,
            2.0 * x * z * z,
            2.0 * x * z * z * z,
            4.0 * x * y,
            4.0 * x * y * z,
            4.0 * x * y * z * z,
            4.0 * x * y * z * z * z,
            6.0 * x * y * y,
            6.0 * x * y * y * z,
            6.0 * x * y * y * z * z,
            6.0 * x * y * y * z * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            3.0 * x * x,
            3.0 * x * x * z,
            3.0 * x * x * z * z,
            3.0 * x * x * z * z * z,
            6.0 * x * x * y,
            6.0 * x * x * y * z,
            6.0 * x * x * y * z * z,
            6.0 * x * x * y * z * z * z,
            9.0 * x * x * y * y,
            9.0 * x * x * y * y * z,
            9.0 * x * x * y * y * z * z,
            9.0 * x * x * y * y * z * z * z,
        ]);
        lhs.append(&mut vec![
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            2.0 * z,
            3.0 * z * z,
            0.0,
            y,
            2.0 * y * z,
            3.0 * y * z * z,
            0.0,
            y * y,
            2.0 * y * y * z,
            3.0 * y * y * z * z,
            0.0,
            y * y * y,
            2.0 * y * y * y * z,
            3.0 * y * y * y * z * z,
            0.0,
            2.0 * x,
            4.0 * x * z,
            6.0 * x * z * z,
            0.0,
            2.0 * x * y,
            4.0 * x * y * z,
            6.0 * x * y * z * z,
            0.0,
            2.0 * x * y * y,
            4.0 * x * y * y * z,
            6.0 * x * y * y * z * z,
            0.0,
            2.0 * x * y * y * y,
            4.0 * x * y * y * y * z,
            6.0 * x * y * y * y * z * z,
            0.0,
            3.0 * x * x,
            6.0 * x * x * z,
            9.0 * x * x * z * z,
            0.0,
            3.0 * x * x * y,
            6.0 * x * x * y * z,
            9.0 * x * x * y * z * z,
            0.0,
            3.0 * x * x * y * y,
            6.0 * x * x * y * y * z,
            9.0 * x * x * y * y * z * z,
            0.0,
            3.0 * x * x * y * y * y,
            6.0 * x * x * y * y * y * z,
            9.0 * x * x * y * y * y * z * z,
        ]);
        lhs.append(&mut vec![
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            2.0 * z,
            3.0 * z * z,
            0.0,
            2.0 * y,
            4.0 * y * z,
            6.0 * y * z * z,
            0.0,
            3.0 * y * y,
            6.0 * y * y * z,
            9.0 * y * y * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            x,
            2.0 * x * z,
            3.0 * x * z * z,
            0.0,
            2.0 * x * y,
            4.0 * x * y * z,
            6.0 * x * y * z * z,
            0.0,
            3.0 * x * y * y,
            6.0 * x * y * y * z,
            9.0 * x * y * y * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            x * x,
            2.0 * x * x * z,
            3.0 * x * x * z * z,
            0.0,
            2.0 * x * x * y,
            4.0 * x * x * y * z,
            6.0 * x * x * y * z * z,
            0.0,
            3.0 * x * x * y * y,
            6.0 * x * x * y * y * z,
            9.0 * x * x * y * y * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            x * x * x,
            2.0 * x * x * x * z,
            3.0 * x * x * x * z * z,
            0.0,
            2.0 * x * x * x * y,
            4.0 * x * x * x * y * z,
            6.0 * x * x * x * y * z * z,
            0.0,
            3.0 * x * x * x * y * y,
            6.0 * x * x * x * y * y * z,
            9.0 * x * x * x * y * y * z * z,
        ]);
        lhs.append(&mut vec![
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            2.0 * z,
            3.0 * z * z,
            0.0,
            2.0 * y,
            4.0 * y * z,
            6.0 * y * z * z,
            0.0,
            3.0 * y * y,
            6.0 * y * y * z,
            9.0 * y * y * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2.0 * x,
            4.0 * x * z,
            6.0 * x * z * z,
            0.0,
            4.0 * x * y,
            8.0 * x * y * z,
            12.0 * x * y * z * z,
            0.0,
            6.0 * x * y * y,
            12.0 * x * y * y * z,
            18.0 * x * y * y * z * z,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            3.0 * x * x,
            6.0 * x * x * z,
            9.0 * x * x * z * z,
            0.0,
            6.0 * x * x * y,
            12.0 * x * x * y * z,
            18.0 * x * x * y * z * z,
            0.0,
            9.0 * x * x * y * y,
            18.0 * x * x * y * y * z,
            27.0 * x * x * y * y * z * z,
        ]);

        rhs.push(point.value);
        rhs.push(point.dx);
        rhs.push(point.dy);
        rhs.push(point.dz);
        rhs.push(point.dxdy);
        rhs.push(point.dxdz);
        rhs.push(point.dydz);
        rhs.push(point.dxdydz);
    }

    let lhs = Matrix64::from_row_slice(&lhs);

    let rhs = Vector64::from_vec(rhs);

    let lhs = lhs.try_inverse().unwrap();
    let coeffs = lhs * rhs;

    coeffs.into()
}

pub fn interpolate_bicubic(x: Float, y: Float, coefficients: [Float; 16]) -> Float {
    let mut result = 0.0;

    for i in 0..=3 {
        for j in 0..=3 {
            result += coefficients[4 * i + j] * x.powi(i as i32) * y.powi(j as i32);
        }
    }

    result
}

pub fn interpolate_tricubic(x: Float, y: Float, z: Float, coefficients: [Float; 64]) -> Float {
    let mut result = 0.0;

    for i in 0..=3 {
        for j in 0..=3 {
            for k in 0..=3 {
                result += coefficients[16 * i + 4 * j + k]
                    * x.powi(i as i32)
                    * y.powi(j as i32)
                    * z.powi(k as i32);
            }
        }
    }

    result
}
