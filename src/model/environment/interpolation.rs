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

//! Module containing interpolation methods.

use crate::Float;
use nalgebra::{Matrix4, SMatrix, SVector, Vector4};

type Vector8 = SVector<Float, 8>;
type Matrix8 = SMatrix<Float, 8, 8>;

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
pub struct Point2D {
    pub x: Float,
    pub y: Float,
    pub value: Float,
}

#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
pub struct Point3D {
    pub x: Float,
    pub y: Float,
    pub z: Float,
    pub value: Float,
}

/// Function computing bilinear interpolation on 2D surface
/// using polynomial fit from 4 given points and
/// coordinates of interpolated point.
pub fn interpolate_bilinear(x: Float, y: Float, points: [Point2D; 4]) -> Float {
    let lhs = Matrix4::from_row_slice(&[
        1.0,
        points[0].x,
        points[0].y,
        points[0].x * points[0].y,
        1.0,
        points[1].x,
        points[1].y,
        points[1].x * points[1].y,
        1.0,
        points[2].x,
        points[2].y,
        points[2].x * points[2].y,
        1.0,
        points[3].x,
        points[3].y,
        points[3].x * points[3].y,
    ]);

    let rhs = Vector4::from_column_slice(&[
        points[0].value,
        points[1].value,
        points[2].value,
        points[3].value,
    ]);

    let lhs = lhs.try_inverse().unwrap();
    let coeffs = lhs * rhs;

    coeffs[0] + coeffs[1] * x + coeffs[2] * y + coeffs[3] * x * y
}

/// Function computing bilinear interpolation in 3D field
/// using polynomial fit from 8 given points and
/// coordinates of interpolated point.
pub fn interpolate_tilinear(x: Float, y: Float, z: Float, points: [Point3D; 8]) -> Float {
    let lhs = Matrix8::from_row_slice(&[
        1.0,
        points[0].x,
        points[0].y,
        points[0].z,
        points[0].x * points[0].y,
        points[0].x * points[0].z,
        points[0].y * points[0].z,
        points[0].x * points[0].y * points[0].z,
        1.0,
        points[1].x,
        points[1].y,
        points[1].z,
        points[1].x * points[1].y,
        points[1].x * points[1].z,
        points[1].y * points[1].z,
        points[1].x * points[1].y * points[1].z,
        1.0,
        points[2].x,
        points[2].y,
        points[2].z,
        points[2].x * points[2].y,
        points[2].x * points[2].z,
        points[2].y * points[2].z,
        points[2].x * points[2].y * points[2].z,
        1.0,
        points[3].x,
        points[3].y,
        points[3].z,
        points[3].x * points[3].y,
        points[3].x * points[3].z,
        points[3].y * points[3].z,
        points[3].x * points[3].y * points[3].z,
        1.0,
        points[4].x,
        points[4].y,
        points[4].z,
        points[4].x * points[4].y,
        points[4].x * points[4].z,
        points[4].y * points[4].z,
        points[4].x * points[4].y * points[4].z,
        1.0,
        points[5].x,
        points[5].y,
        points[5].z,
        points[5].x * points[5].y,
        points[5].x * points[5].z,
        points[5].y * points[5].z,
        points[5].x * points[5].y * points[5].z,
        1.0,
        points[6].x,
        points[6].y,
        points[6].z,
        points[6].x * points[6].y,
        points[6].x * points[6].z,
        points[6].y * points[6].z,
        points[6].x * points[6].y * points[6].z,
        1.0,
        points[7].x,
        points[7].y,
        points[7].z,
        points[7].x * points[7].y,
        points[7].x * points[7].z,
        points[7].y * points[7].z,
        points[7].x * points[7].y * points[7].z,
    ]);

    let rhs = Vector8::from_column_slice(&[
        points[0].value,
        points[1].value,
        points[2].value,
        points[3].value,
        points[4].value,
        points[5].value,
        points[6].value,
        points[7].value,
    ]);

    let lhs = lhs.try_inverse().unwrap();
    let coeffs = lhs * rhs;

    coeffs[0]
        + coeffs[1] * x
        + coeffs[2] * y
        + coeffs[3] * z
        + coeffs[4] * x * y
        + coeffs[5] * x * z
        + coeffs[6] * y * z
        + coeffs[7] * x * y * z
}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use crate::Float;

    use super::{interpolate_bilinear, interpolate_tilinear, Point2D, Point3D};

    #[test]
    fn bilinear() {
        let p1 = Point2D {
            x: 0.0,
            y: 0.0,
            value: 1.0,
        };

        let p2 = Point2D {
            x: 0.0,
            y: 1.0,
            value: 2.0,
        };

        let p3 = Point2D {
            x: 1.0,
            y: 0.0,
            value: 3.0,
        };

        let p4 = Point2D {
            x: 1.0,
            y: 1.0,
            value: 4.0,
        };

        let r = interpolate_bilinear(0.5, 0.5, [p1, p2, p3, p4]);

        assert_approx_eq!(Float, r, 2.5);
    }

    #[test]
    fn trilinear() {
        let p1 = Point3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            value: 1.0,
        };

        let p2 = Point3D {
            x: 0.0,
            y: 1.0,
            z: 0.0,
            value: 2.0,
        };

        let p3 = Point3D {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            value: 3.0,
        };

        let p4 = Point3D {
            x: 1.0,
            y: 1.0,
            z: 0.0,
            value: 4.0,
        };

        let p5 = Point3D {
            x: 0.0,
            y: 0.0,
            z: 1.0,
            value: 5.0,
        };

        let p6 = Point3D {
            x: 0.0,
            y: 1.0,
            z: 1.0,
            value: 6.0,
        };

        let p7 = Point3D {
            x: 1.0,
            y: 0.0,
            z: 1.0,
            value: 7.0,
        };

        let p8 = Point3D {
            x: 1.0,
            y: 1.0,
            z: 1.0,
            value: 8.0,
        };

        let r = interpolate_tilinear(0.5, 0.5, 0.5, [p1, p2, p3, p4, p5, p6, p7, p8]);

        assert_approx_eq!(Float, r, 4.5);
    }
}
