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

//! (TODO: What it is)
//!
//! (Why it is neccessary)

use super::interpolation::{Point2D, Point3D};
use crate::Float;
use ndarray::{concatenate, s, Array2, Array3, Axis, Zip};

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub fn compute_2d_points(
    surface: &Array2<Float>,
    x: &Array2<Float>,
    y: &Array2<Float>,
) -> Array2<Point2D> {
    // this finite-difference computation uses average of deltaX, deltaY
    // on both sides on the stencil
    // y coordinates are against the matrix direction (making things even wierder)
    let dx = (&surface.slice(s![2.., 1..-1]) - &surface.slice(s![0..-2, 1..-1]))
        / (&x.slice(s![2.., 1..-1]) - &x.slice(s![0..-2, 1..-1]));

    let dy = (&surface.slice(s![1..-1, 0..-2]) - &surface.slice(s![1..-1, 2..]))
        / (&y.slice(s![1..-1, 0..-2]) - &y.slice(s![1..-1, 2..]));

    let dxdy = (&surface.slice(s![2.., 0..-2])
        - &surface.slice(s![2.., 2..])
        - &surface.slice(s![0..-2, 0..-2])
        + &surface.slice(s![0..-2, 2..]))
        / ((&x.slice(s![2.., 1..-1]) - &x.slice(s![0..-2, 1..-1]))
            * (&y.slice(s![1..-1, 0..-2]) - &y.slice(s![1..-1, 2..])));

    let mut points: Array2<Point2D> = Array2::default(dxdy.raw_dim());

    Zip::from(&mut points)
        .and(&surface.slice(s![1..-1, 1..-1]))
        .and(&x.slice(s![1..-1, 1..-1]))
        .and(&y.slice(s![1..-1, 1..-1]))
        .for_each(|p, &v, &x, &y| {
            p.value = v;
            p.x = x;
            p.y = y;
        });

    Zip::from(&mut points)
        .and(&dx)
        .and(&dy)
        .and(&dxdy)
        .for_each(|p, &dx, &dy, &dxdy| {
            p.dx = dx;
            p.dy = dy;
            p.dxdy = dxdy;
        });

    points
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub fn compute_3d_points(
    field: &Array3<Float>,
    x: &Array3<Float>,
    y: &Array3<Float>,
    z: &Array3<Float>,
) -> Array3<Point3D> {
    // this finite-difference computation uses average of deltaX, deltaY, deltaZ
    // on both sides on the stencil
    // y coordinates are against the matrix direction (making things even wierder)

    // dx
    let dx = (&field.slice(s![0..-1, 2.., 1..-1]) - &field.slice(s![0..-1, 0..-2, 1..-1]))
        / (&x.slice(s![0..-1, 2.., 1..-1]) - &x.slice(s![0..-1, 0..-2, 1..-1]));

    // dy
    let dy = (&field.slice(s![0..-1, 1..-1, 0..-2]) - &field.slice(s![0..-1, 1..-1, 2..]))
        / (&y.slice(s![0..-1, 1..-1, 0..-2]) - &y.slice(s![0..-1, 1..-1, 2..]));

    // dxdy
    let dxdy = (&field.slice(s![0..-1, 2.., 0..-2])
        - &field.slice(s![0..-1, 2.., 2..])
        - &field.slice(s![0..-1, 0..-2, 0..-2])
        + &field.slice(s![0..-1, 0..-2, 2..]))
        / ((&x.slice(s![0..-1, 2.., 1..-1]) - &x.slice(s![0..-1, 0..-2, 1..-1]))
            * (&y.slice(s![0..-1, 1..-1, 0..-2]) - &y.slice(s![0..-1, 1..-1, 2..])));

    // dz
    let dz = (&field.slice(s![2.., 1..-1, 1..-1]) - &field.slice(s![0..-2, 1..-1, 1..-1]))
        / (&z.slice(s![2.., 1..-1, 1..-1]) - &z.slice(s![0..-2, 1..-1, 1..-1]));

    let dz_bottom = (-3.0 * &field.slice(s![0, 1..-1, 1..-1])
        + 4.0 * &field.slice(s![1, 1..-1, 1..-1])
        - &field.slice(s![2, 1..-1, 1..-1]))
        / (&z.slice(s![2, 1..-1, 1..-1]) - &z.slice(s![0, 1..-1, 1..-1]));

    let dz_bottom = dz_bottom.insert_axis(Axis(0));
    let dz = concatenate![Axis(0), dz_bottom, dz];

    // dxdz
    let dxdz = (&field.slice(s![2.., 2.., 1..-1])
        - &field.slice(s![0..-2, 2.., 1..-1])
        - &field.slice(s![2.., 0..-2, 1..-1])
        + &field.slice(s![0..-2, 0..-2, 1..-1]))
        / ((&x.slice(s![1..-1, 2.., 1..-1]) - &x.slice(s![1..-1, 0..-2, 1..-1]))
            * (&z.slice(s![2.., 1..-1, 1..-1]) - &z.slice(s![0..-2, 1..-1, 1..-1])));

    let dxdz_bottom = (-3.0
        * (&field.slice(s![0, 2.., 1..-1]) - &field.slice(s![0, 0..-2, 1..-1]))
        + 4.0 * (&field.slice(s![1, 2.., 1..-1]) - &field.slice(s![1, 0..-2, 1..-1]))
        - (&field.slice(s![2, 2.., 1..-1]) - &field.slice(s![2, 0..-2, 1..-1])))
        / ((&x.slice(s![0, 2.., 1..-1]) - &x.slice(s![0, 0..-2, 1..-1]))
            * (&z.slice(s![2, 1..-1, 1..-1]) - &z.slice(s![0, 1..-1, 1..-1])));

    let dxdz_bottom = dxdz_bottom.insert_axis(Axis(0));
    let dxdz = concatenate![Axis(0), dxdz_bottom, dxdz];

    // dydz
    let dydz = (&field.slice(s![2.., 1..-1, 0..-2])
        - &field.slice(s![0..-2, 1..-1, 0..-2])
        - &field.slice(s![2.., 1..-1, 2..])
        + &field.slice(s![0..-2, 1..-1, 2..]))
        / ((&y.slice(s![1..-1, 1..-1, 0..-2]) - &y.slice(s![1..-1, 1..-1, 2..]))
            * (&z.slice(s![2.., 1..-1, 1..-1]) - &z.slice(s![0..-2, 1..-1, 1..-1])));

    let dydz_bottom = (-3.0
        * (&field.slice(s![0, 1..-1, 0..-2]) - &field.slice(s![0, 1..-1, 2..]))
        + 4.0 * (&field.slice(s![1, 1..-1, 0..-2]) - &field.slice(s![1, 1..-1, 2..]))
        - (&field.slice(s![2, 1..-1, 0..-2]) - &field.slice(s![2, 1..-1, 2..])))
        / ((&y.slice(s![0, 1..-1, 0..-2]) - &y.slice(s![0, 1..-1, 2..]))
            * (&z.slice(s![2, 1..-1, 1..-1]) - &z.slice(s![0, 1..-1, 1..-1])));

    let dydz_bottom = dydz_bottom.insert_axis(Axis(0));
    let dydz = concatenate![Axis(0), dydz_bottom, dydz];

    // dxdydz
    let dxdydz = (&field.slice(s![2.., 2.., 0..-2])
        - &field.slice(s![2.., 2.., 2..])
        - &field.slice(s![2.., 0..-2, 0..-2])
        + &field.slice(s![2.., 0..-2, 2..])
        - &field.slice(s![0..-2, 2.., 0..-2])
        + &field.slice(s![0..-2, 2.., 2..])
        + &field.slice(s![0..-2, 0..-2, 0..-2])
        - &field.slice(s![0..-2, 0..-2, 2..]))
        / ((&x.slice(s![1..-1, 2.., 1..-1]) - &x.slice(s![1..-1, 0..-2, 1..-1]))
            * (&y.slice(s![1..-1, 1..-1, 0..-2]) - &y.slice(s![1..-1, 1..-1, 2..]))
            * (&z.slice(s![2.., 1..-1, 1..-1]) - &z.slice(s![0..-2, 1..-1, 1..-1]))); //

    let dxdydz_bottom = ((-3.0
        * ((&field.slice(s![0, 2.., 0..-2]))
            - (&field.slice(s![0, 2.., 2..]))
            - (&field.slice(s![0, 0..-2, 0..-2]))
            + (&field.slice(s![0, 0..-2, 2..]))))
        + (4.0
            * ((&field.slice(s![1, 2.., 0..-2]))
                - (&field.slice(s![1, 2.., 2..]))
                - (&field.slice(s![1, 0..-2, 0..-2]))
                + (&field.slice(s![1, 0..-2, 2..]))))
        - ((&field.slice(s![2, 2.., 0..-2]))
            - (&field.slice(s![2, 2.., 2..]))
            - (&field.slice(s![2, 0..-2, 0..-2]))
            + (&field.slice(s![2, 0..-2, 2..]))))
        / ((&x.slice(s![0, 2.., 1..-1]) - &x.slice(s![0, 0..-2, 1..-1]))
            * (&y.slice(s![0, 1..-1, 0..-2]) - &y.slice(s![0, 1..-1, 2..]))
            * (&z.slice(s![2, 1..-1, 1..-1]) - &z.slice(s![0, 1..-1, 1..-1])));

    let dxdydz_bottom = dxdydz_bottom.insert_axis(Axis(0));
    let dxdydz = concatenate![Axis(0), dxdydz_bottom, dxdydz];

    // collect points and derivatives into one array

    let mut points: Array3<Point3D> = Array3::default(dxdydz.raw_dim());

    Zip::from(&mut points)
        .and(&field.slice(s![0..-1, 1..-1, 1..-1]))
        .and(&x.slice(s![0..-1, 1..-1, 1..-1]))
        .and(&y.slice(s![0..-1, 1..-1, 1..-1]))
        .and(&z.slice(s![0..-1, 1..-1, 1..-1]))
        .for_each(|p, &v, &x, &y, &z| {
            p.value = v;
            p.x = x;
            p.y = y;
            p.z = z;
        });

    Zip::from(&mut points)
        .and(&dx)
        .and(&dy)
        .and(&dz)
        .for_each(|p, &dx, &dy, &dz| {
            p.dx = dx;
            p.dy = dy;
            p.dz = dz;
        });

    Zip::from(&mut points)
        .and(&dxdy)
        .and(&dxdz)
        .and(&dydz)
        .and(&dxdydz)
        .for_each(|p, &dxdy, &dxdz, &dydz, &dxdydz| {
            p.dxdy = dxdy;
            p.dxdz = dxdz;
            p.dydz = dydz;
            p.dxdydz = dxdydz;
        });

    points
}
