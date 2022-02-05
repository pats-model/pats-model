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

use super::interpolation::Point2D;
use crate::Float;
use ndarray::{s, Array2, Zip};

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub fn compute_points_with_derivatives(
    surface: Array2<Float>,
    x: Array2<Float>,
    y: Array2<Float>,
) -> Array2<Point2D> {
    // this finite-difference computation uses average of deltaX, deltaY
    // on both sides on the stencil
    // y coordinates are against the matrix direction (making things even wierder)
    let dx = &surface.slice(s![2.., 1..-1])
        - &surface.slice(s![0..-2, 1..-1])
            / (&x.slice(s![2.., 1..-1]) - &x.slice(s![0..-2, 1..-1]));

    let dy = &surface.slice(s![1..-1, 0..-2])
        - &surface.slice(s![1..-1, 2..]) / (&y.slice(s![1..-1, 0..-2]) - &y.slice(s![1..-1, 2..]));

    let dxdy = &surface.slice(s![2.., 0..-2])
        - &surface.slice(s![2.., 2..])
        - &surface.slice(s![0..-2, 0..-2])
        + &surface.slice(s![0..-2, 2..])
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