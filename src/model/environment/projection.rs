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

//! Module with methods to do computations
//! of geographical projection used by the model.
//! Closely follows algorithms and instructions in:
//! <https://pubs.er.usgs.gov/publication/pp1395>

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
use float_cmp::approx_eq;

use crate::constants::{WGS84_A, WGS84_E};
use crate::{errors::ProjectionError, Float};

/// Front-facing struct of Lambert Conformal Conic projection.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
pub struct LambertConicConformal {
    lambda_0: Float,
    n: Float,
    big_f: Float,
    rho_0: Float,
}

impl LambertConicConformal {
    /// LCC projection constructor from reference longitude
    /// and two standard parallels.
    /// Defaults the reference latitude to 0.0
    pub fn new(lon_0: Float, lat_1: Float, lat_2: Float) -> Result<Self, ProjectionError> {
        if approx_eq!(Float, lat_1, lat_2) {
            return Err(ProjectionError::IncorrectParams(
                "standard parallels cannot be equal",
            ));
        }

        if !(-180.0..180.0).contains(&lon_0) {
            return Err(ProjectionError::IncorrectParams("longitude out of bounds"));
        }

        if !(-90.0..90.0).contains(&lat_1) || !(-90.0..90.0).contains(&lat_2) {
            return Err(ProjectionError::IncorrectParams("latitude out of bounds"));
        }

        if !lon_0.is_finite() || !lat_1.is_finite() || !lat_2.is_finite() {
            return Err(ProjectionError::IncorrectParams(
                "one of params is not finite",
            ));
        }

        let lat_0: Float = 0.0;

        let phi_0 = lat_0.to_radians();
        let phi_1 = lat_1.to_radians();
        let phi_2 = lat_2.to_radians();

        let t_0 = t(phi_0);
        let t_1 = t(phi_1);
        let t_2 = t(phi_2);
        let m_1 = m(phi_1);
        let m_2 = m(phi_2);

        let n = n(m_1, m_2, t_1, t_2);
        let big_f = big_f(m_1, n, t_1);
        let rho_0 = rho(big_f, t_0, n);

        Ok(LambertConicConformal {
            lambda_0: lon_0.to_radians(),
            n,
            big_f,
            rho_0,
        })
    }

    /// Function to project geographic coordinates
    /// on WGS84 ellipsoid to cartographic coordinates
    /// with previously specified LCC projection.
    pub fn project(&self, lon: Float, lat: Float) -> (Float, Float) {
        let phi = lat.to_radians();
        let lambda = lon.to_radians();

        let t = t(phi);
        let theta = self.n * (lambda - self.lambda_0);
        let rho = rho(self.big_f, t, self.n);

        let x = rho * theta.sin();
        let y = self.rho_0 - rho * theta.cos();

        (x, y)
    }

    /// Function to inversly project cartographic coordinates
    /// on specified LCC projection to geographic coordinates
    /// on WGS84 ellipsoid.
    pub fn inverse_project(&self, x: Float, y: Float) -> (Float, Float) {
        let rho = (self.n.signum()) * (x.powi(2) + (self.rho_0 - y).powi(2)).sqrt();

        let theta;
        {
            // adjusting signs locally for theta
            let sign = self.n.signum();
            let x = x * sign;
            let y = y * sign;
            let rho_0 = self.rho_0 * sign;
            theta = (x / (rho_0 - y)).atan();
        }

        let t = (rho / (WGS84_A * self.big_f)).powf(1.0 / self.n);

        let lambda = (theta / self.n) + self.lambda_0;
        let phi = phi_for_inverse(t);

        (lambda.to_degrees(), phi.to_degrees())
    }
}

fn t(phi: Float) -> Float {
    ((FRAC_PI_4 - 0.5 * phi).tan())
        / (((1.0 - WGS84_E * phi.sin()) / (1.0 + WGS84_E * phi.sin())).powf(WGS84_E / 2.0))
}

fn m(phi: Float) -> Float {
    phi.cos() / (1.0 - (WGS84_E.powi(2) * (phi.sin()).powi(2))).sqrt()
}

fn n(m_1: Float, m_2: Float, t_1: Float, t_2: Float) -> Float {
    (m_1.ln() - m_2.ln()) / (t_1.ln() - t_2.ln())
}

fn big_f(m_1: Float, n: Float, t_1: Float) -> Float {
    m_1 / (n * t_1.powf(n))
}

fn rho(big_f: Float, t: Float, n: Float) -> Float {
    WGS84_A * big_f * t.powf(n)
}

/// To compute the phi for inverse projection
/// truncated infinite series is used with
/// optimisations for reducing trigonometric
/// functions calls.
fn phi_for_inverse(t: Float) -> Float {
    let chi = FRAC_PI_2 - 2.0 * t.atan();

    let big_a = (WGS84_E.powi(2) / 2.0)
        + 5.0 * (WGS84_E.powi(4) / 24.0)
        + (WGS84_E.powi(6) / 12.0)
        + 13.0 * (WGS84_E.powi(8) / 360.0);

    let big_b = 7.0 * (WGS84_E.powi(4) / 48.0)
        + 29.0 * (WGS84_E.powi(6) / 240.0)
        + 811.0 * (WGS84_E.powi(8) / 11520.0);

    let big_c = 7.0 * (WGS84_E.powi(6) / 120.0) + 81.0 * (WGS84_E.powi(8) / 1120.0);

    let big_d = 4279.0 * (WGS84_E.powi(8) / 161_280.0);

    let a_prime = big_a - big_c;
    let b_prime = 2.0 * big_b - 4.0 * big_d;
    let c_prime = 4.0 * big_c;
    let d_prime = 8.0 * big_d;

    let sin_2chi = (2.0 * chi).sin();
    let cos_2chi = (2.0 * chi).cos();

    chi + (sin_2chi
        * (a_prime + (cos_2chi * (b_prime + (cos_2chi * (c_prime + (d_prime * cos_2chi)))))))
}

#[cfg(test)]
mod tests {
    use super::LambertConicConformal;

    #[test]
    fn project() {
        let proj = LambertConicConformal::new(18.0, 30.0, 60.0).unwrap();

        let (lon_0, lat_0) = (18.58973722443749, 54.41412855026378);

        let (x, y) = proj.project(lon_0, lat_0);
        let (lon, lat) = proj.inverse_project(x, y);
        let (xdiff, ydiff) = (lon - lon_0, lat - lat_0);

        assert!(xdiff < 0.000001);
        assert!(ydiff < 0.000001);
    }
}
