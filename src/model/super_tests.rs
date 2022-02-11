//! This is a module for integration tests of the model,
//! but with access to private fields and methods.
//!
//! As multiple model methods are using the environment readed
//! from the GRIB files it would be tedious to write an environment
//! setup for each unit test. So this "super-unit-test" is a workaround
//! for that issue.

use super::configuration::Config;
use super::environment::{EnvFields, Environment};
use std::path::Path;

#[test]
fn pressure_interpolation() {
    let mut cfg = Config::new_from_file(Path::new("./test-data/config.yaml")).unwrap();
    let env = Environment::new(&mut cfg).unwrap();

    let (x, y) = env
        .projection
        .project(cfg.domain.ref_lon, cfg.domain.ref_lat);

    for z in (250..=10_000).step_by(1) {
        let v = env
            .get_field_value(x, y, z as f64, EnvFields::Pressure)
            .unwrap();

        println!("{:>5.1} {:>5.2}", z as f64, v);
    }
}
