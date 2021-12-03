//! This is a module for integration tests of the model,
//! but with access to private fields and methods.
//!
//! As multiple model methods are using the environment readed
//! from the GRIB files it would be tedious to write an environment
//! setup for each unit test. So this "super-unit-test" is a workaround
//! for that issue.

use std::path::Path;

use crate::model::configuration::Config;

#[test]
fn pressure_interpolation() {
    let cfg = Config::new_from_file(Path::new("./test-data/config.yaml"));
}
