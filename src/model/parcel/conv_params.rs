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
