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

//! Module containing methods responsible for
//! buffering data from GRIB input.

mod fields;
mod surfaces;

use eccodes::{CodesHandle, FallibleIterator, KeyType::Str, KeyedMessage, ProductKind::GRIB};

use crate::{errors::InputError, model::configuration};

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub(super) fn collect_fields(input: &configuration::Input) -> Result<Vec<KeyedMessage>, InputError> {
    let mut data_levels: Vec<KeyedMessage> = vec![];

    for file in &input.data_files {
        let handle = CodesHandle::new_from_file(file, GRIB)?;

        let mut data: Vec<KeyedMessage> = handle
            .filter(|msg| {
                Ok(
                    msg.read_key("typeOfLevel")?.value == Str(input.level_type.clone())
                        && (msg.read_key("shortName")?.value == Str("z".to_string())
                            || msg.read_key("shortName")?.value == Str("q".to_string())
                            || msg.read_key("shortName")?.value == Str("t".to_string())
                            || msg.read_key("shortName")?.value == Str("u".to_string())
                            || msg.read_key("shortName")?.value == Str("v".to_string())),
                )
            })
            .collect()?;

        data_levels.append(&mut data);
    }

    if data_levels.is_empty() {
        return Err(InputError::DataNotSufficient(
            "Not enough variables on isobaric levels, check your input data",
        ));
    }

    Ok(data_levels)
}

/// (TODO: What it is)
///
/// (Why it is neccessary)
pub(super) fn collect_surfaces(input: &configuration::Input) -> Result<Vec<KeyedMessage>, InputError> {
    let mut data_levels: Vec<KeyedMessage> = vec![];

    for file in &input.data_files {
        let handle = CodesHandle::new_from_file(file, GRIB)?;

        let mut data: Vec<KeyedMessage> = handle
            .filter(|msg| {
                Ok(
                    msg.read_key("typeOfLevel")?.value == Str("surface".to_string())
                        && (msg.read_key("shortName")?.value == Str("10u".to_string())
                            || msg.read_key("shortName")?.value == Str("10v".to_string())
                            || msg.read_key("shortName")?.value == Str("2t".to_string())
                            || msg.read_key("shortName")?.value == Str("2d".to_string())
                            || msg.read_key("shortName")?.value == Str("sp".to_string())
                            || msg.read_key("shortName")?.value == Str("z".to_string())),
                )
            })
            .collect()?;

        data_levels.append(&mut data);
    }

    if data_levels.is_empty() {
        return Err(InputError::DataNotSufficient(
            "Not enough variables on surface levels, check your input data",
        ));
    }

    Ok(data_levels)
}
