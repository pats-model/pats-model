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

//! Sub-module responsible for handling
//! pressure level data buffering.
use eccodes::{
    codes_handle::{
        CodesHandle,
        KeyType::{FloatArray, Int, Str},
        KeyedMessage,
        ProductKind::GRIB,
    },
    FallibleIterator,
};
use floccus::constants::G;
use log::debug;
use ndarray::{concatenate, s, stack, Array, Array2, Array3, Axis, Zip};

use crate::{
    errors::{EnvironmentError, InputError},
    model::environment::{buffer::find_extent_edge_indices, Environment, LonLat},
    Float,
};

impl Environment {
    /// Function to read pressure level data from GRIB input
    /// in extent covering domain and margins and buffer it.
    /// 
    /// Data from GRIB files can only be read in a whole
    /// GRIB domain, therefore it needs to be truncated.
    /// 
    /// Some useful variables are not provided by
    /// most NWP models, therefore they need to be computed.
    pub(in crate::model::environment) fn buffer_fields(
        &mut self,
        west_south: LonLat<Float>,
        east_north: LonLat<Float>,
    ) -> Result<(), EnvironmentError> {
        debug!(
            "Buffering input fields in extent: N{:.2} S{:.2} E{:.2} W{:.2}",
            east_north.1, west_south.1, east_north.0, west_south.0
        );

        let distinct_lonlats = self.read_distinct_latlons()?;
        let (edge_lats, edge_lons) =
            find_extent_edge_indices(&distinct_lonlats, west_south, east_north);

        self.assign_raw_fields(edge_lons, edge_lats)?;
        self.compute_intermediate_fields();
        self.cast_lonlat_fields_coords(&distinct_lonlats, edge_lons, edge_lats);

        Ok(())
    }

    /// Reads variables on pressure levels from GRIB file
    /// and buffers them in domain + margins extent.
    fn assign_raw_fields(
        &mut self,
        edge_lons: (usize, usize),
        edge_lats: (usize, usize),
    ) -> Result<(), InputError> {
        let input_shape = self.read_input_shape()?;

        self.fields.pressure = self.read_truncated_pressure(edge_lons, edge_lats);

        let geopotential = self.read_raw_field("z", input_shape)?;
        self.fields.height =
            truncate_field_to_extent(&geopotential, edge_lats, edge_lons).mapv(|v| v / G);

        let temperature = self.read_raw_field("t", input_shape)?;
        self.fields.temperature = truncate_field_to_extent(&temperature, edge_lats, edge_lons);

        let u_wind = self.read_raw_field("u", input_shape)?;
        self.fields.u_wind = truncate_field_to_extent(&u_wind, edge_lats, edge_lons);

        let v_wind = self.read_raw_field("v", input_shape)?;
        self.fields.v_wind = truncate_field_to_extent(&v_wind, edge_lats, edge_lons);

        let spec_humidity = self.read_raw_field("q", input_shape)?;
        self.fields.spec_humidity = truncate_field_to_extent(&spec_humidity, edge_lats, edge_lons);

        Ok(())
    }

    /// Reads all values in GRIB file at specified level type
    /// of variable with given `short_name` and collects them
    /// into a 3d array.
    /// 
    /// In GRIB files data on each level is stored as a separate
    /// message. Those need to be collected and only then can be
    /// converted into a 3d array.
    fn read_raw_field(
        &self,
        short_name: &str,
        shape: (usize, usize),
    ) -> Result<Array3<Float>, InputError> {
        let data_levels = self.read_raw_messages(short_name)?;
        let result_data = messages_to_array(data_levels, shape)?;

        Ok(result_data)
    }

    /// Filters and read all GRIB messages that contain
    /// variable with given `short_name` on specified level type.
    fn read_raw_messages(&self, short_name: &str) -> Result<Vec<KeyedMessage>, InputError> {
        let mut data_levels: Vec<KeyedMessage> = vec![];

        for file in &self.input.data_files {
            let handle = CodesHandle::new_from_file(file, GRIB)?;

            let mut data: Vec<KeyedMessage> = handle
                .filter(|msg| {
                    Ok(
                        msg.read_key("shortName")?.value == Str(short_name.to_string())
                            && msg.read_key("typeOfLevel")?.value
                                == Str(self.input.level_type.clone()),
                    )
                })
                .collect()?;

            data_levels.append(&mut data);
        }

        Ok(data_levels)
    }

    /// Creates a 3d array of pressure data of shape
    /// identical to other pressure level fields.
    /// 
    /// When data in GRIB file is provided on pressure levels
    /// the information about pressure at each level is only
    /// stored in message metadata. Thus, information what
    /// pressure levels are available has to be extracted
    /// and then casted to the 3d array expected by the
    /// [`accesser`](super::super::accesser).
    fn read_truncated_pressure(
        &self,
        edge_lons: (usize, usize),
        edge_lats: (usize, usize),
    ) -> Array3<Float> {
        let xy_shape = (
            (edge_lons.1 as isize - edge_lons.0 as isize).abs() as usize + 1,
            (edge_lats.1 as isize - edge_lats.0 as isize).abs() as usize + 1,
        );

        let mut pressure_levels = vec![];

        for level in &self.levels {
            let pressure_level = Array2::from_elem(xy_shape, *level);
            let pressure_level = pressure_level.mapv(|v| (v as Float) * 100.0);
            pressure_levels.push(pressure_level);
        }

        let mut pressure_views = vec![];

        for level in &pressure_levels {
            pressure_views.push(level.view());
        }

        let pressure_levels = ndarray::stack(Axis(0), pressure_views.as_slice()).unwrap();

        pressure_levels
    }

    /// Computes and buffers additional pressure level data from
    /// values previously read from the GRIB file.
    fn compute_intermediate_fields(&mut self) {
        let mut virtual_temperature: Array3<Float> =
            Array3::zeros(self.fields.temperature.raw_dim());

        Zip::from(&mut virtual_temperature)
            .and(&self.fields.temperature)
            .and(&self.fields.spec_humidity)
            .for_each(|tv, &t, &q| {
                *tv = floccus::virtual_temperature::general3(t, q).expect(
                    "Error while computing virtual temperature: variable out of reasonable bounds",
                );
            });

        self.fields.virtual_temp = virtual_temperature;
    }

    /// Buffers longitudes and latitudes of pressure level data gridpoints.
    fn cast_lonlat_fields_coords(
        &mut self,
        distinct_lonlats: &(Vec<Float>, Vec<Float>),
        edge_lons: (usize, usize),
        edge_lats: (usize, usize),
    ) {
        let lats = distinct_lonlats.1[edge_lats.0..=edge_lats.1].to_vec();
        let lons;

        if edge_lons.0 < edge_lons.1 {
            lons = distinct_lonlats.0[edge_lons.0..=edge_lons.1].to_vec();
        } else {
            let left_half = &distinct_lonlats.0[edge_lons.1..];
            let right_half = &distinct_lonlats.0[..=edge_lons.0];

            lons = [left_half, right_half].concat();
        }

        let lons = Array::from_vec(lons);
        let lats = Array::from_vec(lats);

        let lons_view = vec![lons.view(); lats.len()];
        let lats_view = vec![lats.view(); lons.len()];

        let lons = stack(Axis(1), lons_view.as_slice()).unwrap();
        let lats = stack(Axis(0), lats_view.as_slice()).unwrap();

        self.fields.lons = lons;
        self.fields.lats = lats;
    }
}

/// Collects data from GRIB messages on specified level type
/// into a 3d array,
fn messages_to_array(
    data_levels: Vec<KeyedMessage>,
    shape: (usize, usize),
) -> Result<Array3<Float>, InputError> {
    let mut sorted_data_levels = vec![];

    for msg in data_levels {
        let lvl_id = if let Int(id) = msg.read_key("level")?.value {
            id
        } else {
            return Err(InputError::IncorrectKeyType("level"));
        };

        let lvl_vals = if let FloatArray(vals) = msg.read_key("values")?.value {
            vals
        } else {
            return Err(InputError::IncorrectKeyType("values"));
        };

        let lvl_vals = Array2::from_shape_vec(shape, lvl_vals)?;

        sorted_data_levels.push((lvl_id, lvl_vals));
    }

    sorted_data_levels.sort_unstable_by_key(|k| k.0);
    sorted_data_levels.reverse();

    let sorted_data_levels: Vec<Array2<f64>> =
        sorted_data_levels.into_iter().map(|t| t.1).collect();

    let mut result_data = vec![];
    for lvl in &sorted_data_levels {
        result_data.push(lvl.view());
    }

    let result_data = ndarray::stack(Axis(0), result_data.as_slice())?;
    let result_data = result_data.mapv(|v| v as Float);

    Ok(result_data)
}

/// Truncates data on specified level type from GRIB file
/// to cover only the message + margins extent.
fn truncate_field_to_extent(
    raw_field: &Array3<Float>,
    edge_lats: (usize, usize),
    edge_lons: (usize, usize),
) -> Array3<Float> {
    //truncate in NS axis
    let truncated_field = raw_field.slice(s![.., .., edge_lats.0..=edge_lats.1]);

    if edge_lons.0 < edge_lons.1 {
        let truncated_field = truncated_field.slice(s![.., edge_lons.0..=edge_lons.1, ..]);
        return truncated_field.to_owned();
    }

    let left_half = truncated_field.slice(s![.., edge_lons.0.., ..]);
    let right_half = truncated_field.slice(s![.., ..=edge_lons.1, ..]);
    let truncated_field = concatenate![Axis(1), left_half, right_half];
    truncated_field.to_owned()
}
