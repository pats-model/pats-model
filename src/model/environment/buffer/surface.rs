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
//! surface data buffering.

use crate::{
    errors::{EnvironmentError, InputError},
    model::environment::{buffer::find_extent_edge_indices, Environment, LonLat},
    Float,
};
use eccodes::{
    codes_handle::{
        CodesHandle,
        KeyType::{FloatArray, Str},
        KeyedMessage,
        ProductKind::GRIB,
    },
    FallibleIterator,
};
use floccus::constants::G;
use log::debug;
use ndarray::{Array, Array2, Axis, azip, concatenate, s, stack};

impl Environment {
    /// Function to read surface data from GRIB input
    /// in extent covering domain and margins and buffer it.
    /// 
    /// Data from GRIB files can only be read in a whole
    /// GRIB domain, therefore it needs to be truncated.
    /// 
    /// Some useful surface variables are not provided by
    /// most NWP models, therefore they need to be computed.
    pub(in crate::model::environment) fn buffer_surface(
        &mut self,
        west_south: LonLat<Float>,
        east_north: LonLat<Float>,
    ) -> Result<(), EnvironmentError> {
        debug!(
            "Buffering surface fields in extent: N{} S{} E{} W{}",
            east_north.1, west_south.1, east_north.0, west_south.0
        );

        let distinct_lonlats = self.read_distinct_latlons()?;
        let (edge_lats, edge_lons) =
            find_extent_edge_indices(&distinct_lonlats, west_south, east_north);

        self.assign_raw_surfaces(edge_lons, edge_lats)?;
        self.compute_intermediate_surfaces();
        self.cast_lonlat_surface_coords(&distinct_lonlats, edge_lons, edge_lats);

        Ok(())
    }

    /// Reads variables on surface level from GRIB file
    /// and buffers them in domain + margins extent.
    fn assign_raw_surfaces(
        &mut self,
        edge_lons: (usize, usize),
        edge_lats: (usize, usize),
    ) -> Result<(), InputError> {
        let input_shape = self.read_input_shape()?;

        let geopotential = self.read_raw_surface("z", input_shape)?;
        self.surface.height =
            truncate_surface_to_extent(&geopotential, edge_lats, edge_lons).mapv(|v| v / G);

        let pressure = self.read_raw_surface("sp", input_shape)?;
        self.surface.pressure = truncate_surface_to_extent(&pressure, edge_lats, edge_lons);

        let temperature = self.read_raw_surface("2t", input_shape)?;
        self.surface.temperature = truncate_surface_to_extent(&temperature, edge_lats, edge_lons);

        let dewpoint = self.read_raw_surface("2d", input_shape)?;
        self.surface.dewpoint = truncate_surface_to_extent(&dewpoint, edge_lats, edge_lons);

        let u_wind = self.read_raw_surface("10u", input_shape)?;
        self.surface.u_wind = truncate_surface_to_extent(&u_wind, edge_lats, edge_lons);

        let v_wind = self.read_raw_surface("10v", input_shape)?;
        self.surface.v_wind = truncate_surface_to_extent(&v_wind, edge_lats, edge_lons);

        Ok(())
    }

    /// Reads all values in GRIB file at surface level
    /// of variable with given `short_name`.
    fn read_raw_surface(
        &self,
        short_name: &str,
        shape: (usize, usize),
    ) -> Result<Array2<Float>, InputError> {
        let mut data_level: Vec<KeyedMessage> = vec![];

        for file in &self.input.data_files {
            let handle = CodesHandle::new_from_file(file, GRIB)?;

            let mut data: Vec<KeyedMessage> = handle
                .filter(|msg| {
                    Ok(
                        msg.read_key("shortName")?.value == Str(short_name.to_string())
                            && msg.read_key("typeOfLevel")?.value == Str("surface".to_string()),
                    )
                })
                .collect()?;

            data_level.append(&mut data);
        }

        if data_level.is_empty() {
            return Err(InputError::DataNotSufficient(
                "Not enough variables in surface levels, check your input data",
            ));
        }

        let data_level = data_level[0].read_key("values")?.value;
        let data_level = if let FloatArray(v) = data_level {
            v
        } else {
            return Err(InputError::IncorrectKeyType("values"));
        };

        let result_data = Array2::from_shape_vec(shape, data_level)?;
        let result_data = result_data.mapv(|v| v as Float);

        Ok(result_data)
    }

    /// Computes and buffers additional surface data from
    /// values previously read from the GRIB file.
    fn compute_intermediate_surfaces(&mut self) {
        let mut mixing_ratio: Array2<Float> = Array2::zeros(self.surface.temperature.raw_dim());
        let mut sat_mixing_ratio: Array2<Float> = Array2::zeros(self.surface.temperature.raw_dim());
        let mut virtual_temperature: Array2<Float> =
            Array2::zeros(self.surface.temperature.raw_dim());

        azip!((r in &mut mixing_ratio, 
            &p in &self.surface.pressure, 
            &d in &self.surface.dewpoint) 
            *r = floccus::mixing_ratio::accuracy1(d, p)
            .expect("Error while computing mixing ratio: variable out of reasonable bounds"));

        self.surface.mixing_ratio = mixing_ratio;

        azip!((rs in &mut sat_mixing_ratio, 
            &p in &self.surface.pressure, 
            &t in &self.surface.temperature) 
            *rs = floccus::mixing_ratio::accuracy1(t, p)
            .expect("Error while computing saturation mixing ratio: variable out of reasonable bounds"));

        self.surface.sat_mixing_ratio = sat_mixing_ratio;

        azip!((tv in &mut virtual_temperature, 
            &t in &self.surface.temperature, 
            &r in &self.surface.mixing_ratio) 
            *tv = floccus::virtual_temperature::general1(t, r)
            .expect("Error while computing virtual temperature: variable out of reasonable bounds"));

        self.surface.virtual_temp = virtual_temperature;
    }

    /// Buffers longitudes and latitudes of surface data gridpoints.
    fn cast_lonlat_surface_coords(
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

        self.surface.lons = lons;
        self.surface.lats = lats;
    }
}

/// Truncates surface data array from GRIB file to
/// cover only the domain + margins extent.
fn truncate_surface_to_extent(
    raw_field: &Array2<Float>,
    edge_lats: (usize, usize),
    edge_lons: (usize, usize),
) -> Array2<Float> {
    //truncate in NS axis
    let truncated_field = raw_field.slice(s![.., edge_lats.0..=edge_lats.1]);

    if edge_lons.0 < edge_lons.1 {
        let truncated_field = truncated_field.slice(s![edge_lons.0..=edge_lons.1, ..]);
        return truncated_field.to_owned();
    }

    let left_half = truncated_field.slice(s![edge_lons.0.., ..]);
    let right_half = truncated_field.slice(s![..=edge_lons.1, ..]);
    let truncated_field = concatenate![Axis(0), left_half, right_half];
    truncated_field.to_owned()
}
