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

//! Module containg methods for conducting
//! binary search (bisection) of elements closests
//! to searched values in datasets.

use crate::errors::SearchError;

/// Core bisection function, simply an implementation
/// of binary search algorithm adapted to searching values
/// in-between the set items.
/// 
/// Alternatively, `binary_search()` function for slice type could be used,
/// but this function is highly customised to the model needs and there are no
/// apparent advantages of using built-in binary_search() over custom one.
fn binary_search<T: PartialOrd>(array: &[T], x: &T) -> Result<usize, SearchError> {
    if array.is_empty() {
        return Err(SearchError::EmptyArray);
    }

    if x < array.first().unwrap() && x < array.last().unwrap()
        || x > array.first().unwrap() && x > array.last().unwrap()
    {
        return Err(SearchError::OutOfBounds);
    }

    let mut lo = 0;
    let mut hi = array.len() - 1;

    // if the array is sorted descendingly we use a function with reversed signs
    if array.first().unwrap() < array.last().unwrap() {
        while lo < hi {
            let mid = (lo + hi) / 2;

            if array[mid] >= *x {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
    } else {
        while lo < hi {
            let mid = (lo + hi) / 2;

            if array[mid] <= *x {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
    }

    Ok(lo)
}

/// Convienience public method to find a closest value
/// to requested to the left of the searched item.
pub fn find_left_closest<T: PartialOrd>(array: &[T], x: &T) -> Result<usize, SearchError> {
    let found_index = binary_search(array, x)?;

    if array.first().unwrap() < array.last().unwrap() {
        if array[found_index] <= *x {
            Ok(found_index)
        } else {
            Ok(found_index - 1)
        }
    } else if array[found_index] >= *x {
        Ok(found_index)
    } else {
        Ok(found_index - 1)
    }
}

/// Convienience public method to find a closest value
/// to requested to the right of the searched item.
pub fn find_right_closest<T: PartialOrd>(array: &[T], x: &T) -> Result<usize, SearchError> {
    let found_index = binary_search(array, x)?;

    if array.first().unwrap() < array.last().unwrap() {
        if array[found_index] >= *x {
            Ok(found_index)
        } else {
            Ok(found_index - 1)
        }
    } else if array[found_index] <= *x {
        Ok(found_index)
    } else {
        Ok(found_index - 1)
    }
}
