# Reading GRIB input

Values of GRIB data are interpolated using pseudo-tricubic method.
Initially, trilinear interpolation was used, but lack of continuous derivative of interpolated
values was suboptimal. The used method is pseudo-tricubic as it does not interpolates values
using one three-dimensional function but using a sequence of three one-dimensional functions.
The result is effectively the same as real tricubic, but much easier to implement and of
comparable performance.
