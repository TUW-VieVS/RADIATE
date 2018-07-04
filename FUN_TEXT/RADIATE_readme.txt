Readme for program RADIATE

by Armin Hofmeister
 
version of 30.09.2015

-------------------------------------------------
-------------------------------------------------
I. Important requirements
-------------------------------------------------
-------------------------------------------------


-------------------------------------------------
1. Grib-file with numerical weather model data
-------------------------------------------------

a) The grib-file must contain only 1 Epoch and must be named according to the epoch yyyymmddhh, where y is for the year digits, m is for the month digits, d is for the day digits
   and hh is for the hour digits of the epoch at exactly hh:00:00. Currently the file extension is expected to be ".armin".
   
b) The contained data must cover a global grid starting at 90° latitude and 0° longitude with intervals of [90°:-90°] for latitude and [0°:360°[ longitude.
   The input of other intervals, e.g. [-180°,180°] for longitude are not possible.
   The starting value for latitude and longitude must be strictly adhered to, otherwise the undulation grid would not fit. An error message will appear if these requirements are
   not met.

c) The number of pressure levels in the file must be exactly 25.

d) The file must contain records of geopotential Z, temperature T and specific humidity Q at each pressure level for a global grid.

e) The horizontal resolution of the grib-file in latitude and longitude is in principle variable as long as equivalently resolved undulation grid-files are available.
   Currently resolutions of
		    1°x1°
		  0.5°x0.5°
		 0.25°x0.25°
		0.125°x0.125°
		  0.1°x0.1°
   are possible as for these undulation-files have been created.

   
-------------------------------------------------
2. General specifics
-------------------------------------------------

a) The station names contained in the station catalogue and AZEL-file must not contain a space inside the name, otherwise this will lead to an abortion of the program.


-------------------------------------------------
3. Fortran specifics
-------------------------------------------------

a) The session name must use the convention of the VLBI sessions, so the azel-file must be named azel_[14 characters].txt.
   So, after "azel_" the session name is placed with 14 characters. It is important to really use exactly 14 characters, otherwise the program fails due to erroneous filenames.
   Example for a session name like in VLBI: 11SEP15_N004
   
   
   
   
   
   
   