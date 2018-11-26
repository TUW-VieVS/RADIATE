
# RADIATE

RADIATE is a ray-tracing software written in Fortran. Ray-traced delays can be calculated for observations in the microwave as well as in the optical frequency range. The computation is based on numerical weather models and the output files contain several tropospheric parameters


## Software requirements

RADIATE is run from the Linux command line (Terminal). To compile the software, the *gfortran 5 compiler* (or newer version) needs to be installed.


## Dependencies

RADIATE is an independent software package of VieVS.


## Installation

For installation the following steps are required:
* Create a directory, for example 'RADIATE'
* Clone the VieVS module RADIATE into the new directory using

```
$ mkdir RADIATE
$ cd RADIATE
$ git clone https://github.com/TUW-VieVS/RADIATE.git
```


## Use the VieVS ray-tracer ##

The Fortran version of the VieVS ray-tracer is to be operated from the Linux command line (Terminal). The procedure is generally divided into two parts: the compilation and the actual ray-tracing. For compilation, the *gfortran 5 compiler* or higher must be installed. On most modern Linux systems this is installed by default. If not, it can be downloaded [here](https://gcc.gnu.org/wiki/GFortran). Apart from the command line, the ray-tracer can also be run from Windows using an IDE like e.g. Microsoft Visual Studio together with the Intel Fortran Compiler. This makes adapting the script much more comfortable. However, for operational purposes the use of gfortran is recommended due to higher computational speed.

In the following, there is a step-by-step description of how to create ray-traced delays.  
1. Before calculation, the Numerical Weather Models (NWM) in text format of all desired epochs must be stored in *DATA/GRIB/* without subdirectories. These text files have to strictly follow a formatting, which is described in the section "Required format of the NWM text files" below. We are aware that the creation of these text format NWMs might be one of the major challenges of the usage of the ray-tracer.
2. Open a Terminal window and navigate to the directory *FUN_TEXT/*, where all required functions of the ray-tracer are stored. These functions read the text-file versions of the NWM.
3. As a first step of the ray-tracing, an executable file has to be created, that is, all scripts have to be compiled. This is handled through a batch file, in which the compilation command is placed twice. To execute the batch file, type `./creategf5`. This will produce the executable file *radiate*.
4. Now the actual ray-tracing is executed. There are several options to accomplish this, which are listed in the section "Input arguments to the ray-tracer". It is possible to perform the ray-tracing either for single sessions, or for a number of sessions at once. All sessions are referenced through their respective azel filename:
5. Single sessions: Type `./radiate azel_yyyymmddhh_UNI.txt stations_file.ell`. The specification of the session name and the station coordinates file is the minimum requirement, all other options will get their default values. If further options are to be specified, they have to be appended following a blank.
6. Multiple sessions: In order to evaluate multiple sessions at once, batch files have to be created. One such batch file is `./domore`, which evaluates all sessions for which azel files are available in *DATA/AZEL/*. Data in subdirectories is not considered. Changes in the input options to the ray-tracer have to be done inside the batch file. 
7. The results of the ray-tracing are stored in *RESULTS/*.


### Input arguments to the ray-tracer ###

Following the command which runs the executable file `./radiate`, a set of input arguments must be passed over the command line, the first two of which are mandatory. If no further input arguments are specified, the default settings of each are applied. The first argument must be the session name and the second argument must be the station coordinates file, while the order in which the optional input arguments are placed is negligible.

1. Name of azel-file of the session which shall be processed in the format *azel_yyyymmddhh_UNI.txt*, e.g. `azel_2014071200_UNI.txt`
2. Name of the respective station coordinates file:
   - `vlbi.ell`: containing the ellipsoidal coordinates of all VLBI stations (~240)
   - `slr.ell`: containing the ellipsoidal coordinates of SLR stations (~170)
   - `gnss.ell`: containing the ellipsoidal coordinates of most IGS stations (~460)
   - `doris.ell`: containing the ellipsoidal coordinates of all DORIS stations (~200)
   - `gridpoint_coord_5x5`: containing the ellipsoidal coordinates of all grid points of a global 5°x5° grid (2592)
   - `gridpoint_coord_1x1`: containing the ellipsoidal coordinates of all grid points of a global 1°x1° grid (64800)
    
- Method for the vertical interpolation of the meteorological data:
  - `-profilewise` (default)
  - `-gridwise`
- Ray-tracing method (approach for defining intersection points and delays):
  - `-pwl` (default)
  - `-ref_pwl` (only realized in the MATLAB version!)
  - `-Thayer` (only realized in the MATLAB version!)
- Wavelength of observations (choose between ray-tracing for microwave or for optical techniques):
  - `-microwave` (default)
  - `-optical` (implemented wavelength: 532nm)
- Time interpolation method (mode of assigning the observations to epochs):
  - `-one_epoch_per_obs`: each observation is assigned to only that grib-file epoch which it is closest to the observation
  - `-two_epochs_per_obs` (default): each observation is assigned to the two surrounding grib-file epochs, and is linearly interpolated from them. This approach is indeed more precise, however it demands double calculation time. Note: if an observation is exactly at a NWM epoch, the ray-tracing is still carried out for the neigboring NWM as well, although the second epolog entry will not alter the time interpolation result.
- Specify whether a .trp file shall be created and saved to *RESULTS/* in addition to the standardly produced .radiate file:
  - `-trp` (default)
  - `-notrp`
- Specify whether the error log shall be saved to *DATA/ERROR_LOG*:
  - `-errorlog` (default)
  - `-noerrorlog`
- Specify whether the produced session index- and epolog-files shall be cleaned up at the end:
  - `-cleanup`
  - `-nocleanup` (default)
- Specify how sessions are evaluated:
  - `-readAzel` (default): the input settings (stations, azimuths, elevations) for the ray-tracing of the respective session have to be stored in an .azel file, containing a line for every single observation. The .azel file must be stored in *DATA/AZEL/*. For all the observations contained in it, ray-traced delays will be produced.
  - `-createUniAzel`: in case ray-traced delays for uniformly distributed azimuths and constant elevations are desired, a further option is possible. The input settings (number of uniformly distributed azimuths, list of elevations) have to be specified in *DATA/INPUT/azel_spec.txt*. Thus, a virtual .azel file is created internally (without being saved as a .txt. file) for all stations from the selected station coordinates file. This means that no .azel file is required in*DATA/AZEL/*.
  
  
### Required format of the NWM text files ###

The VieVS ray-tracer requires the input text files containing the NWM information to follow a strict formatting. In the following, this format is described. In the directory //DATA/GRIB/sample//, there are some example NWMs in the required text format. Please mind that these files may be too large to be opened with standard text editors like Notepad, Wordpad or Notepad++! You may want to use a powerful text viewer such as [lister](https://www.ghisler.com/lister/) for this purpose.

The resolution of the Numerical Weather Models MUST be 1°x1° in the horizontal and 25 pressure levels in the vertical. The 3 quantities which are of interest are:
- geopotential height Z [m]
- specific humidity Q 
- temperature T [K]

##### Header Lines: #####
1. `181 360`: number of latitudes, number of longitudes (hardcoded!)
2. `90.0000000 -90.0000000`: first and last latitudes (hardcoded!)
3. `0.00000000 359.000000`: first and last longitudes (hardcoded!)
4. `1.00000000 1.00000000`: interval of longitudes, interval of latitudes (hardcoded!)
5. `25`: number of pressure levels (hardcoded!)
6. list of pressure levels in [hPa], starting with the top pressure level and ending with the bottom pressure level
7. `3`: number of parameters (hardcoded!)
8. `yyyy mm dd hh mm ss ii`: date of the NWM; *ii* represents the forecast step, but it can contain any number

##### Data lines: #####
In total 75 lines, each value in the format: value (number of pressure level, latitude, longitude), with pressure level #1 being the highest pressure level (in most cases 1 hPa).

![grib_txt file](https://user-images.githubusercontent.com/45286008/49023114-bb073400-f196-11e8-8788-42924bf7085f.PNG)



### How to bring ECMWF NWMs into text format ###

In general, all kinds of NWMs can be converted into the text format described in the section above. At TU Wien, we download NWM data from the European Center for Medium-range Weather Forecasts (ECMWF), which comes in so-called .grib format. However, as this .grib files are not publicly available, we cannot distribute them, which is why all the steps written above are necessary in order to perform the ray-tracing. We are aware that the creation of these text format NWMs might be one of the major challenges of the usage of the ray-tracer.

In case of ECMWF data, we create the .grib files online for the specified height levels, latitudes and longitudes and the meteorological parameters geopotential height, specific humidity and temperature. With the ECMWF-owned Fortran tool *ecCodes* we decipher the binary .grib file and use the internal function `grib_filter` to extract the data columns and write them into the text file using Bash. In case you wish to have a detailed description of how to decipher .grib files from ECMWF, please contact us by email.



### Troubleshooting ###

In particular the compilation part is prone to errors. In the following, there is a (incomplete) list of possible error sources, depending on whether the error happened during the compilation, or during the actual ray-tracing.

- Errors during the compilation:
  - Running the compilation batch file *creategf5* does not work, consequently no executable file is produced, but no error message is displayed in the command line: comment out the first line in *creategf5*, saying `scl enable devtoolset-4 bash`. Most likely you will have to execute this command alone in the command line. 
   - Running the compilation batch file *creategf5* does not work, and loads of errors are displayed in the command line, saying that a certain .mod file is missing: apparently there is any defective .mod file in the folder. Delete all .mod files by typing `rm *.mod` into the command line, and then execute *creategf5* at least twice. A set of further error messages will appear, but in the end the compilation should work without any problems.
- Errors during the ray-tracing:
  - In case a certain input argument to the ray-tracer is wrong or missing, the output of the ray-tracer to the command window is generally self-explaining.
  
  
  
  ### Important hints for compilation with gfortran ###

- A line must not be longer than 132 letters, otherwise it will be truncated. This can be prevented if the flag `-ffree-line-length-none` is used in the compilation following the `gfortran` command.
- Numeric values in types (= structures) which shall get the value 0, must be initialized with this value. In Microsoft VS (with Intel Fortran Compiler) this is not necessary, however in gfortran it is mandatory.
- Divisions by 0 do not yield NaN as in Microsoft VS, but an error.
- To all real or double precision numbers, *d0* MUST be appended.
- Logical comparisons must be *.eqv.* instead of *==*.



### Documented problems of the ray-tracer ###

- in very very rare cases (in 4 of >16000 VLBI sessions in total), but e.g. in session *azel_03JUL16XF.txt* at station CHICHI10, the ray-tracer produces an error. More precisely, in line 541 of module_getRayTrace2D_pwl_global.f90 at some point the location dependent elevation angle *theta* gets *NaN*, perhaps because for some reason the underlying refractivities change too rapidly within the height levels. So far it could not be found out why exactly this is happening, even less how to solve it.





