
# RADIATE

RADIATE is a ray-tracing software written in Fortran. Ray-traced delays can be calculated for observations in the microwave as well as in the optical frequency range. The computation is based on numerical weather models and the output files contain several tropospheric parameters


## Software requirements

RADIATE is run from the Linux command line (Terminal). To compile the software, 'gfortran 5 compiler' (or newer version) needs to be installed.


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

The Fortran version of the VieVS ray-tracer is to be operated from the Linux command line (Terminal). The procedure is generally divided into two parts: the compilation and the actual ray-tracing. For compilation, the "gfortran 5 compiler" or higher must be installed. On most modern Linux systems this is installed by default. If not, it can be downloaded [here](https://gcc.gnu.org/wiki/GFortran). Apart from the command line, the ray-tracer can also be run from Windows using an IDE like e.g. Microsoft Visual Studio together with the Intel Fortran Compiler. This makes adapting the script much more comfortable. However, for operational purposes the use of gfortran is recommended due to higher computational speed.

In the following, there is a step-by-step description of how to create ray-traced delays.  
1. Before calculation, the Numerical Weather Models (NWM) in text format of all desired epochs must be stored in *DATA/GRIB/* without subdirectories. These text files have to strictly follow a formatting, which is described in the section "Required format of the NWM in text format" below.
2. Open a Terminal window and navigate to the directory *FUN_TEXT/*, where all required functions of the ray-tracer are stored. These functions read the text-file versions of the NWM.
3. As a first step of the ray-tracing, an executable file has to be created, that is, all scripts have to be compiled. This is handled through a batch file, in which the compilation command is placed twice. To execute the batch file, type `./creategf5`. This will produce the executable file *radiate*.
4. Now the actual ray-tracing is executed. There are several options to accomplish this, which are listed in the subsection // Input options //. It is possible to perform the ray-tracing either for single sessions, or for a number of sessions at once. All sessions are referenced through their respective azel filename:
5. Single sessions: Type ''./radiate azel_yyyymmddhh_UNI.txt stations_file.ell''. The specification of the session name and the station coordinates file is the minimum requirement, all other options will get their default values. If further options are to be specified, they have to be appended following a blank.
6. Multiple sessions: In order to evaluate multiple sessions at once, batch files have to be created. One such batch file is '' ./domore '', which evaluates all sessions for which azel files are available in // DATA/AZEL/ //. Data in subdirectories is not considered. Changes in the input options to the ray-tracer have to be done inside the batch file. 
7. The results of the ray-tracing are stored in // RESULTS/ //.




