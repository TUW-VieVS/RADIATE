
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

In the following, there are directions on how to download, install and use the ray-tracer.


