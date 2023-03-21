# XIOS testbed

Software to test [XIOS](https://forge.ipsl.jussieu.fr/ioserver/) output capabilities.

## Build code

The Software comes with four example build environments which can be adapted to suit other environments. Currently the example build environments are for two supercomputer centres [archer2](https://www.archer2.ac.uk) and [jasmin](https://jasmin.ac.uk).

There are four external libraries needed to build the code (plus all of their dependencies)

 - [NetCDF-FORTRAN](https://github.com/Unidata/netcdf-fortran) (compulsory)
 - [XIOS](https://forge.ipsl.jussieu.fr/ioserver/) (compulsory)
 - [Udunits-2](https://github.com/Unidata/UDUNITS-2) (optional)
 - [Eccodes](https://confluence.ecmwf.int/display/ECC/Releases) (optional)

When reading an input netCDF file the code will use the metadata to determine which axis type a given netCDF coordinate refers to (X,Y,Z or T). One method to do this is to examine the coordinates *units* attribute, this method can only be used when the code has been built with the udunits-2 library. This is the only usage of the udunits-2 library by the code and would only be needed if there is no other way of determining the axis type for the netCDF coordinate.

The eccodes library is used when the code has been setup to regrid the output data onto a Gaussian grid, otherwise this library is not needed.

Once the env_setup and Makefile files have been modified to suit the local environment, the xios test code can be compiled by running

```
source ./env_setup
make
```

One of the build environments uses ESDM (Earth System Data Middleware) on jasmin, using ESDM will require extensive work to port to a new environment. Details about ESDM can be found [here](https://github.com/ESiWACE/esdm), additionally a special version of [NetCDF-FORTRAN](https://github.com/ESiWACE/esdm-netcdf-4.6.2-old) is needed to support ESDM NetCDF output from XIOS.

## Run code

In each example build environment there are two job files used to run the code as a batch job on a computer resource. One job is to test writing ensemble data output using XIOS and the other is to test using XIOS to regrid the data. Some modification of these files will be needed to run on different computer systems.

There are three input files created in each of the job files, two Fortran namelist files needed by xios_test and a XML file needed by XIOS. One namelist file defines the input parameters for xios_test and the second defines the input netCDF file. There is one example netCDF file provided on a 192x145x17 grid, it should be possible to use other input NetCDF files with different grid sizes providing the horizontal grid is a regular latitude/longitude grid.

The final part of the job file runs xios_test and XIOS in MPMD mode, i.e. there are either one xios_test executable and one XIOS executable running for the regrid data job or for the ensemble job there are number of ensemble members lots of xios_test executables running plus one XIOS executable.

