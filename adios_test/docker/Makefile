.SUFFIXES : .f90 .F90

VPATH = ../source

LIBS = 

ADIOS_CFLAGS = $(shell adios2-config --c-flags)
ADIOS_FFLAGS = $(shell adios2-config --fortran-flags)
ADIOS_LIBS = $(shell adios2-config --fortran-libs)

NETCDF_CFLAGS = $(shell nc-config --cflags)
NETCDF_FFLAGS = $(shell nc-config --fflags)
NETCDF_LIBS = $(shell nc-config --flibs)

UDUNITS_CFLAGS = 
UDUNITS_FFLAGS = 
UDUNITS_LIBS = 

ECCODES_CFLAGS = 
ECCODES_FFLAGS = 
ECCODES_LIBS = 

FC=mpif90
F90C=mpif90
CC=mpicc
CFLAGS = $(ADIOS_CFLAGS) $(NETCDF_CFLAGS) $(UDUNITS_CFLAGS) $(ECCODES_CFLAGS)
F90FLAGS = $(ADIOS_FFLAGS) $(NETCDF_FFLAGS) $(UDUNITS_FFLAGS) $(ECCODES_FFLAGS)
CPPFLAGS = 
LDFLAGS = $(ADIOS_LIBS) $(NETCDF_LIBS) $(UDUNITS_LIBS) $(ECCODES_LIBS) $(LIBS)
OBJ = adios_test_c.o adios_test_mod.o regrid_test_mod.o ens_test_mod.o adios_test.o

.f90.o:
	$(F90C) $(F90FLAGS) -c $<

.F90.o:
	$(F90C) $(F90FLAGS) $(CPPFLAGS) -c $<

adios_test: $(OBJ)
	$(F90C) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o *.mod *.smod adios_test
