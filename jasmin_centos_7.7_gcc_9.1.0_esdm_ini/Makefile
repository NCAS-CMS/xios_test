.SUFFIXES : .f90 .F90

VPATH = ../source

LIBS = -lstdc++

ESDM_PATH = /home/users/jcole/esdm/install

XIOS_DIR = /home/users/jcole/software/xios_trunk_r1870/config_GCC_LINUX_ESDM_dev
XIOS_LIBS = -L $(XIOS_DIR)/lib -lxios

HDF5_DIR = /home/users/kunkel/install/hdf5/lib
NETCDF_DIR = $(ESDM_PATH)/lib
NETCDF_LIBS = -L $(NETCDF_DIR) -Wl,-rpath=$(NETCDF_DIR) -L $(HDF5_DIR) -Wl,-rpath=$(HDF5_DIR) -lnetcdff -lnetcdf -lhdf5_hl -lhdf5

UDUNITS_DIR = 
UDUNITS_LIBS = 

ECCODES_DIR = 
ECCODES_LIBS = 

FC=mpif90
F90C=mpif90
CC=mpicc
F90FLAGS = -I $(XIOS_DIR)/inc -I $(ESDM_PATH)/include
CPPFLAGS = -I $(ESDM_PATH)/include
LDFLAGS = $(XIOS_LIBS) $(NETCDF_LIBS) $(UDUNITS_LIBS) $(ECCODES_LIBS) $(LIBS)
OBJ = xios_test_c.o xios_test_mod.o regrid_test_mod.o ens_test_mod.o xios_test.o

.f90.o:
	$(F90C) $(F90FLAGS) -c $<

.F90.o:
	$(F90C) $(F90FLAGS) $(CPPFLAGS) -c $<

xios_test: $(OBJ)
	$(F90C) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o *.mod xios_test
