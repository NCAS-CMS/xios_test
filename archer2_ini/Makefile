.SUFFIXES : .f90 .F90

VPATH = ../source

LIBS = -lstdc++

XIOS_DIR = /work/n02/n02/jwc/XIOS/extract/xios_trunk_r2245
XIOS_LIBS = -L $(XIOS_DIR)/lib -lxios

NETCDF_DIR = 
NETCDF_LIBS = 

UDUNITS_DIR = 
UDUNITS_LIBS = 

#ECCODES_DIR = /work/y07/y07/umshared/lib/cce-8.5.8/eccodes/2.10.0
#ECCODES_LIBS = -L $(ECCODES_DIR)/lib -leccodes_f90 -leccodes
ECCODES_DIR =
ECCODES_LIBS =


F90C=ftn
CC=cc
F90FLAGS = -I $(XIOS_DIR)/inc -g
#F90FLAGS = -I $(XIOS_DIR)/inc -hfp0 -hflex_mp=strict
#CPPFLAGS = -DUSE_ECCODES -DUSE_UDUNITS
CPPFLAGS =
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
