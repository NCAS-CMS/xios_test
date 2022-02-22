PROGRAM xios_test

USE, INTRINSIC :: ISO_FORTRAN_ENV

USE :: util_mod
USE :: timing_mod
USE :: input_grid_mod
USE :: read_netcdf_file_mod
USE :: mpp_mod
USE :: mpp_grid_mod
USE :: xios_mod
USE :: regrid_mod
USE :: ens_mod

USE :: xios
USE :: netcdf

IMPLICIT NONE

INTEGER, PARAMETER :: file_len = 1024
INTEGER :: it,nt,i

CHARACTER(LEN=file_len) :: filename
CHARACTER(LEN=nf90_max_name) :: varname_xios,varname_xios_copy

REAL(KIND=REAL64), ALLOCATABLE :: data(:,:,:)
REAL(KIND=REAL64), ALLOCATABLE :: local_data(:,:,:)

LOGICAL :: do_regrid,do_ens

CALL init_run()
CALL get_timing('After init_run')

CALL read_netcdf_file(filename,varname,data)
CALL get_timing('After read_netcdf_file')

CALL write_grid_info()

CALL decomp_grid_2d(nlon_in,nlat_in,lon_in,lat_in,lon_in_local,lat_in_local, &
                    start_in_x,size_in_x,start_in_y,size_in_y, &
                    pe_index_in_x,pe_index_in_y)
CALL get_timing('After decomp_grid input')

CALL scatter_data(data,local_data,start_in_x,size_in_x,start_in_y,size_in_y)
CALL get_timing('After scatter_data input')

IF (do_regrid) THEN
  CALL create_regrid()
  CALL get_timing('After create_regrid')
END IF

CALL xios_def()
CALL get_timing('After xios_def')

IF (do_regrid) THEN
  CALL xios_def_regrid()
  CALL get_timing('After xios_def_regrid')
END IF

IF (do_ens) THEN
  CALL xios_def_ens()
  CALL get_timing('After xios_def_ens')
END IF

CALL xios_clone_field(ncopies,do_regrid,do_ens,varname_xios)

CALL xios_close_context_definition()
CALL get_timing('After xios_close_context_definition')

DO it=1,nt
  IF (verbose > 0) THEN
    WRITE(message,*)'Timestep ',it
    CALL get_timing(message)
  END IF
  CALL xios_update_calendar(it)
  CALL xios_send_field(varname_xios,local_data)
  DO i=1,ncopies
    WRITE(varname_xios_copy,'(A,I0)') TRIM(varname_xios)//'_copy_',i
    CALL xios_send_field(varname_xios_copy,local_data)
  END DO
END DO
CALL get_timing('After time loop')

CALL xios_final()
CALL get_timing('After xios_final')

CALL cleanup()

IF (do_regrid) THEN
  CALL cleanup_regrid()
END IF

CALL get_timing('End of xios_test')

CLOSE(UNIT=outunit)

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE init_run()

IMPLICIT NONE

INTEGER :: ierr

CHARACTER(LEN=file_len) :: infile_data_base
CHARACTER(LEN=file_len) :: infile_data
CHARACTER(LEN=file_len) :: outfile

NAMELIST /input/ do_regrid,do_ens,infile_data_base, &
                 verbose,write_output,nprocx,nprocy,nt,ncopies
NAMELIST /input_data/ filename,varname,varname_xios

! Read input from namelist files

do_regrid = .FALSE.
do_ens = .FALSE.
infile_data_base = 'input_data'
verbose = 0
write_output = .FALSE.
ncopies = 0

OPEN(UNIT=nmlunit,FILE=infile)
READ(UNIT=nmlunit,NML=input,IOSTAT=ierr)
IF (ierr > 0) THEN
  WRITE(error_unit,*)'Error reading input namelist file, IOSTAT = ',ierr
  CALL abort()
END IF
CLOSE(UNIT=nmlunit)

! Initialise XIOS and MPI
IF (do_ens) THEN
  CALL get_ensemble_member_number()
END IF
CALL xios_init(do_ens,ensemble_member)

IF (do_ens) THEN
  CALL read_ens_nml()
  WRITE(outfile, '(A,I0,A,I0)') 'output_mem',ensemble_member,'_pe',mype
ELSE
  WRITE(outfile, '(A,I0)') 'output_pe',mype
END IF
OPEN(UNIT=outunit,FILE=outfile)

CALL get_timing('Start of xios_test', first=.TRUE.)
WRITE(UNIT=outunit,NML=input)
IF (do_ens) THEN
  WRITE(UNIT=outunit,NML=input_ens)
END IF

CALL check_input()

IF (do_ens) THEN
  WRITE(UNIT=outunit,FMT=*) 'Ensemble member number ',ensemble_member
  ! Each ensemble member has their own input grid file
  WRITE(infile_data, '(A,I0)') TRIM(infile_data_base)//'_mem',ensemble_member
ELSE
  infile_data = TRIM(infile_data_base)
END IF
WRITE(UNIT=outunit,FMT=*) 'Reading input grid info from file ',TRIM(infile_data)

OPEN(UNIT=nmlunit,FILE=infile_data)
READ(UNIT=nmlunit,NML=input_data,IOSTAT=ierr)
IF (ierr > 0) THEN
  WRITE(error_unit,*)'Error reading input_data namelist file, IOSTAT = ',ierr
  CALL abort()
END IF
CLOSE(UNIT=nmlunit)
WRITE(UNIT=outunit,NML=input_data)

END SUBROUTINE init_run

!-------------------------------------------------------------------------------

SUBROUTINE write_grid_info()

IMPLICIT NONE

INTEGER :: iy,iz

IF (mype == 0) THEN

  WRITE(outunit,*) 'Input lon/lat size = ',nlon_in,nlat_in
  IF (verbose > 0) THEN
    WRITE(outunit,*) 'Input longitude = ',lon_in
    WRITE(outunit,*) 'Input latitude = ',lat_in
  END IF
  WRITE(outunit,*) 'Input data ',TRIM(varname),' = '
  iz = 1
  iy = 1
  WRITE(outunit,*) data(1,iy,iz),data(2,iy,iz),data(3,iy,iz),' ... ', &
             data(nlon_in-2,iy,iz),data(nlon_in-1,iy,iz),data(nlon_in,iy,iz)
  iy = 2
  WRITE(outunit,*) data(1,iy,iz),data(2,iy,iz),data(3,iy,iz),' ... ', &
             data(nlon_in-2,iy,iz),data(nlon_in-1,iy,iz),data(nlon_in,iy,iz)
  iy = nlat_in-1
  WRITE(outunit,*) data(1,iy,iz),data(2,iy,iz),data(3,iy,iz),' ... ', &
             data(nlon_in-2,iy,iz),data(nlon_in-1,iy,iz),data(nlon_in,iy,iz)
  iy = nlat_in
  WRITE(outunit,*) data(1,iy,iz),data(2,iy,iz),data(3,iy,iz),' ... ', &
             data(nlon_in-2,iy,iz),data(nlon_in-1,iy,iz),data(nlon_in,iy,iz)
END IF

END SUBROUTINE write_grid_info

!-------------------------------------------------------------------------------

SUBROUTINE cleanup()

IMPLICIT NONE

IF (ALLOCATED(data)) DEALLOCATE(data)
IF (ALLOCATED(local_data)) DEALLOCATE(local_data)
IF (ALLOCATED(lon_in)) DEALLOCATE(lon_in)
IF (ALLOCATED(lat_in)) DEALLOCATE(lat_in)
IF (ALLOCATED(lev_in)) DEALLOCATE(lev_in)
IF (ALLOCATED(size_in_x)) DEALLOCATE(size_in_x)
IF (ALLOCATED(size_in_y)) DEALLOCATE(size_in_y)
IF (ALLOCATED(start_in_x)) DEALLOCATE(start_in_x)
IF (ALLOCATED(start_in_y)) DEALLOCATE(start_in_y)
IF (ALLOCATED(pe_index_in_x)) DEALLOCATE(pe_index_in_x)
IF (ALLOCATED(pe_index_in_y)) DEALLOCATE(pe_index_in_y)
IF (ALLOCATED(lon_in_local)) DEALLOCATE(lon_in_local)
IF (ALLOCATED(lat_in_local)) DEALLOCATE(lat_in_local)

END SUBROUTINE cleanup

END PROGRAM xios_test
