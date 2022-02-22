
!===============================================================================

MODULE util_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: infile = 'xios_test.in'
INTEGER, PARAMETER :: outunit = 7
INTEGER, PARAMETER :: nmlunit = 8
INTEGER :: verbose
INTEGER :: ncopies
LOGICAL :: write_output

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE upcase(string)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(INOUT) :: string

INTEGER :: i,ilen,iaup,ialo,idum

ilen = LEN(string)
iaup = ICHAR('A')
ialo = ICHAR('a')
DO i = 1, ilen
  idum = ICHAR(string(i:i))
  IF (idum .ge. ialo .and. idum .le. ialo + 25) THEN
     string(i:i) = CHAR(idum - ialo + iaup)
  END IF
END DO

END SUBROUTINE upcase

!-------------------------------------------------------------------------------

SUBROUTINE locase(string)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(INOUT) :: string

INTEGER :: i,ilen,iaup,ialo,idum

ilen = LEN(string)
iaup = ICHAR('A')
ialo = ICHAR('a')
DO i = 1, ilen
  idum = ICHAR(string(i:i))
  IF (idum .ge. iaup .and. idum .le. iaup + 25) THEN
     string(i:i) = CHAR(idum - iaup + ialo)
  END IF
END DO

END SUBROUTINE locase

END MODULE util_mod

!===============================================================================

MODULE fort2c_udunits

#ifdef USE_UDUNITS

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR,C_PTR,C_INT,C_NULL_CHAR,C_NULL_PTR, &
                                       C_LOC,C_ASSOCIATED

IMPLICIT NONE

TYPE(C_PTR) :: utsystem

INTERFACE
  TYPE(C_PTR) FUNCTION read_xml_c(path) BIND (c, NAME="read_xml")

  IMPORT :: C_PTR

  TYPE(C_PTR), VALUE :: path

  END FUNCTION
END INTERFACE

INTERFACE
  INTEGER(C_INT) FUNCTION are_convertible(utsystem,unit1,unit2) BIND (c)

  IMPORT :: C_INT,C_CHAR,C_PTR

  TYPE(C_PTR), VALUE :: utsystem
  CHARACTER(LEN=1,KIND=C_CHAR) :: unit1(*),unit2(*)

  END FUNCTION
END INTERFACE

INTERFACE
  SUBROUTINE ut_free_system(utsystem) BIND (c)

  IMPORT :: C_PTR

  TYPE(C_PTR), VALUE :: utsystem

  END SUBROUTINE
END INTERFACE

CONTAINS

!-------------------------------------------------------------------------------

INTEGER FUNCTION read_xml(path)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: path

CHARACTER(LEN=1,KIND=C_CHAR), ALLOCATABLE, TARGET :: c_path(:)
TYPE(C_PTR) :: path_ptr

read_xml = 0

IF (PRESENT(path)) THEN
  ALLOCATE(c_path(LEN_TRIM(PATH)+1))
  c_path = TRIM(path)//C_NULL_CHAR
  path_ptr = C_LOC(c_path)
  utsystem = read_xml_c(path_ptr)
  DEALLOCATE(c_path)
ELSE
  utsystem = read_xml_c(C_NULL_PTR)
END IF

IF (.NOT. C_ASSOCIATED(utsystem)) read_xml = 1

END FUNCTION read_xml

!-------------------------------------------------------------------------------

SUBROUTINE free_system()

IMPLICIT NONE

CALL ut_free_system(utsystem)

END SUBROUTINE free_system

!-------------------------------------------------------------------------------

LOGICAL FUNCTION ispressure(unit)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: unit

CHARACTER(LEN=LEN_TRIM(unit)+1,KIND=C_CHAR) :: c_unit
CHARACTER(LEN=*,KIND=C_CHAR), PARAMETER :: pa_unit = "Pa"//C_NULL_CHAR

c_unit = TRIM(unit)//C_NULL_CHAR

ispressure = are_convertible(utsystem,c_unit,pa_unit) /= 0

END FUNCTION ispressure

!-------------------------------------------------------------------------------

LOGICAL FUNCTION istime(unit)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: unit

CHARACTER(LEN=LEN_TRIM(unit)+1,KIND=C_CHAR) :: c_unit
CHARACTER(LEN=*,KIND=C_CHAR), PARAMETER :: sec_unit = "Second"//C_NULL_CHAR

c_unit = TRIM(unit)//C_NULL_CHAR

istime = are_convertible(utsystem,c_unit,sec_unit) /= 0

END FUNCTION istime

#else

IMPLICIT NONE

CONTAINS

!-------------------------------------------------------------------------------

INTEGER FUNCTION read_xml(path)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: path

read_xml = 1

END FUNCTION read_xml

!-------------------------------------------------------------------------------

SUBROUTINE free_system()

IMPLICIT NONE

END 

!-------------------------------------------------------------------------------

LOGICAL FUNCTION ispressure(unit)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: unit

ispressure = .FALSE.

END FUNCTION ispressure

!-------------------------------------------------------------------------------

LOGICAL FUNCTION istime(unit)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: unit

istime = .FALSE.

END FUNCTION istime

#endif

END MODULE fort2c_udunits

!===============================================================================

MODULE input_grid_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE :: netcdf

IMPLICIT NONE

CHARACTER(LEN=nf90_max_name) :: varname
INTEGER :: nlon_in, nlat_in, nlev_in
REAL(KIND=REAL64), ALLOCATABLE :: lon_in(:),lat_in(:),lev_in(:)

INTEGER, ALLOCATABLE :: size_in_x(:),size_in_y(:)
INTEGER, ALLOCATABLE :: start_in_x(:),start_in_y(:)
INTEGER, ALLOCATABLE :: pe_index_in_x(:),pe_index_in_y(:)
REAL(KIND=REAL64), ALLOCATABLE :: lon_in_local(:),lat_in_local(:)

END MODULE input_grid_mod

!===============================================================================

MODULE mpp_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE :: mpi
USE :: util_mod

IMPLICIT NONE

INTEGER(KIND=INT32) :: comm
INTEGER(KIND=INT32) :: nproc,mype
INTEGER :: nprocx,nprocy
INTEGER :: iprocx,iprocy

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE check_input()

IF (verbose > 0) WRITE(outunit,*)'nproc,mype =',nproc,mype

IF (nproc /= nprocx*nprocy) THEN

  WRITE(error_unit,*)'ERROR grid domain decomposition doesn''t match number of mpi tasks'
  WRITE(error_unit,*)nprocx,nprocy,nproc

  CALL abort()

END IF

iprocx = MOD(mype,nprocx)
iprocy = mype/nprocx
IF (verbose > 0) WRITE(outunit,*)'PE ',mype,' iprocx,iprocy = ',iprocx,iprocy

END SUBROUTINE check_input

END MODULE mpp_mod


!===============================================================================

MODULE mpp_grid_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE :: mpp_mod
USE :: input_grid_mod
USE :: util_mod

IMPLICIT NONE

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE decomp_grid_0d(nlon,nlat,lon,lat,lon_local,lat_local, &
                          start_x,size_x,start_y,size_y)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nlon,nlat
INTEGER, INTENT(INOUT), ALLOCATABLE :: size_x(:),size_y(:)
INTEGER, INTENT(INOUT), ALLOCATABLE :: start_x(:),start_y(:)

REAL(KIND=REAL64), INTENT(IN) :: lon(:),lat(:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: lon_local(:),lat_local(:)

INTEGER :: npts_x,npts_y,start
INTEGER :: irest,iproc

ALLOCATE(size_x(0:nprocx-1),size_y(0:nprocy-1))
ALLOCATE(start_x(0:nprocx-1),start_y(0:nprocy-1))

start_x = nlon+1
start_y = nlat+1
size_x = 0
size_y = 0
IF (mype == 0) THEN
  start_x(iprocx) = 1
  start_y(iprocy) = 1
  size_x(iprocx) = nlon
  size_y(iprocy) = nlat
END IF

IF (verbose > 0) THEN
  WRITE(outunit,*)'PE ',mype,' local start,size lon = ',start_x(iprocx),size_x(iprocx)
  WRITE(outunit,*)'PE ',mype,' local start,size lat = ',start_y(iprocy),size_y(iprocy)
END IF

ALLOCATE(lon_local(size_x(iprocx)))
ALLOCATE(lat_local(size_y(iprocy)))

IF (mype == 0) THEN
  lon_local = lon
  lat_local = lat
  IF (verbose > 0) THEN
    WRITE(outunit,*)'PE ',mype,' local lon data = ',lon_local
    WRITE(outunit,*)'PE ',mype,' local lat data = ',lat_local
  END IF
END IF

END SUBROUTINE decomp_grid_0d

!-------------------------------------------------------------------------------

SUBROUTINE decomp_grid_2d(nlon,nlat,lon,lat,lon_local,lat_local, &
                          start_x,size_x,start_y,size_y, &
                          pe_index_x,pe_index_y)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nlon,nlat
INTEGER, INTENT(INOUT), ALLOCATABLE :: size_x(:),size_y(:)
INTEGER, INTENT(INOUT), ALLOCATABLE :: start_x(:),start_y(:)
INTEGER, INTENT(INOUT), ALLOCATABLE :: pe_index_x(:),pe_index_y(:)

REAL(KIND=REAL64), INTENT(IN) :: lon(:),lat(:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: lon_local(:),lat_local(:)

INTEGER :: npts_x,npts_y,start,ipt
INTEGER :: irest,iproc
INTEGER :: prow_s,prow_n

ALLOCATE(size_x(0:nprocx-1),size_y(0:nprocy-1))
ALLOCATE(start_x(0:nprocx-1),start_y(0:nprocy-1))
ALLOCATE(pe_index_x(nlon),pe_index_y(nlat))

start = 1
npts_x = nlon / nprocx
irest = nlon - (nprocx*npts_x)
DO iproc = 0,nprocx-1
  start_x(iproc) = start
  !size_x(iproc) = npts_x
  !IF (irest > 0) THEN
  !  size_x(iproc) = size_x(iproc) + 1
  !  irest = irest - 1
  !END IF
  ! Decompose grid the same way as the UM (irest must be even)
  ! Work out the decomposition in the East-West direction. As far as
  ! possible each processor has the same local row length. However, if
  ! this is not possible, the extra points are distributed symetrically
  ! such that each processor has the same number of points as the
  ! processor on the opposite side of the globe.
  IF (iproc < nprocx/2) THEN
    IF (iproc < irest/2) THEN
      size_x(iproc) = npts_x + 1
    ELSE
      size_x(iproc) = npts_x
    END IF
  ElSE
    IF (iproc - (nprocx/2) < irest/2) THEN
      size_x(iproc) = npts_x + 1
    ELSE
      size_x(iproc) = npts_x
    END IF
  END IF
  start = start + size_x(iproc)
  DO ipt=start_x(iproc),start_x(iproc)+size_x(iproc)-1
    pe_index_x(ipt)=iproc
  END DO
END DO

IF (verbose > 0) WRITE(outunit,*)'PE ',mype,' local start,size lon = ',start_x(iprocx),size_x(iprocx)

ALLOCATE(lon_local(size_x(iprocx)))
lon_local = lon(start_x(iprocx):start_x(iprocx)+size_x(iprocx)-1)

IF (verbose > 0) WRITE(outunit,*)'PE ',mype,' local lon data = ',lon_local

start = 1
npts_y = nlat / nprocy
irest = nlat - (nprocy*npts_y)
!DO iproc = 0,nprocy-1
!  start_y(iproc) = start
!  size_y(iproc) = npts_y
!  IF (irest > 0) THEN
!    size_y(iproc) = size_y(iproc) + 1
!    irest = irest - 1
!  END IF
!  start = start + size_y(iproc)
!  DO ipt=start_y(iproc),start_y(iproc)+size_y(iproc)-1
!    pe_index_y(ipt)=iproc
!  END DO
!END DO

! Decompose grid the same way as the UM
! Work out the decomposition in the North-South direction. As far as
! possible each processor has the same number of rows. However, if this
! is not possible, the extra rows are distributed thus:
! - an extra row is given to the Northern most processor
! - the remaining extra rows are distributed symetrically around the
!   equator, starting at the processor(s) closest to the equator.

DO iproc = 0,nprocy-1
  size_y(iproc) = npts_y
END DO

IF (irest  >  0) THEN
  ! give Southern most processors an extra row
  size_y(0)=npts_y+1
  irest=irest-1
END IF

! Set up pointers to processor rows to which we will add extra rows
! to. These start around the equator, and will work out towards
! the poles.

IF (MOD(nprocy,2)  ==  0) THEN  ! Even number of NS processors
  prow_s=nprocy/2
  prow_n=prow_s+1
ELSE  ! Odd number of NS processors
  prow_s=(nprocy/2)+1
  prow_n=prow_s
END IF

DO WHILE (irest  >  0)

  IF (prow_n  ==  prow_s) THEN
    size_y(prow_n-1)=npts_y+1
    irest=irest-1
  ELSE
    size_y(prow_s-1)=npts_y+1
    irest=irest-1
    IF (irest  >=  1) THEN
      size_y(prow_n-1)=npts_y+1
      irest=irest-1
    END IF
  END IF

  prow_s=MAX(1,prow_s-1)
  prow_n=MIN(nprocy,prow_n+1)

END DO

DO iproc=0,nprocy-1
  start_y(iproc)=start
  start=start+size_y(iproc)
  DO ipt=start_y(iproc),start_y(iproc)+size_y(iproc)-1
    pe_index_y(ipt)=iproc
  END DO
END DO

IF (verbose > 0) WRITE(outunit,*)'PE ',mype,' local start,size lat = ',start_y(iprocy),size_y(iprocy)

ALLOCATE(lat_local(size_y(iprocy)))
lat_local = lat(start_y(iprocy):start_y(iprocy)+size_y(iprocy)-1)

IF (verbose > 0) WRITE(outunit,*)'PE ',mype,' local lat data = ',lat_local
IF (verbose > 0) FLUSH(outunit)

END SUBROUTINE decomp_grid_2d

!-------------------------------------------------------------------------------

SUBROUTINE decomp_grid_reduce1(nlat,ntot, &
                               lon,lat, &
                               bnds_lon,bnds_lat, &
                               lon_local,lat_local, &
                               bnds_lon_local,bnds_lat_local, &
                               start_x,size_x)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nlat,ntot
INTEGER, INTENT(INOUT), ALLOCATABLE :: start_x(:)
INTEGER, INTENT(INOUT), ALLOCATABLE :: size_x(:)

REAL(KIND=REAL64), INTENT(IN) :: lon(:),lat(:)
REAL(KIND=REAL64), INTENT(IN) :: bnds_lon(:,:),bnds_lat(:,:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: lon_local(:),lat_local(:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: bnds_lon_local(:,:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: bnds_lat_local(:,:)

INTEGER :: count1,count2,i
INTEGER(KIND=INT32) :: ierror

REAL(KIND=REAL64) :: lon_step,lat_step,lon1,lat1
REAL(KIND=REAL64) :: lonmin,lonmax,latmin,latmax

ALLOCATE(size_x(0:nproc-1))
ALLOCATE(start_x(0:nproc-1))

lon_step = lon_in(2) - lon_in(1)
lat_step = lat_in(2) - lat_in(1)
IF (verbose > 0) WRITE(outunit,*) 'lon_step,lat_step = ',lon_step,lat_step

lonmin = lon_in_local(1) - 0.5*lon_step
lonmax = lon_in_local(SIZE(lon_in_local)) + 0.5*lon_step
latmin = lat_in_local(1) - 0.5*lat_step
latmax = lat_in_local(SIZE(lat_in_local)) + 0.5*lat_step
IF (verbose > 0) WRITE(outunit,*) 'lonmin,lonmax,latmin,latmax = ',lonmin,lonmax,latmin,latmax

count1 = 0
DO i=1,ntot
  lon1 = lon(i)
  lat1 = lat(i)
  IF (lon1 >= lonmax) THEN
    lon1 = lon1 - 360.0
  ELSE IF (lon1 < lonmin) THEN
    lon1 = lon1 + 360.0
  END IF
  IF ((lon1 >= lonmin .AND. lon1 < lonmax) .AND. &
      (lat1 >= latmin .AND. lat1 < latmax)) THEN
    count1 = count1 + 1
    IF (verbose > 0) THEN
      WRITE(outunit,*) '*** on local grid = ',lon1,lonmin,lonmax
      WRITE(outunit,*) '*** on local grid = ',lat1,latmin,latmax
    ENDIF
  ELSE
    IF (verbose > 0) THEN
      WRITE(outunit,*) 'not on local grid = ',lon1,lonmin,lonmax
      WRITE(outunit,*) 'not on local grid = ',lat1,latmin,latmax
    ENDIF
  END IF
END DO
IF (verbose > 0) WRITE(outunit,*)'PE ',mype,' nlat,ntot,count1 = ',nlat,ntot,count1

CALL mpi_allgather(count1,1,MPI_INTEGER,size_x,1,MPI_INTEGER,comm,ierror)
IF (verbose > 1) THEN
  DO i=0,nproc-1
    WRITE(outunit,'(A,I0,2X,I0)') 'iproc,size_x = ',i,size_x(i)
  END DO
END IF
IF (verbose > 0) WRITE(outunit,'(A,I0,2X,I0)') 'SUM(size_x) = ',SUM(size_x),ntot
IF (SUM(size_x) /= ntot) THEN
  WRITE(error_unit,*) 'Sum of local gridpoint values invalid ',SUM(size_x),ntot
  CALL abort()
END IF

start_x(0) = 1
DO i=1,nproc-1
   start_x(i) = start_x(i-1) + size_x(i-1)
END DO
IF (verbose > 1) THEN
  DO i=0,nproc-1
    WRITE(outunit,'(A,I0,2X,I0)') 'iproc,start_x = ',i,start_x(i)
  END DO
END IF

ALLOCATE(lon_local(count1))
ALLOCATE(lat_local(count1))
ALLOCATE(bnds_lon_local(4,count1))
ALLOCATE(bnds_lat_local(4,count1))

count2 = 0
DO i=1,ntot
  lon1 = lon(i)
  lat1 = lat(i)
  IF (lon1 >= lonmax) THEN
    lon1 = lon1 - 360.0
  ELSE IF (lon1 < lonmin) THEN
    lon1 = lon1 + 360.0
  END IF
  IF ((lon1 >= lonmin .AND. lon1 < lonmax) .AND. &
      (lat1 >= latmin .AND. lat1 < latmax)) THEN
    count2 = count2 + 1
    IF (count2 > count1) THEN
      WRITE(error_unit,*) 'Number of local gridpoints exceeded ',count1,count2
      CALL abort()
    END IF
    lon_local(count2) = lon(i)
    lat_local(count2) = lat(i)
    bnds_lon_local(:,count2) = bnds_lon(:,i)
    bnds_lat_local(:,count2) = bnds_lat(:,i)
    IF (verbose > 2) WRITE(outunit,*)'lat,lon = ',lon(i),lat(i)
  END IF
END DO
IF (count2 /= count1) THEN
  WRITE(error_unit,*) 'Number of local gridpoints not consistent ',count1,count2
  CALL abort()
END IF

IF (verbose > 0) FLUSH(outunit)

END SUBROUTINE decomp_grid_reduce1

!-------------------------------------------------------------------------------

SUBROUTINE decomp_grid_reduce2(nlat,ntot, &
                               lon,lat, &
                               bnds_lon,bnds_lat, &
                               lon_local,lat_local, &
                               bnds_lon_local,bnds_lat_local, &
                               start_x,size_x)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nlat,ntot
INTEGER, INTENT(INOUT), ALLOCATABLE :: start_x(:)
INTEGER, INTENT(INOUT), ALLOCATABLE :: size_x(:)

REAL(KIND=REAL64), INTENT(IN) :: lon(:),lat(:)
REAL(KIND=REAL64), INTENT(IN) :: bnds_lon(:,:),bnds_lat(:,:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: lon_local(:),lat_local(:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: bnds_lon_local(:,:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: bnds_lat_local(:,:)

INTEGER :: count1,count2,i,j,ii
INTEGER :: pe_col,pe_row,proc
INTEGER :: which_proc(ntot)
INTEGER(KIND=INT32) :: ierror

REAL(KIND=REAL64) :: base_xi1,delta_xi1
REAL(KIND=REAL64) :: base_xi2,delta_xi2

ALLOCATE(size_x(0:nproc-1))
ALLOCATE(start_x(0:nproc-1))

count1 = 0
base_xi1 = 0.0
delta_xi1 = 360.0/nlon_in
base_xi2 = -90.0
delta_xi2 = 180.0/nlat_in

DO ii=1,ntot
  i = NINT(((lon(ii) - base_xi1)/delta_xi1) + 0.5)
  IF (i < 1) THEN
    i = 1
  ELSE IF (i > nlon_in) THEN
    i = nlon_in
  END IF

  j = NINT(((lat(ii) - base_xi2)/delta_xi2) + 0.5)
  IF (j < 1) THEN
    j = 1
  ELSE IF (j > nlat_in) THEN
    j = nlat_in
  END IF

  pe_col = pe_index_in_x(i)
  pe_row = pe_index_in_y(j)
  proc = pe_row*nprocx + pe_col
  which_proc(ii) = proc
  IF (proc == mype) THEN
    count1 = count1 + 1
  END IF
END DO
IF (verbose > 0) WRITE(outunit,*)'PE ',mype,' nlat,ntot,count1 = ',nlat,ntot,count1

CALL mpi_allgather(count1,1,MPI_INTEGER,size_x,1,MPI_INTEGER,comm,ierror)
IF (verbose > 1) THEN
  DO i=0,nproc-1
    WRITE(outunit,'(A,I0,2X,I0)') 'iproc,size_x = ',i,size_x(i)
  END DO
END IF
IF (verbose > 0) WRITE(outunit,'(A,I0,2X,I0)') 'SUM(size_x) = ',SUM(size_x),ntot

start_x(0) = 1
DO i=1,nproc-1
   start_x(i) = start_x(i-1) + size_x(i-1)
END DO
IF (verbose > 1) THEN
  DO i=0,nproc-1
    WRITE(outunit,'(A,I0,2X,I0)') 'iproc,start_x = ',i,start_x(i)
  END DO
END IF

ALLOCATE(lon_local(count1))
ALLOCATE(lat_local(count1))
ALLOCATE(bnds_lon_local(4,count1))
ALLOCATE(bnds_lat_local(4,count1))

count2 = 0
DO ii=1,ntot
  IF (which_proc(ii) == mype) THEN
    count2 = count2 + 1
    lon_local(count2) = lon(ii)
    lat_local(count2) = lat(ii)
    bnds_lon_local(:,count2) = bnds_lon(:,ii)
    bnds_lat_local(:,count2) = bnds_lat(:,ii)
    IF (verbose > 2) WRITE(outunit,*)'lat,lon = ',lon(ii),lat(ii)
  END IF
END DO

IF (verbose > 0) FLUSH(outunit)

END SUBROUTINE decomp_grid_reduce2

!-------------------------------------------------------------------------------

SUBROUTINE scatter_data(data,local_data,start_x,size_x,start_y,size_y)

IMPLICIT NONE

REAL(KIND=REAL64), INTENT(IN) :: data(:,:,:)
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: local_data(:,:,:)

INTEGER, INTENT(IN) :: size_x(0:),size_y(0:)
INTEGER, INTENT(IN) :: start_x(0:),start_y(0:)

INTEGER :: nx,ny,nz
INTEGER :: iz

nx = size_x(iprocx)
ny = size_y(iprocy)
nz = nlev_in

ALLOCATE(local_data(nx,ny,nz))

DO iz=1,nz
  CALL scatter_data_2D(data(:,:,iz),local_data(:,:,iz), &
                       start_in_x,size_in_x,start_in_y,size_in_y)
END DO

END SUBROUTINE scatter_data

!-------------------------------------------------------------------------------

SUBROUTINE scatter_data_2D(data,local_data,start_x,size_x,start_y,size_y)

IMPLICIT NONE

REAL(KIND=REAL64), INTENT(IN) :: data(:,:)
REAL(KIND=REAL64), INTENT(INOUT) :: local_data(:,:)

INTEGER, INTENT(IN) :: size_x(0:),size_y(0:)
INTEGER, INTENT(IN) :: start_x(0:),start_y(0:)

REAL(KIND=REAL64),ALLOCATABLE :: data_reorder(:)

INTEGER :: nx,ny,ii,jx,jy,iproc,jproc,proc
INTEGER(KIND=INT32) :: local_data_size,sendcounts(nproc),displs(nproc),ierror

nx = size_x(iprocx)
ny = size_y(iprocy)
local_data_size = nx*ny

! Copy each PE chunk of 2D data into contiguous 1D array suitable for mpi_scatterv
IF (mype == 0) THEN
  ALLOCATE(data_reorder(SIZE(data)))
  ii=1
  displs(1) = 0
  DO proc=0,nproc-1
    iproc = MOD(proc,nprocx)
    jproc = proc/nprocx
    sendcounts(proc+1) = size_x(iproc)*size_y(jproc)
    IF (proc < nproc-1) THEN
      displs(proc+2) = displs(proc+1) + sendcounts(proc+1)
    END IF

    DO jy=start_y(jproc),start_y(jproc)+size_y(jproc)-1
      DO  jx=start_x(iproc),start_x(iproc)+size_x(iproc)-1
        data_reorder(ii) = data(jx,jy)
        ii = ii + 1
      END DO
    END DO

  END DO
END IF

CALL mpi_scatterv(data_reorder,sendcounts,displs,MPI_DOUBLE_PRECISION, &
                  local_data,local_data_size,MPI_DOUBLE_PRECISION,0,comm,ierror)

IF (mype == 0) THEN
  DEALLOCATE(data_reorder)
END IF

END SUBROUTINE scatter_data_2D

END MODULE mpp_grid_mod

!===============================================================================

MODULE timing_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE :: mpp_mod

IMPLICIT NONE

CHARACTER(LEN=256) :: message

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE get_timing(message, rank, first)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: message
INTEGER(KIND=INT32), OPTIONAL, INTENT(IN) :: rank
LOGICAL, OPTIONAL, INTENT(IN) :: first

INTEGER :: datetime(8)
LOGICAL :: local_first
REAL, SAVE ::time0
REAL :: time1
CHARACTER(LEN=23) :: cdatetime

IF (PRESENT(rank)) THEN
  IF (rank /= mype) RETURN
END IF

IF (PRESENT(first)) THEN
  local_first = first
ELSE
  local_first = .FALSE.
END IF

CALL DATE_AND_TIME(values=datetime)
WRITE(cdatetime,'(i4.4,A1,i2.2,A1,i2.2,A1,i2.2,A1,i2.2,A1,i2.2,A1,i3.3)') &
    datetime(1),'/',datetime(2),'/',datetime(3),' ', &
    datetime(5),':',datetime(6),':',datetime(7),'.',datetime(8)

IF (local_first) THEN
   CALL CPU_TIME(time0)
   WRITE(outunit,*)cdatetime,' : ',TRIM(message)
ELSE
   CALL CPU_TIME(time1)
   WRITE(outunit,*)cdatetime,' : ',TRIM(message),' elasped cpu time = ',time1-time0
END IF

FLUSH(outunit)

END SUBROUTINE get_timing

END MODULE timing_mod

!===============================================================================

MODULE read_netcdf_file_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_NULL_CHAR
USE :: netcdf
USE :: mpi
USE :: mpp_mod
USE :: input_grid_mod
USE :: fort2c_udunits

IMPLICIT NONE

INTEGER, PARAMETER :: max_attname = LEN('standard_name')

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE read_netcdf_file(filename,varname,data)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: filename,varname
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: data(:,:,:)

INTEGER(KIND=INT32) :: ncid,status

IF (mype == 0) THEN
  WRITE(outunit,*) 'Opening netCDF file ',TRIM(filename)
  status = nf90_open(filename, nf90_nowrite, ncid)
  IF (status /= nf90_noerr) CALL netcdf_error(status,filename)
END IF

CALL get_dim_data(ncid,varname,'X',nlon_in,lon_in)
CALL get_dim_data(ncid,varname,'Y',nlat_in,lat_in)
CALL get_dim_data(ncid,varname,'Z',nlev_in,lev_in)

CALL get_data(ncid,varname,nlon_in,nlat_in,nlev_in,data)

IF (mype == 0) THEN
  WRITE(outunit,*) 'Closing netCDF file ',TRIM(filename)
  status = nf90_close(ncid)
  IF (status /= nf90_noerr) CALL netcdf_error(status,filename)
END IF

END SUBROUTINE read_netcdf_file

!-------------------------------------------------------------------------------

SUBROUTINE get_dim_data(ncid,varname,axis,len,data)

IMPLICIT NONE

INTEGER(KIND=INT32), INTENT(IN) :: ncid
CHARACTER(LEN=*), INTENT(IN) :: varname
CHARACTER(LEN=1), INTENT(IN) :: axis
INTEGER(KIND=INT32), INTENT(OUT) :: len
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: data(:)

INTEGER(KIND=INT32) :: status,dimid,varid,ndims,dimvarid
INTEGER(KIND=INT32) :: start(1),count(1)

LOGICAL :: found_axis
CHARACTER(LEN=1) :: dimaxis
CHARACTER(LEN=nf90_max_name) :: dimname

IF (mype == 0) THEN
  status = nf90_inq_varid(ncid, varname, varid)
  IF (status /= nf90_noerr) CALL netcdf_error(status,varname)

  status = nf90_inquire_variable(ncid, varid, ndims=ndims)
  IF (status /= nf90_noerr) CALL netcdf_error(status,varname)

  found_axis = .FALSE.
  DO dimid=1,ndims
    status = nf90_inquire_dimension(ncid,dimid,len=len,name=dimname)
    IF (status /= nf90_noerr) CALL netcdf_error(status,varname)

    CALL get_dim_axis(ncid,dimname,dimaxis)
    IF (dimaxis == axis) THEN
      WRITE(outunit,*)'Found dimension ',TRIM(dimname),' for ',axis,' axis size = ',len
      found_axis = .TRUE.
      EXIT
    END IF
  END DO
  IF (.NOT. found_axis) THEN
    WRITE(error_unit,*)'No dimension variable found for axis ',axis
    len = 1
  END IF
END IF

CALL mpi_bcast(len,1,MPI_INTEGER,0,comm,status)
ALLOCATE(data(len))

IF (mype == 0) THEN
  IF (found_axis) THEN
    status = nf90_inq_varid(ncid,dimname,dimvarid)
    IF (status /= nf90_noerr) CALL netcdf_error(status,dimname)
    
    start(1) = 1
    count(1) = len
    status = nf90_get_var(ncid,dimvarid,data,start,count)
  ELSE
    data(1) = 0.0
  END IF
END IF

CALL mpi_bcast(data,len,MPI_DOUBLE_PRECISION,0,comm,status)

END SUBROUTINE get_dim_data

!-------------------------------------------------------------------------------

SUBROUTINE get_data(ncid,name,nx,ny,nz,data)

IMPLICIT NONE

INTEGER(KIND=INT32), INTENT(IN) :: ncid
CHARACTER(LEN=*), INTENT(IN) :: name
INTEGER(KIND=INT32), INTENT(IN) :: nx,ny,nz
REAL(KIND=REAL64), INTENT(INOUT), ALLOCATABLE :: data(:,:,:)

INTEGER(KIND=INT32) :: status,varid
INTEGER(KIND=INT32) :: start(3),count(3)
INTEGER :: iz

ALLOCATE(data(nx,ny,nz))

IF (mype == 0) THEN
  status = nf90_inq_varid(ncid,name,varid)
  IF (status /= nf90_noerr) CALL netcdf_error(status,name)

  DO iz=1,nz
    start(1) = 1
    start(2) = 1
    start(3) = iz
    count(1) = nx
    count(2) = ny
    count(3) = iz
    status = nf90_get_var(ncid,varid,data(:,:,iz),start,count)
  END DO
END IF

END SUBROUTINE get_data

!-------------------------------------------------------------------------------

SUBROUTINE get_dim_axis(ncid,dimname,axis)

IMPLICIT NONE

INTEGER(KIND=INT32), INTENT(IN) :: ncid
CHARACTER(LEN=*), INTENT(IN) :: dimname
CHARACTER(LEN=1), INTENT(OUT) :: axis

INTEGER(KIND=INT32) :: status,dimid,varid
INTEGER :: i
CHARACTER(LEN=1024) :: value
CHARACTER(LEN=max_attname), PARAMETER :: attnames(*) = &
      [CHARACTER(LEN=max_attname) :: 'axis','standard_name','units','positive']

axis = ' '

status = nf90_inq_dimid(ncid,dimname,dimid)
IF (status /= nf90_noerr) CALL netcdf_error(status,dimname)

status = nf90_inq_varid(ncid,dimname,varid)
IF (status /= nf90_noerr) CALL netcdf_error(status,dimname)

DO i = 1,SIZE(attnames)
  status = nf90_get_att(ncid,varid,TRIM(attnames(i)),value)
  IF (status == nf90_noerr) THEN
    ! Remove trailing \0 from files written by c code
    IF (value(LEN_TRIM(value):LEN_TRIM(value)) == C_NULL_CHAR) &
      value(LEN_TRIM(value):LEN_TRIM(value)) = ' '
    SELECT CASE (TRIM(attnames(i)))
    CASE ('axis')
      CALL upcase(value)
      SELECT CASE (TRIM(value))
      CASE ('X','Y','Z','T')
        axis = value(1:1)
        EXIT
      END SELECT
    CASE ('standard_name')
      SELECT CASE (TRIM(value))
      CASE ('longitude')
        axis = 'X'
        EXIT
      CASE ('latitude')
        axis = 'Y'
        EXIT
      CASE ('air_pressure','altitude','height','depth', &
            'atmosphere_hybrid_height_coordinate', &
            'atmosphere_hybrid_sigma_pressure_coordinate', &
            'atmosphere_sigma_coordinate', &
            'atmosphere_ln_pressure_coordinate', &
            'atmosphere_sleve_coordinate', &
            'model_level_number', &
            'air_potential_temperature','ertel_potential_vorticity')
        axis = 'Z'
        EXIT
      CASE ('time')
        axis = 'T'
        EXIT
      END SELECT
    CASE ('units')
      SELECT CASE (TRIM(value))
      CASE ('degrees_east','degree_east','degrees_E', &
            'degree_E','degreesE','degreeE')
        axis = 'X'
        EXIT
      CASE ('degrees_north','degree_north','degrees_N', &
            'degree_N','degreesN','degreeN')
        axis = 'Y'
        EXIT
      END SELECT
      IF (read_xml() /= 0) EXIT
      IF (ispressure(value)) THEN
        axis = 'Z'
        CALL free_system()
        EXIT
      END IF
      IF (istime(value)) THEN
        axis = 'T'
        CALL free_system()
        EXIT
      END IF
      CALL free_system()
      CALL locase(value)
      IF (INDEX(value,' since ') > 0) THEN
        axis = 'T'
        EXIT
      END IF
    CASE ('positive')
      CALL locase(value)
      SELECT CASE (TRIM(value))
      CASE ('up','down')
        axis = 'Z'
        EXIT
      END SELECT
    END SELECT
  ELSE IF (status /= nf90_enotatt) THEN
    CALL netcdf_error(status,TRIM(attnames(i)))
  END IF
END DO

END SUBROUTINE get_dim_axis

!-------------------------------------------------------------------------------

SUBROUTINE netcdf_error(status,message)

IMPLICIT NONE

INTEGER(KIND=INT32), INTENT(IN) :: status
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

IF (PRESENT(message)) THEN
  WRITE(error_unit,*) TRIM(message)//' : ',TRIM(nf90_strerror(status))
ELSE
  WRITE(error_unit,*) TRIM(nf90_strerror(status))
END IF

CALL abort()

END SUBROUTINE netcdf_error

END MODULE read_netcdf_file_mod

!===============================================================================

MODULE xios_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE :: xios
USE :: mpp_mod
USE :: timing_mod
USE :: input_grid_mod

IMPLICIT NONE

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE xios_init(do_ens,ensemble_member)

IMPLICIT NONE

LOGICAL, INTENT(IN) :: do_ens
INTEGER(KIND=INT32), INTENT(IN) :: ensemble_member

INTEGER(KIND=INT32) :: key,ens_comm,ierror

CALL xios_initialize('client',return_comm=comm)

CALL xios_context_initialize('xios_test',comm)

IF (do_ens) THEN
  key = 0
  CALL mpi_comm_split(comm,ensemble_member,key,ens_comm,ierror)
  comm = ens_comm
END IF

CALL mpi_comm_size(comm,nproc,ierror)
CALL mpi_comm_rank(comm,mype,ierror)

IF (ierror /= MPI_SUCCESS) THEN
  CALL xios_finalize()
  WRITE(error_unit,*)'Error initialising MPI, ierror = ',ierror
  CALL abort()
END IF

END SUBROUTINE xios_init

!-------------------------------------------------------------------------------

SUBROUTINE xios_def()

IMPLICIT NONE

INTEGER(KIND=INT32) :: ierror

TYPE(xios_duration) :: timestep

CALL get_timing('Start of xios_def')

CALL xios_set_domain_attr('input_domain',type='rectilinear', &
                          ni_glo=nlon_in,nj_glo=nlat_in, &
                          ibegin=start_in_x(iprocx)-1,ni=size_in_x(iprocx), &
                          jbegin=start_in_y(iprocy)-1,nj=size_in_y(iprocy), &
                          data_dim=2, &
                          lonvalue_1d=lon_in_local, &
                          latvalue_1d=lat_in_local)
CALL get_timing('After xios_set_domain_attr input grid')

! 1 hour timestep
timestep%second = 3600
CALL xios_define_calendar('Gregorian',timestep)
CALL get_timing('After xios_define_calendar')

CALL get_timing('End of xios_def')

END SUBROUTINE xios_def

!-------------------------------------------------------------------------------

SUBROUTINE xios_final()

IMPLICIT NONE

CALL xios_context_finalize()

CALL xios_finalize()

END SUBROUTINE xios_final

!-------------------------------------------------------------------------------

SUBROUTINE xios_clone_field(ncopies,do_regrid,do_ens,field_id)

IMPLICIT NONE

INTEGER, INTENT(IN) :: ncopies
LOGICAL, INTENT(IN) :: do_regrid,do_ens
CHARACTER(LEN=*), INTENT(IN) :: field_id

INTEGER, PARAMETER :: max_size = 256
INTEGER :: i
LOGICAL :: do_regrid_local,do_ens_local,do_ens_mean_local
LOGICAL :: grid_ref_exists,check_if_active_exists
LOGICAL :: check_if_active
CHARACTER(LEN=max_size) :: grid_in,grid_regrid,grid_ens_mean
CHARACTER(LEN=max_size) :: field_regrid_id,field_copy_id,field_regrid_copy_id
CHARACTER(LEN=max_size) :: long_name

TYPE(xios_fieldgroup) :: field_def_handle
TYPE(xios_field) :: field_handle
TYPE(xios_file) :: file_handle,file_handle_ens

IF (ncopies == 0) RETURN

do_regrid_local = do_regrid
do_ens_local = do_ens
do_ens_mean_local = do_ens

IF (xios_is_valid_file('output')) THEN
  CALL xios_get_handle('output',file_handle)
ELSE
  WRITE(outunit,*)'Unknown XIOS file element "output"'
  do_ens_local = .FALSE.
  do_regrid_local = .FALSE.
END IF

IF (.NOT. (do_regrid_local .OR. do_ens_local .OR. do_ens_mean_local)) RETURN

CALL xios_get_handle('field_definition',field_def_handle)

CALL xios_is_defined_field_attr(field_id,grid_ref=grid_ref_exists, &
                                         check_if_active=check_if_active_exists)
IF (grid_ref_exists) THEN
  CALL xios_get_field_attr(field_id,grid_ref=grid_in)
  WRITE(outunit,*)'Grid reference for field ',TRIM(field_id),' is ',TRIM(grid_in)
ELSE
  WRITE(outunit,*)'Field ',TRIM(field_id),' has no grid reference'
  RETURN
END IF
IF (check_if_active_exists) THEN
  CALL xios_get_field_attr(field_id,check_if_active=check_if_active)
END IF

IF (do_regrid_local) THEN
  field_regrid_id = TRIM(field_id)//'_regrid'
  IF (xios_is_valid_field(field_regrid_id)) THEN
    CALL xios_get_field_attr(field_regrid_id,grid_ref=grid_regrid)
    WRITE(outunit,*)'Grid reference for field ',TRIM(field_regrid_id),' is ',TRIM(grid_regrid)
  ELSE
    do_regrid_local = .FALSE.
    WRITE(outunit,*)'Regrid field ',TRIM(field_regrid_id),' unknown. No regrid field copies'
  END IF
END IF

IF (do_ens_mean_local) THEN
  IF (xios_is_valid_file('output_ens')) THEN
    CALL xios_get_handle('output_ens',file_handle_ens)
  ELSE
    do_ens_mean_local = .FALSE.
    WRITE(outunit,*)'Unknown XIOS file element "output_ens"'
  END IF

  grid_ens_mean = 'grid_ens_mean'
  IF (xios_is_valid_grid(grid_ens_mean)) THEN
    WRITE(outunit,*)'Grid for ensemble mean is ',TRIM(grid_ens_mean)
  ELSE
    do_ens_mean_local = .FALSE.
    WRITE(outunit,*)'Unknown grid for ensemble mean ',TRIM(grid_ens_mean)
  END IF
END IF

IF (.NOT. (do_regrid_local .OR. do_ens_local .OR. do_ens_mean_local)) RETURN

DO i=1,ncopies
  WRITE(field_copy_id,'(A,I0)') TRIM(field_id)//'_copy_',i
  WRITE(long_name,'(A,I0,A)') 'Copy ',i,' of '//TRIM(field_id)

  CALL xios_add_child(field_def_handle,field_handle,field_copy_id)
  CALL xios_set_attr(field_handle,grid_ref=grid_in,long_name=long_name)
  IF (check_if_active_exists) THEN
    CALL xios_set_attr(field_handle,check_if_active=check_if_active)
  END IF

  IF (do_regrid_local .OR. do_ens_local) THEN
    CALL xios_add_child(file_handle,field_handle)
    CALL xios_set_attr(field_handle,field_ref=field_copy_id,operation='instant')
  END IF
  
  IF (do_regrid_local) THEN
    WRITE(field_regrid_copy_id,'(A,I0)') TRIM(field_regrid_id)//'_copy_',i
    WRITE(long_name,'(A,I0,A)') 'Copy ',i,' of '//TRIM(field_regrid_id)

    CALL xios_add_child(field_def_handle,field_handle,field_regrid_copy_id)
    CALL xios_set_attr(field_handle,field_ref=field_copy_id, &
                       grid_ref=grid_regrid,long_name=long_name)
    
    CALL xios_add_child(file_handle,field_handle)
    CALL xios_set_attr(field_handle,field_ref=field_regrid_copy_id,operation='instant')
  END IF

  IF (do_ens_mean_local) THEN
    CALL xios_add_child(file_handle_ens,field_handle)
    CALL xios_set_attr(field_handle,field_ref=field_copy_id,operation='instant', &
                       name=TRIM(field_copy_id)//'_ens',grid_ref=grid_ens_mean)
  END IF
END DO

END SUBROUTINE xios_clone_field

END MODULE xios_mod
