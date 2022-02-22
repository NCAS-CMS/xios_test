
!===============================================================================

MODULE regrid_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE, INTRINSIC :: ISO_C_BINDING
USE :: xios
USE :: mpp_mod
USE :: mpp_grid_mod
USE :: util_mod

IMPLICIT NONE

INTERFACE
  INTEGER FUNCTION c_grib_get_gaussian_latitudes(trunc, lats)                  &
    BIND(c,NAME="grib_get_gaussian_latitudes")

  IMPORT ::  C_LONG, C_DOUBLE

  INTEGER(KIND=C_LONG), VALUE :: trunc
  REAL(KIND=C_DOUBLE) :: lats(*)

  END FUNCTION
END INTERFACE

INTEGER :: nlon_regrid,nlat_regrid,ntot_regrid,nvertex
LOGICAL :: lreduce
REAL(KIND=REAL64), ALLOCATABLE :: lon_regrid(:),lat_regrid(:)
REAL(KIND=REAL64), ALLOCATABLE :: bnds_lon_regrid(:,:),bnds_lat_regrid(:,:)

INTEGER, ALLOCATABLE :: size_regrid_x(:),size_regrid_y(:)
INTEGER, ALLOCATABLE :: start_regrid_x(:),start_regrid_y(:)
INTEGER, ALLOCATABLE :: pe_index_regrid_x(:),pe_index_regrid_y(:)
REAL(KIND=REAL64), ALLOCATABLE :: lon_regrid_local(:),lat_regrid_local(:)
REAL(KIND=REAL64), ALLOCATABLE :: bnds_lon_regrid_local(:,:)
REAL(KIND=REAL64), ALLOCATABLE :: bnds_lat_regrid_local(:,:)

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE create_regrid()

IMPLICIT NONE

INTEGER :: ierr
INTEGER :: grid_type

NAMELIST /input_regrid/ grid_type,nlat_regrid,nlon_regrid

grid_type = 1

OPEN(UNIT=nmlunit,FILE=infile)
READ(UNIT=nmlunit,NML=input_regrid,IOSTAT=ierr)
IF (ierr > 0) THEN
  WRITE(error_unit,*)'Error reading input_regrid namelist file, IOSTAT = ',ierr
  CALL abort()
END IF
CLOSE(UNIT=nmlunit)
WRITE(UNIT=outunit,NML=input_regrid)

SELECT CASE (grid_type)
CASE (1)
  WRITE(UNIT=outunit,FMT=*)'Regridding to Octahedral reduced Gaussian grid'
  CALL create_regrid_orgg()
CASE (2)
  WRITE(UNIT=outunit,FMT=*)'Regridding to Octahedral reduced regular grid'
  CALL create_regrid_orrg()
CASE (3)
  WRITE(UNIT=outunit,FMT=*)'Regridding to reduced Gaussian grid'
  CALL create_regrid_gg()
CASE (4)
  WRITE(UNIT=outunit,FMT=*)'Regridding to regular T grid'
  CALL create_regrid_t()
CASE (5)
  WRITE(UNIT=outunit,FMT=*)'Regridding to regular UV grid'
  CALL create_regrid_uv()
CASE DEFAULT
  WRITE(UNIT=error_unit,FMT=*)'Grid type number ',grid_type,' unknown'
  CALL abort()
END SELECT

IF (lreduce) THEN
  CALL decomp_grid_reduce2(nlat_regrid,ntot_regrid, &
                           lon_regrid,lat_regrid, &
                           bnds_lon_regrid,bnds_lat_regrid, &
                           lon_regrid_local,lat_regrid_local, &
                           bnds_lon_regrid_local,bnds_lat_regrid_local, &
                           start_regrid_x,size_regrid_x)
ELSE
  CALL decomp_grid_2d(nlon_regrid,nlat_regrid,lon_regrid,lat_regrid, &
                      lon_regrid_local,lat_regrid_local, &
                      start_regrid_x,size_regrid_x,start_regrid_y,size_regrid_y, &
                      pe_index_regrid_x,pe_index_regrid_y)
END IF

IF (mype == 0) THEN
  IF (lreduce) THEN
    WRITE(outunit,*) 'Regrid total/lat size = ',ntot_regrid,nlat_regrid
  ELSE
    WRITE(outunit,*) 'Regrid lon/lat size = ',nlon_regrid,nlat_regrid
  END IF
  IF (verbose > 1) THEN
    WRITE(outunit,*) 'Regrid longitude = ',lon_regrid
    WRITE(outunit,*) 'Regrid latitude = ',lat_regrid
  END IF
END IF

END SUBROUTINE create_regrid

!-------------------------------------------------------------------------------

SUBROUTINE xios_def_regrid()

IMPLICIT NONE

IF (lreduce) THEN
  CALL xios_set_domain_attr('regrid_domain',type='gaussian', &
  !CALL xios_set_domain_attr('regrid_domain',type='unstructured', &
                            ni_glo=ntot_regrid, &
                            ibegin=start_regrid_x(mype)-1, &
                            ni=size_regrid_x(mype), &
                            data_dim=1, &
                            lonvalue_1d=lon_regrid_local, &
                            latvalue_1d=lat_regrid_local, &
                            nvertex=nvertex, &
                            bounds_lon_1d=bnds_lon_regrid_local, &
                            bounds_lat_1d=bnds_lat_regrid_local)
ELSE
  CALL xios_set_domain_attr('regrid_domain',type='rectilinear', &
                            ni_glo=nlon_regrid,nj_glo=nlat_regrid, &
                            ibegin=start_regrid_x(iprocx)-1, &
                            ni=size_regrid_x(iprocx), &
                            jbegin=start_regrid_y(iprocy)-1, &
                            nj=size_regrid_y(iprocy), &
                            data_dim=2, &
                            lonvalue_1d=lon_regrid_local, &
                            latvalue_1d=lat_regrid_local)

END IF

END SUBROUTINE xios_def_regrid

!-------------------------------------------------------------------------------

SUBROUTINE cleanup_regrid()

IF (ALLOCATED(lon_regrid)) DEALLOCATE(lon_regrid)
IF (ALLOCATED(lat_regrid)) DEALLOCATE(lat_regrid)
IF (ALLOCATED(size_regrid_x)) DEALLOCATE(size_regrid_x)
IF (ALLOCATED(size_regrid_y)) DEALLOCATE(size_regrid_y)
IF (ALLOCATED(start_regrid_x)) DEALLOCATE(start_regrid_x)
IF (ALLOCATED(start_regrid_y)) DEALLOCATE(start_regrid_y)
IF (ALLOCATED(pe_index_regrid_x)) DEALLOCATE(pe_index_in_x)
IF (ALLOCATED(pe_index_regrid_y)) DEALLOCATE(pe_index_in_y)
IF (ALLOCATED(lon_regrid_local)) DEALLOCATE(lon_regrid_local)
IF (ALLOCATED(lat_regrid_local)) DEALLOCATE(lat_regrid_local)
IF (ALLOCATED(bnds_lon_regrid)) DEALLOCATE(bnds_lon_regrid)
IF (ALLOCATED(bnds_lat_regrid)) DEALLOCATE(bnds_lat_regrid)

END SUBROUTINE cleanup_regrid

!-------------------------------------------------------------------------------

SUBROUTINE create_regrid_uv()

IMPLICIT NONE

REAL(KIND=REAL64) :: lon0,lat0,lon_step,lat_step
INTEGER :: i

lreduce = .FALSE.

ALLOCATE(lon_regrid(nlon_regrid))
ALLOCATE(lat_regrid(nlat_regrid))

lon_step = 360.0_REAL64 / nlon_regrid
lat_step = 180.0_REAL64 / (nlat_regrid-1)
lon0 = 0.0_REAL64
lat0 = -90.0_REAL64

DO i=1,nlon_regrid
  lon_regrid(i) = lon0 + (i-1)*lon_step
END DO

DO i=1,nlat_regrid
  lat_regrid(i) = lat0 + (i-1)*lat_step
END DO

END SUBROUTINE create_regrid_uv

!-------------------------------------------------------------------------------

SUBROUTINE create_regrid_t()

IMPLICIT NONE

REAL(KIND=REAL64) :: lon0,lat0,lon_step,lat_step
INTEGER :: i

lreduce = .FALSE.

ALLOCATE(lon_regrid(nlon_regrid))
ALLOCATE(lat_regrid(nlat_regrid))

lon_step = 360.0_REAL64 / nlon_regrid
lat_step = 180.0_REAL64 / nlat_regrid
lon0 = 0.5_REAL64 * lon_step
lat0 = -90.0_REAL64 + 0.5_REAL64 * lat_step

DO i=1,nlon_regrid
  lon_regrid(i) = lon0 + (i-1)*lon_step
END DO

DO i=1,nlat_regrid
  lat_regrid(i) = lat0 + (i-1)*lat_step
END DO

END SUBROUTINE create_regrid_t

!-------------------------------------------------------------------------------

SUBROUTINE create_regrid_gg()

IMPLICIT NONE

REAL(KIND=REAL64) :: lon0,lon_step
INTEGER :: i

lreduce = .FALSE.
nlon_regrid = 2*nlat_regrid

ALLOCATE(lon_regrid(nlon_regrid))
ALLOCATE(lat_regrid(nlat_regrid))

lon_step = 360.0_REAL64 / nlon_regrid
lon0 = 0.0_REAL64

DO i=1,nlon_regrid
  lon_regrid(i) = lon0 + (i-1)*lon_step
END DO

CALL get_lat(nlat_regrid,lat_regrid)

END SUBROUTINE create_regrid_gg

!-------------------------------------------------------------------------------

SUBROUTINE create_regrid_orgg()

IMPLICIT NONE

REAL(KIND=REAL64) :: lon0,lon_step
!REAL(KIND=REAL64), ALLOCATABLE :: lat(:)
REAL(KIND=C_DOUBLE), ALLOCATABLE :: lat(:)
INTEGER :: i,j,ii,jj0,jj,nloni
INTEGER(KIND=INT32) nlat

lreduce = .TRUE.
nvertex = 4
ntot_regrid = nlat_regrid*(nlat_regrid+18)
nlat = nlat_regrid

ALLOCATE(lon_regrid(ntot_regrid))
ALLOCATE(lat_regrid(ntot_regrid))
ALLOCATE(lat(0:nlat_regrid+1))
ALLOCATE(bnds_lon_regrid(nvertex,ntot_regrid))
ALLOCATE(bnds_lat_regrid(nvertex,ntot_regrid))

CALL get_lat(nlat_regrid,lat(1:nlat))
lat(0) = 2.0_REAL64*lat(1)-lat(2)
lat(nlat_regrid+1) = -lat(0)
IF (mype == 0 .AND. verbose > 0) THEN
  WRITE(outunit,*)'lat(0) = ',lat(0)
END IF

IF (write_output) OPEN(UNIT=11, FILE='output')

lon0 = 0.0_REAL64
ii=1
jj0=ntot_regrid
DO i=1,nlat_regrid/2
  nloni = 4*i + 16
  lon_step = 360.0_REAL64 / nloni
  jj0 = jj0 - nloni
  jj = jj0+1

  IF (write_output) THEN
    WRITE(11,*) 'latitude = ',lat(i),-lat(i)
    WRITE(11,*) 'number of longitudes = ',nloni
    WRITE(11,*) 'lon_step = ',lon_step
    WRITE(11,*)
  END IF

  DO j=1,nloni
    lat_regrid(ii) = lat(i)
    bnds_lat_regrid(1,ii) = 0.5_REAL64*(lat(i) + lat(i+1))
    bnds_lat_regrid(2,ii) = bnds_lat_regrid(1,ii)
    bnds_lat_regrid(3,ii) = 0.5_REAL64*(lat(i-1) + lat(i))
    bnds_lat_regrid(4,ii) = bnds_lat_regrid(3,ii)

    lon_regrid(ii) = lon0 + (j-1)*lon_step
    bnds_lon_regrid(1,ii) = lon0 + (j-1.5_REAL64)*lon_step
    bnds_lon_regrid(2,ii) = lon0 + (j-0.5_REAL64)*lon_step
    bnds_lon_regrid(3,ii) = bnds_lon_regrid(2,ii)
    bnds_lon_regrid(4,ii) = bnds_lon_regrid(1,ii)

    lat_regrid(jj) = -lat_regrid(ii)
    bnds_lat_regrid(1,jj) = -bnds_lat_regrid(4,ii)
    bnds_lat_regrid(2,jj) = -bnds_lat_regrid(3,ii)
    bnds_lat_regrid(3,jj) = -bnds_lat_regrid(2,ii)
    bnds_lat_regrid(4,jj) = -bnds_lat_regrid(1,ii)

    lon_regrid(jj) = lon_regrid(ii)
    bnds_lon_regrid(1,jj) = bnds_lon_regrid(4,ii)
    bnds_lon_regrid(2,jj) = bnds_lon_regrid(3,ii)
    bnds_lon_regrid(3,jj) = bnds_lon_regrid(2,ii)
    bnds_lon_regrid(4,jj) = bnds_lon_regrid(1,ii)

    IF (write_output) THEN
      WRITE(11,*) 'Point ',j,i,ii,' in northern hemisphere'
      WRITE(11,*) 'centre vertex = ',lon_regrid(ii),lat_regrid(ii)
      WRITE(11,*) 'boundary vertex 1 = ',bnds_lon_regrid(1,ii),bnds_lat_regrid(1,ii)
      WRITE(11,*) 'boundary vertex 2 = ',bnds_lon_regrid(2,ii),bnds_lat_regrid(2,ii)
      WRITE(11,*) 'boundary vertex 3 = ',bnds_lon_regrid(3,ii),bnds_lat_regrid(3,ii)
      WRITE(11,*) 'boundary vertex 4 = ',bnds_lon_regrid(4,ii),bnds_lat_regrid(4,ii)
      WRITE(11,*)
    
      WRITE(11,*) 'Point ',j,nlat_regrid-i+1,jj,' in southern hemisphere'
      WRITE(11,*) 'centre vertex = ',lon_regrid(jj),lat_regrid(jj)
      WRITE(11,*) 'boundary vertex 1 = ',bnds_lon_regrid(1,jj),bnds_lat_regrid(1,jj)
      WRITE(11,*) 'boundary vertex 2 = ',bnds_lon_regrid(2,jj),bnds_lat_regrid(2,jj)
      WRITE(11,*) 'boundary vertex 3 = ',bnds_lon_regrid(3,jj),bnds_lat_regrid(3,jj)
      WRITE(11,*) 'boundary vertex 4 = ',bnds_lon_regrid(4,jj),bnds_lat_regrid(4,jj)
      WRITE(11,*)
    END IF

    ii=ii+1
    jj=jj+1
  END DO
END DO

IF (write_output) CLOSE(UNIT=11)

DEALLOCATE(lat)

END SUBROUTINE create_regrid_orgg

!-------------------------------------------------------------------------------

SUBROUTINE create_regrid_orrg()

IMPLICIT NONE

REAL(KIND=REAL64) :: lon0,lat0,lon_step,lat_step
REAL(KIND=REAL64), ALLOCATABLE :: lat(:)
INTEGER :: i,j,ii,jj0,jj,nloni

lreduce = .TRUE.
nvertex = 4
ntot_regrid = nlat_regrid*(nlat_regrid+18)

ALLOCATE(lon_regrid(ntot_regrid))
ALLOCATE(lat_regrid(ntot_regrid))
ALLOCATE(lat(0:nlat_regrid+1))
ALLOCATE(bnds_lon_regrid(nvertex,ntot_regrid))
ALLOCATE(bnds_lat_regrid(nvertex,ntot_regrid))

lat_step = 180.0_REAL64 / nlat_regrid
lat0 = 90.0_REAL64 - 0.5_REAL64 * lat_step

DO i=1,nlat_regrid
  lat(i) = lat0 - (i-1)*lat_step
END DO

lat(0) = 2.0_REAL64*lat(1)-lat(2)
lat(nlat_regrid+1) = -lat(0)
IF (mype == 0 .AND. verbose > 0) THEN
  WRITE(outunit,*)'lat(0) = ',lat(0)
END IF

IF (write_output) OPEN(UNIT=11, FILE='output')

lon0 = 0.0_REAL64
ii=1
jj0=ntot_regrid
DO i=1,nlat_regrid/2
  nloni = 4*i + 16
  lon_step = 360.0_REAL64 / nloni
  jj0 = jj0 - nloni
  jj = jj0+1

  IF (write_output) THEN
    WRITE(11,*) 'latitude = ',lat(i),-lat(i)
    WRITE(11,*) 'number of longitudes = ',nloni
    WRITE(11,*)
  END IF

  DO j=1,nloni
    lat_regrid(ii) = lat(i)
    bnds_lat_regrid(1,ii) = 0.5*(lat(i) + lat(i+1))
    bnds_lat_regrid(2,ii) = bnds_lat_regrid(1,ii)
    bnds_lat_regrid(3,ii) = 0.5*(lat(i-1) + lat(i))
    bnds_lat_regrid(4,ii) = bnds_lat_regrid(3,ii)

    lon_regrid(ii) = lon0 + (j-1)*lon_step
    bnds_lon_regrid(1,ii) = lon0 + (j-1.5_REAL64)*lon_step
    bnds_lon_regrid(2,ii) = lon0 + (j-0.5_REAL64)*lon_step
    bnds_lon_regrid(3,ii) = bnds_lon_regrid(2,ii)
    bnds_lon_regrid(4,ii) = bnds_lon_regrid(1,ii)

    lat_regrid(jj) = -lat_regrid(ii)
    bnds_lat_regrid(1,jj) = -bnds_lat_regrid(4,ii)
    bnds_lat_regrid(2,jj) = -bnds_lat_regrid(3,ii)
    bnds_lat_regrid(3,jj) = -bnds_lat_regrid(2,ii)
    bnds_lat_regrid(4,jj) = -bnds_lat_regrid(1,ii)

    lon_regrid(jj) = lon_regrid(ii)
    bnds_lon_regrid(1,jj) = bnds_lon_regrid(4,ii)
    bnds_lon_regrid(2,jj) = bnds_lon_regrid(3,ii)
    bnds_lon_regrid(3,jj) = bnds_lon_regrid(2,ii)
    bnds_lon_regrid(4,jj) = bnds_lon_regrid(1,ii)

    IF (write_output) THEN
      WRITE(11,*) 'Point ',j,i,ii,' in northern hemisphere'
      WRITE(11,*) 'centre vertex = ',lon_regrid(ii),lat_regrid(ii)
      WRITE(11,*) 'boundary vertex 1 = ',bnds_lon_regrid(1,ii),bnds_lat_regrid(1,ii)
      WRITE(11,*) 'boundary vertex 2 = ',bnds_lon_regrid(2,ii),bnds_lat_regrid(2,ii)
      WRITE(11,*) 'boundary vertex 3 = ',bnds_lon_regrid(3,ii),bnds_lat_regrid(3,ii)
      WRITE(11,*) 'boundary vertex 4 = ',bnds_lon_regrid(4,ii),bnds_lat_regrid(4,ii)
      WRITE(11,*)
    
      WRITE(11,*) 'Point ',j,nlat_regrid-i+1,jj,' in southern hemisphere'
      WRITE(11,*) 'centre vertex = ',lon_regrid(jj),lat_regrid(jj)
      WRITE(11,*) 'boundary vertex 1 = ',bnds_lon_regrid(1,jj),bnds_lat_regrid(1,jj)
      WRITE(11,*) 'boundary vertex 2 = ',bnds_lon_regrid(2,jj),bnds_lat_regrid(2,jj)
      WRITE(11,*) 'boundary vertex 3 = ',bnds_lon_regrid(3,jj),bnds_lat_regrid(3,jj)
      WRITE(11,*) 'boundary vertex 4 = ',bnds_lon_regrid(4,jj),bnds_lat_regrid(4,jj)
      WRITE(11,*)
    END IF

    ii=ii+1
    jj=jj+1
  END DO
END DO

IF (write_output) CLOSE(UNIT=11)

DEALLOCATE(lat)

END SUBROUTINE create_regrid_orrg

!-------------------------------------------------------------------------------

SUBROUTINE get_lat(nlat, lats)

IMPLICIT NONE

INTEGER (KIND=INT32), INTENT(IN) :: nlat  ! Number of latitude points
                                          ! from pole to pole

REAL(KIND=C_DOUBLE), INTENT(INOUT) :: lats(:) ! Array of latitude values from
                                              ! south to north pole
!
!  Local variables
!
INTEGER(KIND=C_LONG) :: trunc ! trunc is the Gaussian number (or order)
                              ! i.e. Number of parallels between a pole
                              ! and the equator
INTEGER :: ierr               ! non-zero on error

trunc = nlat/2
! nlat must be even
IF (nlat /= 2*trunc) THEN
  ierr = 101
  WRITE(error_unit,*) 'Number of latitudes must be even nlat = ',nlat
  CALL abort()
END IF

#ifdef USE_ECCODES
ierr =  c_grib_get_gaussian_latitudes(trunc, lats)
#else
WRITE(error_unit,*) 'ECCODES library needed to calculate Gaussian latitudes'
ierr = 102
#endif
diff 
IF (ierr /= 0) THEN
  WRITE(error_unit,*) 'Error calculating Gaussian latitudes, ierr = ',ierr
  CALL abort()
END IF

! Convert lats from north-south to south-north orientation
! lats is symetrical about the equator
lats = -lats

END SUBROUTINE get_lat

END MODULE regrid_mod
