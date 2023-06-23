
!===============================================================================

MODULE ens_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE :: util_mod
!USE :: xios

IMPLICIT NONE

INTEGER :: num_ens
INTEGER(KIND=INT32) :: ensemble_member

NAMELIST /input_ens/ num_ens

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE read_ens_nml()

IMPLICIT NONE

INTEGER :: ierr

num_ens = 0

OPEN(UNIT=nmlunit,FILE=infile)
READ(UNIT=nmlunit,NML=input_ens,IOSTAT=ierr)
IF (ierr > 0) THEN
  WRITE(error_unit,*)'Error reading input_ens namelist file, IOSTAT = ',ierr
  CALL abort()
END IF
CLOSE(UNIT=nmlunit)

END SUBROUTINE read_ens_nml

!-------------------------------------------------------------------------------

SUBROUTINE get_ensemble_member_number()

IMPLICIT NONE

INTEGER :: ierr

CHARACTER(LEN=32) :: c_ensmember

CALL GET_ENVIRONMENT_VARIABLE('ENS_MEMBER',c_ensmember,STATUS=ierr)
IF (ierr > 0) THEN
  WRITE(error_unit,*) 'Error reading ENS_MEMBER envionment variable'
  WRITE(error_unit,*) 'Error status = ',ierr
  CALL abort()
END IF

READ(c_ensmember, '(I5)') ensemble_member

END SUBROUTINE get_ensemble_member_number

!-------------------------------------------------------------------------------

SUBROUTINE xios_def_ens()

IMPLICIT NONE

INTEGER(KIND=INT32) :: n = 1
INTEGER(KIND=INT32) :: prec = 8

!CALL xios_set_axis_attr("ensemble", n_glo=num_ens, begin=ensemble_member, n=n, prec=prec)

END SUBROUTINE xios_def_ens

END MODULE ens_mod
