!******************************************************************************
! PURPOSE: OUTPUT_NETCDF.F90 - Routines for outputting NCOM_GEM data to a
!          NetCDF file.
! NOTES:   Non-ADT module.
! HISTORY: 2010/04/26, Todd Plessel, plessel.todd@epa.gov, Created.
!******************************************************************************

MODULE OUTPUT_NETCDF_GD

!  USE netcdf

  IMPLICIT NONE

#ifdef _MPI
!  include 'mpif.h'
#endif

  ! Private



PUBLIC CREATE_FILE, OPEN_FILE, &
       WRITE_DATA, CLOSE_FILE, FLUSH_FILE, &
       WRITE_EXTRA_DATA, OUTPUT_NETCDF_GD_allocate

PRIVATE
CONTAINS

  ! Public

  ! Commands:

  ! CREATE_FILE: Create output NetCDF file with given header data.
  ! In a concurrent program, only one process should call this routine
  ! then call CLOSE_FILE then
  ! worker processes should call OPEN_FILE, WRITE_DATA, CLOSE_FILE.
  !
  SUBROUTINE CREATE_FILE()
    IMPLICIT NONE
    RETURN
  END SUBROUTINE CREATE_FILE



  ! OPEN_FILE: Open existing output NetCDF file for shared writing.
  ! In a concurrent program, each worker process should call
  ! OPEN_FILE, WRITE_DATA, CLOSE_FILE.
  ! FIRST_TIMESTEP is 0-based timestep number to begin writing.
  !
  SUBROUTINE OPEN_FILE()
    IMPLICIT NONE
    RETURN
  END SUBROUTINE OPEN_FILE



  ! CLOSE_FILE: Close output file.
  !
  SUBROUTINE CLOSE_FILE()
    IMPLICIT NONE
    RETURN
  END SUBROUTINE CLOSE_FILE


  ! WRITE_DATA: Write current timestep data to output file.
  ! Called by non-concurrent programs to write all data per timestep.
  !

  SUBROUTINE WRITE_DATA()
    IMPLICIT NONE
    RETURN
  END SUBROUTINE WRITE_DATA


  ! WRITE_EXTRA_DATA: Write current timestep data to output file.
  ! Called by concurrent programs to write all data per timestep
  !
  SUBROUTINE WRITE_EXTRA_DATA()
    IMPLICIT NONE
    RETURN
  END SUBROUTINE WRITE_EXTRA_DATA



  ! FLUSH_FILE: Flush buffers to disk in case of crash.
  !
  SUBROUTINE FLUSH_FILE()
    IMPLICIT NONE
    RETURN
  END SUBROUTINE FLUSH_FILE


  ! Private



Subroutine OUTPUT_NETCDF_GD_allocate
  IMPLICIT NONE
  RETURN
END Subroutine OUTPUT_NETCDF_GD_allocate

END MODULE OUTPUT_NETCDF_GD

