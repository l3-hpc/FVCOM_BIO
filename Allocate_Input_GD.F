Subroutine Allocate_Input_GD()

USE Model_dim
USE OUTPUT_NETCDF_GD
USE EUT
USE InRemin
USE MASS_BALANCE_GD
USE OUTPUT

IMPLICIT NONE

         nf = 23
         EXTRA_VARIABLES = 24
         STATE_VARIABLES = nf
         call OUTPUT_NETCDF_GD_allocate()
         call eut_allocate()
         call InRemin_allocate()
         call MASS_BALANCE_GD_allocate()

return
END SUBROUTINE Allocate_Input_GD
