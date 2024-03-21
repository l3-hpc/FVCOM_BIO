Module INPUT_VARS_GD

USE Model_dim, ONLY: nf

IMPLICIT NONE

!--Run Specifics---------------
real*8 DayS_8, DayE_8
!--Switches in GEM---------
integer Which_Fluxes(9)
integer Which_irradiance
!----River Params-----------------------------------
real rcNO3
real rcNH4
real rcPO4
real rcSi
!----Other including Boundary Conditions------------
real KH_coeff 
integer Which_Outer_BC
real m_OM_init,m_OM_bc,m_OM_sh

END MODULE INPUT_VARS_GD
