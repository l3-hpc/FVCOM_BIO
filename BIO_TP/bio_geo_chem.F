!=======================================================================
!  GoMDOM for FVCOM, bio-geo-chem subroutine
!=======================================================================
SUBROUTINE bio_geo_chem(B_N,T_B,IRRAD)
!------------------------------------------------------------------------------
!
!FVCOM  Define global data.
USE mod_1d, only: SPP,NTT,KBV,T_STEP

!Definitions from other subroutines
!mod_bio_3D.F:        DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS

IMPLICIT NONE

!GOMDOM
!INTEGER*8, INTENT(IN) :: TC_8
REAL(SPP), INTENT(IN) :: IRRAD
REAL(SPP), INTENT(IN) :: T_B(KBV)
REAL(SPP), INTENT(INOUT) :: B_N(KBV,NTT)
REAL(SPP) :: PAR
REAL(SPP) :: SETRATE(NTT)
INTEGER :: i,j,k,nz,myi

CALL NUTRIENTS

!-----------------------------------------------------------------------------
RETURN

END SUBROUTINE bio_geo_chem
