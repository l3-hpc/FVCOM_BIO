!=======================================================================
!  GoMDOM for FVCOM, bio-geo-chem subroutine
!=======================================================================
SUBROUTINE bio_geo_chem(B_N,T_B,S_B,IRRAD,Wind)
!------------------------------------------------------------------------------
!
!FVCOM  Define global data.
USE mod_1d, only: SPP,NTT,KBV,W_N,N_0

!Definitions from other subroutines
!mod_bio_3D.F:        DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS

IMPLICIT NONE

!GOMDOM
REAL(SPP), INTENT(IN) :: IRRAD,Wind
REAL(SPP), INTENT(IN) :: T_B(KBV),S_B(KBV)
REAL(SPP), INTENT(INOUT) :: B_N(KBV,NTT)

call sinking(B_N)

!-----------------------------------------------------------------------------
RETURN

END SUBROUTINE bio_geo_chem
