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

!GOMDOM
USE Model_dim
USE STATES
USE INPUT_VARS, ONLY: Read_Solar,ws
USE INPUT_VARS_GD, ONLY : Which_Fluxes,Which_irradiance
USE EUT
USE Which_Flux
USE InRemin
USE FLAGS, ONLY : DoDroop

IMPLICIT NONE

!GOMDOM
!INTEGER*8, INTENT(IN) :: TC_8
REAL(SPP), INTENT(IN) :: IRRAD
REAL(SPP), INTENT(IN) :: T_B(KBV)
REAL(SPP), INTENT(INOUT) :: B_N(KBV,NTT)
REAL(SPP) :: PAR
REAL(SPP) :: SETRATE(nf)
INTEGER :: i,j,k,nz,myi
INTEGER :: icell, jcell
REAL(SPP) :: rivLoadConvFactor


!Light
!PARFRAC     =  0.43_SPP
!--- Add light model
! PAR(:,:,1) = Rad(:,:)  * PARFRAC
!IRRAD0=-SWRAD(I)*PARFRAC*SPRO/RAMP, from ../mod_bio_3D.F
  PAR = IRRAD 

!--- Add river loads.
!  river...blah

!---Start Loop...
!example of previous call:
!  CALL DIATOMS(f(myi,j,k,:),DTM(i,j,k,:),T(i,j,k),PAR(i,j,k),Vol(i,j,k),dTime,i,j,k)
!changes:
!  T_BIO defined in mod_bio_3D as T_BIO(K)=T1(I,K)
!  DTM moved to eut
!
!         DO K=1,KBV
!write(6,*) "ZOO before=",BIO_N(1,jzoo),BIO_N(19,jzoo)
             DTM=0.
             CALL ZOO(B_N,T_B)                !  Zooplankton kinetics
             CALL DIATOMS(B_N,T_B,PAR) !         ! Diatom kinetics
             CALL GREENS(B_N,T_B,PAR) !         ! Greens kinetics
             CALL CARBON(B_N,T_B)         ! Carbon (detritus) kinetics
             CALL PHOSPH(B_N,T_B)                      !  !         Phosphorous kinetics
             CALL NITROG(B_N,T_B)         !  !         Nitrogen kinetics
             CALL SILICA(B_N,T_B)             !  !         Silica kinetics
             CALL DISSOLVED_OXYGEN(B_N,T_B)   !         !  Dissolved Oxygen
             DO K=1,KBV
              B_N(K,:)=B_N(K,:)+DTM(K,:)*T_STEP
             END DO !K=1,KBV
             B_N(:,JTR) = IRRAD

!write(6,*) "T_STEP", T_STEP
!write(6,*) "ZOO after=",BIO_N(1,jzoo),BIO_N(19,jzoo)
!stop

!-----------------------------------------------------------------------------
RETURN

END SUBROUTINE bio_geo_chem
