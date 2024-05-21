!=======================================================================
!  GoMDOM for FVCOM, bio-geo-chem subroutine
!=======================================================================
SUBROUTINE bio_geo_chem(B_N,T_B,S_B,IRRAD,Wind)
!------------------------------------------------------------------------------
!
!FVCOM  Define global data.
USE mod_1d, only: SPP,NTT,KBV,KBVM1,T_STEP,DELTA_D,BIO_NAME

!Definitions from other subroutines
!mod_bio_3D.F:        DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS

!GOMDOM
USE STATES
USE EUT
!USE FLAGS, ONLY : DoDroop

IMPLICIT NONE

!GOMDOM
REAL(SPP), INTENT(IN) :: IRRAD,Wind
REAL(SPP), INTENT(IN) :: T_B(KBV),S_B(KBV)
REAL(SPP), INTENT(INOUT) :: B_N(KBV,NTT)
REAL(SPP) :: PAR(KBV),x
!REAL(SPP) :: SETRATE(NTT)
INTEGER :: k,j

!Light
!PARFRAC     =  0.43_SPP
!--- Add light model
!IRRAD0=-SWRAD(I)*PARFRAC*SPRO/RAMP, from ../mod_bio_3D.F
PAR=0.
!Need to add a light model (so PAR depends on depth)
!IRRAD is input, outputs PAR
if(IRRAD.gt.tiny(x)) then
  call gd_light(B_N,S_B,IRRAD,DELTA_D,PAR)
endif
!example of previous call:
!  CALL DIATOMS(f(myi,j,k,:),DTM(i,j,k,:),T(i,j,k),PAR(i,j,k),Vol(i,j,k),dTime,i,j,k)
             DTM=0.
             !write(6,*) "T,S",T_B,S_B 
             !write(6,*) "Before WQEM",B_N(1,:)
             CALL ZOO(B_N,T_B)                !  Zooplankton kinetics
             CALL DIATOMS(B_N,T_B,PAR) !         ! Diatom kinetics
             CALL GREENS(B_N,T_B,PAR) !         ! Greens kinetics
             CALL CARBON(B_N,T_B)         ! Carbon (detritus) kinetics
             CALL PHOSPH(B_N,T_B)                      !  !         Phosphorous kinetics
             CALL NITROG(B_N,T_B)         !  !         Nitrogen kinetics
             CALL SILICA(B_N,T_B)             !  !         Silica kinetics
             CALL DISSOLVED_OXYGEN(B_N,T_B)   !         !  Dissolved Oxygen
             DO K=1,KBVM1
             ! write(6,*) "k,T_STEP,DTM",k,T_STEP,DTM(K,:)
              B_N(K,:)=B_N(K,:)+DTM(K,:)*T_STEP
             END DO !K=1,KBVM1
             !write(6,*) "Before sinking",B_N(1,:)
             call sinking(B_N)
             !write(6,*) "after sinking",B_N(1,:)

             !For printing, to check PAR 
             B_N(:,JTR) = PAR(:) 

!if(IRRAD.gt.600) then
!  DO K=1,KBVM1
!   write(6,*) "DTM",DTM(K,:)*T_STEP
!  END DO !K=1,KBVM1
!  stop
!endif             

!-----------------------------------------------------------------------------
RETURN

END SUBROUTINE bio_geo_chem
