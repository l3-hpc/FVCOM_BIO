SUBROUTINE gd_light(f,S,Rad,dz,PAR)
!------------------------------------------------------------------------------
!FVCOM  Define global data.
USE mod_1d, only: SPP,NTT,KBV,KBVM1

USE STATES
USE EUT
USE FLAGS, ONLY: KDWD

IMPLICIT NONE

REAL(SPP), INTENT(INOUT) :: f(KBV,NTT)
REAL(SPP), INTENT(IN)    :: S(KBV),Rad,dz(KBV)
REAL(SPP), INTENT(OUT) :: PAR(KBV)
REAL(SPP) :: KESS(KBV),SAL_TERM,CHL_TERM,POC_TERM
REAL(SPP) :: IATTOP, IATBOT(KBV),OPTDEPTH,Rad_Watts
INTEGER :: k,nz

!
!------------------------------------------------------------------------------
!
!FVCOM uses umol/m2/s (at least, I think it does...need to check)
!Rad_Watts = Rad/2.77e14*PARfac
!if(Read_Solar.ne.2) then
! Rad_Watts = Rad_Watts*0.48
!endif
!if(Read_Solar.eq.2) Rad_Watts = Rad/4.57
Rad_Watts = Rad
nz = KBVM1

!write(6,*) "gd_light,Rad,KE,CCHLG,KDWD",Rad,KE,CCHLG,KDWD
!write(6,*) "nz",nz
!write(6,*) "S",S
!write(6,*) "DIA",f(:,JDIA)

!GoMDOM LIGHT MODEL, No Wind Speed
         if (KDWD .eq. 1) then
           do k = 1, nz
             KESS(k) = KE + KECHL * (f(k,JDIA) / CCHLD + f(k,JGRE) / CCHLG)
           enddo
         else
           !write(6,*) "KDWD,k",KDWD,k
           do k = 1, nz
             !write(6,*) k
             SAL_TERM = 1.084E-06 * (S(k)**4)
             !write(6,*) "sal_term",SAL_TERM 
             IF ((f(k,JDIA) + f(k,JGRE)) < 1.0E-07) THEN
                   CHL_TERM = 0.0
             ELSE
                   CHL_TERM = 0.2085 * LOG( (f(k,JDIA) * 1.0E6 / CCHLD) + &
                            & (f(k,JGRE) * 1.0E6 / CCHLG) )
             ENDIF
             !write(6,*) "chl_term",CHL_TERM 
             POC_TERM = 0.7640 * SQRT( (f(k,JLOC) * 1.0E3) +  &
                      & (f(k,JROC) * 1.0E3) + (f(k,JZOO) * 1.0E3) )
             !write(6,*) "poc_term",POC_TERM 
             KESS(k) = ( ( -0.10 * (-0.5606 - SAL_TERM + CHL_TERM + POC_TERM) ) &
                     &  + 1 ) ** (1.0/(-0.10))
             !write(6,*) "KESS",KESS(k)

             enddo
         endif

      DO k = 1,1
         IATTOP    =  Rad_Watts
         OPTDEPTH  =  KESS(k) * dz(k)
         IATBOT(k) =  IATTOP  * EXP(-OPTDEPTH)
         PAR(k)   =  (IATTOP - IATBOT(k)) / OPTDEPTH
         !write(6,*) "k,PAR",k,PAR
      END DO
      IF(k>=2) THEN
       DO k = 2,nz
         IATTOP    =  IATBOT(k-1)
         OPTDEPTH  =  KESS(k) * dz(k) 
         IATBOT(k) =  IATTOP * EXP(-OPTDEPTH)
         PAR(k) =  (IATTOP - IATBOT(k)) / OPTDEPTH
         !write(6,*) "k,PAR",k,PAR
       END DO
      ENDIF


!Check above, I think FVCOM's PAR matches WQEM
! PAR = PAR*4.57
!write(6,*) "dz",dz
!write(6,*) "KESS",KESS
!write(6,*) "OPTDEPTH",OPTDEPTH
!write(6,*) "PAR",PAR
!
!-----------------------------------------------------------------------------


RETURN

END SUBROUTINE gd_light 
