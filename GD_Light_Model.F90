SUBROUTINE GD_Light_Model(f,S,Rad,PAR,dz)
!------------------------------------------------------------------------------
!
!Input: Rad in Watts...actually, just convert it to Watts for this one
!so Rad is in photons/cm2/s or whatever.

USE Model_dim
USE STATES
USE EUT
USE INPUT_VARS, ONLY:Read_Solar
USE LIGHT_VARS, ONLY:PARfac
USE FLAGS, ONLY: KDWD

IMPLICIT NONE

REAL, INTENT(INOUT) :: f(0:(myim+1),jm,km,nf)
REAL, INTENT(IN)    :: S(im,jm,km),Rad(im,jm)
REAL, INTENT(IN) :: dz(im,jm,km)
REAL, INTENT(OUT) :: PAR(im,jm,km)
REAL :: KESS(km),SAL_TERM,CHL_TERM,POC_TERM
REAL :: IATTOP, IATBOT(im,jm,km),OPTDEPTH,Rad_Watts(im,jm)
INTEGER :: i,j,k,nz,myi

!
!------------------------------------------------------------------------------
!

Rad_Watts = Rad/2.77e14*PARfac

if(Read_Solar.ne.2) then
 Rad_Watts = Rad_Watts*0.48
endif

if(Read_Solar.eq.2) Rad_Watts = Rad/4.57

!GoMDOM LIGHT MODEL, No Wind Speed
 do j = 1,jm
     myi = 1
     do i = myi_start,myi_end
        nz = nza(i,j)
       if(nz.gt.0) then
         if (KDWD .eq. 1) then
           do k = 1, nz
             KESS(k) = KE + KECHL * (f(myi,j,k,JDIA) / CCHLD + f(myi,j,k,JGRE) / CCHLG)
           enddo
         else
           do k = 1, nz
             SAL_TERM = 1.084E-06 * (S(i,j,k)**4)
    
             IF ((f(myi,j,k,JDIA) + f(myi,j,k,JGRE)) < 1.0E-07) THEN
                   CHL_TERM = 0.0
             ELSE
                   CHL_TERM = 0.2085 * LOG( (f(myi,j,k,JDIA) * 1.0E6 / CCHLD) + &
                            & (f(myi,j,k,JGRE) * 1.0E6 / CCHLG) )
             ENDIF
    
             POC_TERM = 0.7640 * SQRT( (f(myi,j,k,JLOC) * 1.0E3) +  &
                      & (f(myi,j,k,JROC) * 1.0E3) + (f(myi,j,k,JZOO) * 1.0E3) )
    
             KESS(k) = ( ( -0.10 * (-0.5606 - SAL_TERM + CHL_TERM + POC_TERM) ) &
                     &  + 1 ) ** (1.0/(-0.10))
           enddo
         endif

      DO k = 1,1
         IATTOP    =  Rad_Watts(i,j)
         OPTDEPTH  =  KESS(k) * dz(i,j,k)
         IATBOT(i,j,k) =  IATTOP  * EXP(-OPTDEPTH)
         PAR(i,j,k)   =  (IATTOP - IATBOT(i,j,k)) / OPTDEPTH
      END DO
      IF(k>=2) THEN
       DO k = 2,nz
         IATTOP    =  IATBOT(i,j,k-1)
         OPTDEPTH  =  KESS(k) * dz(i,j,k) 
         IATBOT(i,j,k) =  IATTOP * EXP(-OPTDEPTH)
         PAR(i,j,k) =  (IATTOP - IATBOT(i,j,k)) / OPTDEPTH
       END DO
      ENDIF
       endif !End of if(nza(i,j) statement
       myi = myi + 1
   enddo      ! end of do i block do loop
 enddo      ! end of do j block do loop


 PAR = PAR*4.57

!
!-----------------------------------------------------------------------------


RETURN

END SUBROUTINE GD_Light_Model 
