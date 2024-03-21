SUBROUTINE ZOO(BIO_N,T_BIO)
!DTM moved to eut
!------------------------------------------------------------------------------
!- Created: 03/05/99  J. Pauer & M. Settles
!- Revised:  09/02/14  W. Melendez :    Added parallelization code.
!------------------------------------------------------------------------------
USE mod_1d, ONLY: SPP,KBV,NTT !T_BIO(K)=T1(I,K), BIO_N(K,J)
USE STATES
USE EUT

IMPLICIT NONE
!
!------------------------------------------------------------------------------
!
REAL(SPP), INTENT(IN) :: BIO_N(KBV,NTT),T_BIO(KBV)
INTEGER :: K
REAL(SPP)  :: ZPLANK     ! Zooplankton concentration
REAL(SPP)  :: PHYTO      ! Phytoplankton concentration 
REAL(SPP)  :: ZGRAZ      ! Factor in zoo grazing rate 
REAL(SPP)  :: ZGRTH      ! Growth rate of zooplankton
REAL(SPP)  :: TEMP


DO K=1,KBV

!
!------------------------------------------------------------------------------
!
TEMP = T_BIO(K)

!------------------------------------------------------------------------------
! Zoo concentration
!------------------------------------------------------------------------------
    ZPLANK = BIO_N(K,JZOO)
                                  
!------------------------------------------------------------------------------
! Total phytoplankton concentration
!------------------------------------------------------------------------------
    PHYTO  = BIO_N(K,JDIA) + BIO_N(K,JGRE)    
               
    IF (PHYTO > 1.0E-20_SPP) THEN

!------------------------------------------------------------------------------
! Zooplankton grazing factor
!------------------------------------------------------------------------------
        ZGRAZ = PHYTO / (KSZ +  PHYTO) * CGZ * ZTHET**(TEMP - TZREF) &
              * ZPLANK

!------------------------------------------------------------------------------
! Diatom predation rate
!------------------------------------------------------------------------------
        PRD(k) = ZGRAZ * BIO_N(K,JDIA) / PHYTO                

!------------------------------------------------------------------------------
! Greens predation rate
!------------------------------------------------------------------------------
        PRG(k) = ZGRAZ * BIO_N(K,JGRE) / PHYTO  
              
    ELSE    

        ZGRAZ     = 0.0_SPP
        PRD(k) = 0.0_SPP
        PRG(k) = 0.0_SPP

    ENDIF

!------------------------------------------------------------------------------
! Zooplankton "growth"
! Higher predation (mortality) term in the zooplankton kinetics was 
! changed from a linear to a quadratic dependence on the concentration.
! Phil DePetro proposed the following equation for higher order predation
!
!     ZOOGRAZ = (ZDTH * (ZPLANK * 10^6)^2) / 10^6
!
! where ZDTH has units of 1/(mg/m^3 s).
! 
! The factor of 10^6 that is multipliying the zooplankton concentration 
! converts the units from kg/m^3 to mg/m^3. The factor of 10^6 in the 
! denominator converts the units back to kg/m^3. All this unit conversion 
! scheme involving the 10^6 conversion factor looks too arbitrary and from
! a mathematical standpoint it is simply acting as a type of coefficient 
! that scales up the mortality rate.  It's kind of hard to justify a unit
! conversion to mg/m^3 in the equation above and not some other units such 
! as ug/m^3. The main drive behind the use of this particular unit converstion 
! is that it makes the zooplankton population more stable and that it prevents 
! the phytoplankton from undergoing boom and bust cycles. With this in mind, 
! we decided to take a closer look at the equation above to see whether it 
! could be rewritten in a way that could be seen as being less dependent on the
! choice of unit conversions, and we found a way of doing it. The basic idea
! was to make the 10^6 conversion factor "disappear" from the equation in
! order to give it a cleaner look.  This was accomplished by noticing the 
! following:
!
!     (ZPLANK * 10^6)^2)/10^6 = ZPLANK * (10^6)^2/10^6 = ZPLANK * 10^6
!          
! Grouping 10^6 with ZDTH we get that
!                   
!     ZOOGRAZ = (ZDTH * 10^6) * ZPLANK^2
!
! This is equivalent to 
!
!     ZOOGRAZ = ZDTH * ZPLANK^2
!
! where ZDTH now has units of 1/(kg/m^3 s).  This means that in the input 
! deck users will have to specify ZDTH values that are much larger than 
! the ones that they used to input before.  For example: instead of setting 
! ZDTH to 1.0E-07, it should now be set to 1.0E-01 to be consistent with the 
! units of 1/(kg/m^3 s).
!
! See Cerco and Noel (2002) for more on details in using a quadratic form
! for zooplankton higher predation.
!
!------------------------------------------------------------------------------
    ZGRTH = (GREFF * ZGRAZ) - (ZDTH * ZPLANK * ZPLANK)                      

#    if defined (DEBUG)
write(6,*) "In ZOO"
write(6,*) "TEMP=",TEMP
write(6,*) "k,DTM=",k,DTM(k,JZOO)
write(6,*) "ZPLANK=",ZPLANK
write(6,*) "PHYTO=",PHYTO
write(6,*) "ZGRAZ=",ZGRAZ
write(6,*) "endif"
write(6,*) "ZGRTH=",ZGRTH
write(6,*) "k,DTM",k,DTM(k,JZOO)
write(6,*) "exit ZOO"
#    endif

!------------------------------------------------------------------------------
! Change in zooplankton mass
!------------------------------------------------------------------------------
    DTM(K,JZOO) = DTM(K,JZOO) +  ZGRTH  

ENDDO
!
!------------------------------------------------------------------------------
!

RETURN

END SUBROUTINE ZOO
