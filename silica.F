SUBROUTINE SILICA(BIO_N,T_BIO)
!------------------------------------------------------------------------------
!-   Purpose and Methods : Kinetic sources and sinks of silica are
!-                         computed in this subroutine.
!-
!-   Created:    July 28, 1992
!-   Authors:    Carl F. Cerco      : Water quality scheme
!-               Raymond S. Chapman : Numerical solution scheme
!-               Thomas M. Cole     : Computer algorithms & coding
!-               Hydroqual          : Sediment compartment
!-
!-               Water Quality Modeling Group
!-               U.S. Army Corps of Engineers
!-               Waterways Experiment Station
!-               Vicksburg, Mississippi 39180
!------------------------------------------------------------------------------
!
USE mod_1D, ONLY: SPP,KBV,NTT
USE STATES
USE EUT

IMPLICIT NONE

!
!------------------------------------------------------------------------------
!
REAL(SPP), INTENT(IN) :: BIO_N(KBV,NTT),T_BIO(KBV)
INTEGER :: K

REAL(SPP) TEMP


REAL(SPP) :: SU                         ! Unavailable silica
REAL(SPP) :: SAD                        ! Silica in diatoms

REAL(SPP) :: PSD                        ! Change in particulate silica dissolution
REAL(SPP) :: DAS                        ! Change in available silica
REAL(SPP) :: PUS                        ! Change in particulate unavailable silica
                                   ! metabolism and predation

DO K=1,KBV

!
!------------------------------------------------------------------------------
!
TEMP = T_BIO(K)

!------------------------------------------------------------------------------
! Particulate biogenic silica
!------------------------------------------------------------------------------
   SU = BIO_N(K,JSU)   

!------------------------------------------------------------------------------
! Diatom silica
!------------------------------------------------------------------------------
   SAD = ASCD * BIO_N(K,JDIA)                            

!------------------------------------------------------------------------------
! Particulate biogenic silica dissolution
!------------------------------------------------------------------------------
   PSD = KSUA * SU * EXP(KTSUA * (TEMP - TRSUA))        

!------------------------------------------------------------------------------
! Available silica: predation/production
!------------------------------------------------------------------------------
   DAS = FSAP * PRD(K) * ASCD - (PD(K) * SAD)             

!------------------------------------------------------------------------------
! Unavailable silica: metabolism/predation
!------------------------------------------------------------------------------
   PUS = BMD(K) * SAD + (1.0 - FSAP) * PRD(K) * ASCD      

!------------------------------------------------------------------------------
! Dissolved available silica mass derivative
!------------------------------------------------------------------------------
   DTM(K,JSA)  = DTM(K,JSA) + PSD+DAS  

!------------------------------------------------------------------------------
! Particulate biogenic silica mass derivative
!------------------------------------------------------------------------------
   DTM(K,JSU) = DTM(K,JSU) + PUS - PSD   

!
!------------------------------------------------------------------------------
!
ENDDO

RETURN

END SUBROUTINE SILICA

