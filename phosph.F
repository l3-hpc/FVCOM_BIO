SUBROUTINE PHOSPH(BIO_N,T_BIO)
!------------------------------------------------------------------------------
!-
!-   $Id: phosph.F90,v 1.0.6.1 2014/08/26 22:54:04 wmelende Exp wmelende $
!-   $Locker: wmelende $
!-
!-   Purpose and Methods : Kinetic sources and sinks of phosphorus are
!-                         computed in this subroutine.
!-
!-   Inputs  :   None
!-   Outputs :   None
!-   Controls:   None
!-
!-   Calls:      REPORT
!-   Called by:  EUTRO
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
!-
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
REAL(SPP) :: ALGP              ! Algal phosphate                         
REAL(SPP) :: DOP               ! Dissolved organic P mass derivative
REAL(SPP) :: POPL              ! Labile detrital P produced by
                          ! algal metabolism and predation
REAL(SPP) :: POPR              ! Refractory detrital P produced by
                          ! algal metabolism and predation

REAL(SPP) :: PO4D              ! Change in phosphate from diatom source/sink
REAL(SPP) :: PO4G              ! Change in phosphate from greens source/sink
REAL(SPP) :: MNLDOP            ! Dissolved P mineralization mass derivative
REAL(SPP) :: HDRLPOP           ! Labile detrital P hydrolysis mass derivative
REAL(SPP) :: HDRRPOP           ! Refractory detrital P hydrolysis mass derivative

REAL(SPP) :: ZDEATH            ! Zooplankton mortality concentration
REAL(SPP) :: PO4ZOO            ! Phosphate released by zoo. mortality
REAL(SPP) :: DOPZOO            ! Dissolved organic P released by zoo. mortality
REAL(SPP) :: LPOPZOO           ! Labile organic P released by zoo. mortality
REAL(SPP) :: RPOPZOO           ! Refractory organic P released by zoo. mortality


REAL(SPP) :: KDOP              ! Dissolved organic phosphorous mineralization rate [s^-1]
REAL(SPP) :: KPOPL             ! Labile organic phosphorous hydrolysis rate [s^-1]
REAL(SPP) :: KPOPR             ! Labile organic phosphorous hydrolysis rate [s^-1]

REAL(SPP) :: FRAC_SRP          ! Fraction of phosphorous that is in the silica
REAL(SPP) :: FRAC_DOP          ! Fraction of phosphorous that is dissolved

REAL(SPP) :: ALG               ! Total phytoplankton concentration (dia+gre)


!
!------------------------------------------------------------------------------
!
DO K=1,KBV

 
!
!  Contribution from diatoms
!
   ALGP = APCP * BIO_N(K,JDIA)
   IF (BIO_N(K,JSRP) > 0.0) THEN
       IF (AVFRAC > 0.0) THEN
          FRAC_SRP = BIO_N(K,JSRP)/(BIO_N(K,JSRP) + AVFRAC*BIO_N(K,JDOP))
       ELSE
          FRAC_SRP = 1.0
       ENDIF
   ELSE
       FRAC_SRP = 0.0
   ENDIF
   PO4D = (FPID*BMD(K)-PD(K)*FRAC_SRP)*ALGP + &
           FPIP*(1.0-GREFF)*PRD(K)*APCP
   IF (AVFRAC*BIO_N(K,JDOP) > 0.0)THEN
      FRAC_DOP = AVFRAC*BIO_N(K,JDOP)/(BIO_N(K,JSRP) + AVFRAC*BIO_N(K,JDOP))
   ELSE
      FRAC_DOP = 0.0
   ENDIF
   DOP = (FPDD*BMD(K) - PD(K)*FRAC_DOP)*ALGP + & 
         FPDP*(1.0-GREFF)*PRD(K)*APCP
   POPL = FPLD*BMD(K)*ALGP + FPLP*(1.0-GREFF)*PRD(K)*APCP
   POPR = FPRD*BMD(K)*ALGP + FPRP*(1.0-GREFF)*PRD(K)*APCP
   DTM(K,JSRP) = DTM(K,JSRP) + PO4D 
   DTM(K,JDOP) = DTM(K,JDOP) + DOP 
   DTM(K,JLOP) = DTM(K,JLOP) + POPL 
   DTM(K,JROP) = DTM(K,JROP) + POPR

!
!  Contribution from greens
!
   ALGP = APCP * BIO_N(K,JGRE)
   IF (BIO_N(K,JSRP) > 0.0) THEN
      IF (AVFRAC > 0.0) THEN
         FRAC_SRP = BIO_N(K,JSRP)/(BIO_N(K,JSRP) + AVFRAC*BIO_N(K,JDOP))
      ELSE
         FRAC_SRP = 1.0
      ENDIF
   ELSE
      FRAC_SRP = 0.0
   ENDIF
   PO4G = (FPIG*BMG(K)-PG(K)*FRAC_SRP)*ALGP + & 
           FPIP*(1.0-GREFF)*PRG(K)*APCP
   IF (AVFRAC*BIO_N(K,JDOP) > 0.0) THEN
      FRAC_DOP = AVFRAC*BIO_N(K,JDOP)/(BIO_N(K,JSRP) + AVFRAC*BIO_N(K,JDOP))
   ELSE
      FRAC_DOP = 0.0
   ENDIF
   DOP  = (FPDG*BMG(K) - PG(K)*FRAC_DOP)*ALGP + &
           FPDP*(1.0-GREFF)*PRG(K)*APCP
   POPL = FPLG*BMG(K)*ALGP + FPLP*(1.0-GREFF)*PRG(K)*APCP
   POPR = FPRG*BMG(K)*ALGP + FPRP*(1.0-GREFF)*PRG(K)*APCP
   DTM(K,JSRP) = DTM(K,JSRP) + PO4G
   DTM(K,JDOP) = DTM(K,JDOP) + DOP 
   DTM(K,JLOP) = DTM(K,JLOP) + POPL
   DTM(K,JROP) = DTM(K,JROP) + POPR
                                                            
!------------------------------------------------------------------------------
!  Zooplankton 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Zooplankton mortality: it was changed from a linear to quadratic dependence
! on the zooplankton concentration. See the zoo.F90 subroutine for more
! details about this change.
! [ZDTH] = 1/(kg/m^3 s)
!------------------------------------------------------------------------------
   ZDEATH = ZDTH * BIO_N(K,JZOO) * BIO_N(K,JZOO)  

!------------------------------------------------------------------------------
! PO4 component
!------------------------------------------------------------------------------
   PO4ZOO  = FPIZ*APCP*ZDEATH   

!------------------------------------------------------------------------------
! DOP component
!------------------------------------------------------------------------------
   DOPZOO  = FPDZ*APCP*ZDEATH                                    

!------------------------------------------------------------------------------
! Labile particulate component
!------------------------------------------------------------------------------
   LPOPZOO = FPLZ*APCP*ZDEATH                                    

!------------------------------------------------------------------------------
!  Refractory particulate component
!------------------------------------------------------------------------------
   RPOPZOO = FPRZ*APCP*ZDEATH                               

!------------------------------------------------------------------------------
!  PO4 time derivative
!------------------------------------------------------------------------------
   DTM(K,JSRP) = DTM(K,JSRP) + PO4ZOO   

!------------------------------------------------------------------------------
!  DOP time derivative
!------------------------------------------------------------------------------
   DTM(K,JDOP) = DTM(K,JDOP) + DOPZOO             

!------------------------------------------------------------------------------
!  Labile particulate time derivative
!------------------------------------------------------------------------------
   DTM(K,JLOP) = DTM(K,JLOP) + LPOPZOO            

!------------------------------------------------------------------------------
!  Refractory particulate time derivative
!------------------------------------------------------------------------------
   DTM(K,JROP) = DTM(K,JROP) + RPOPZOO         

!------------------------------------------------------------------------------
!  Mineralization and hydrolysis
!------------------------------------------------------------------------------

   ALG = BIO_N(K,JGRE)+BIO_N(K,JDIA)
   
   KDOP  = KDP + KDPALG*ALG*KHP/(KHP+BIO_N(K,JSRP))
   KPOPL = KLP + KLPALG*ALG*KHP/(KHP+BIO_N(K,JSRP))
   KPOPR = KRP + KRPALG*ALG*KHP/(KHP+BIO_N(K,JSRP))
   
   MNLDOP  = KDOP*FTMNL(K)*BIO_N(K,JDOP)
   HDRLPOP = KPOPL*FTHDR(K)*BIO_N(K,JLOP)
   HDRRPOP = KPOPR*FTHDR(K)*BIO_N(K,JROP)
   DTM(K,JSRP) = DTM(K,JSRP) + MNLDOP 
   DTM(K,JDOP) = DTM(K,JDOP) + HDRRPOP + HDRLPOP - MNLDOP
   DTM(K,JLOP) = DTM(K,JLOP) - HDRLPOP 
   DTM(K,JROP) = DTM(K,JROP) - HDRRPOP 

!
!------------------------------------------------------------------------------
!
ENDDO

RETURN

END SUBROUTINE PHOSPH

