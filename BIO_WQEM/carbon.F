SUBROUTINE CARBON(BIO_N,T_BIO)
!------------------------------------------------------------------------------
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
!-
!- Revised: 09/02/14  W. Melendez :  Added parallelization code.
!------------------------------------------------------------------------------
!
USE mod_1D, ONLY: SPP,KBV,KBVM1,NTT
USE STATES
USE EUT
USE FLAGS, ONLY: DO_DO2

IMPLICIT NONE

!
!------------------------------------------------------------------------------
!
REAL(SPP), INTENT(IN) :: BIO_N(KBV,NTT),T_BIO(KBV)
INTEGER :: K
REAL(SPP) TEMP

REAL(SPP) :: DOCMET            ! DOC production by algal metabolism
REAL(SPP) :: LPOCMET           ! LPOC production by algal metabolism
REAL(SPP) :: RPOCMET           ! RPOC production by algal metabolism
REAL(SPP) :: DOCPRD            ! DOC production by algal predation
REAL(SPP) :: RPDCPRD           ! Refractory detritus production by algal predation
REAL(SPP) :: LPDCPRD           ! Labile detritus production by algal predation


REAL(SPP) :: KDOC              ! Maximum mineralization rate of DOC [sec^-1]
REAL(SPP) :: KRPOC             ! Maximum hydrolysis rate of refractory POC [sec^-1]
REAL(SPP) :: KLPOC             ! Maximum hydrolysis rate of labile POC [sec^-1]

REAL(SPP) :: MNLDOC            ! DOC mineralization mass derivative
REAL(SPP) :: HDRRPOC           ! Refractory PDC hydrolysis mass derivative
REAL(SPP) :: HDRLPOC           ! Labile PDC hydrolysis mass derivative

REAL(SPP) :: DOCZOO            ! DOC produced by zooplankton mortality
REAL(SPP) :: LPOCZOO           ! Labile detritus produced by zooplankton mortality
REAL(SPP) :: RPOCZOO           ! Refractory detritus produced by zooplankton 
                          ! mortality

REAL(SPP) :: ZDEATH            ! Zooplankton mortality concentration

REAL(SPP) :: DIAT              ! Diatom carbon concentration
REAL(SPP) :: GRE               ! Green algae carbon concentration
REAL(SPP) :: ALG               ! Phytoplankton concentration (dia+gre)

!------------------------------------------------------------------------------
! Denitrification rate in units of kg/(m^3 * s)
!------------------------------------------------------------------------------
REAL(SPP) :: DENITRIFICATION_C   

DO K=1,KBVM1

!
!------------------------------------------------------------------------------
!
TEMP=T_BIO(K)
!------------------------------------------------------------------------------
!  Diatoms
!------------------------------------------------------------------------------
    
!------------------------------------------------------------------------------
! Diatom carbon concentration
!------------------------------------------------------------------------------
   DIAT = BIO_N(K,JDIA)                         

!------------------------------------------------------------------------------
! DOC from diatom metabolism
!------------------------------------------------------------------------------
   DOCMET  = FCDD * BMD(K) * DIAT           

!------------------------------------------------------------------------------
! LPOC from diatom metabolism
!------------------------------------------------------------------------------
   LPOCMET = FCLD * BMD(K) * DIAT           

!------------------------------------------------------------------------------
! RPOC from diatom metabolism
!------------------------------------------------------------------------------
   RPOCMET = FCRD * BMD(K) * DIAT           

!------------------------------------------------------------------------------
! DOC from diatom predation
!------------------------------------------------------------------------------
   DOCPRD  = FCDP * (1.0 - GREFF) * PRD(K)  

!------------------------------------------------------------------------------
! Refractory particulate from diatom predation
!------------------------------------------------------------------------------
   RPDCPRD = FCRP * (1.0 - GREFF) * PRD(K)  
                                               
!------------------------------------------------------------------------------
! Labile particulate from diatom predation 
!------------------------------------------------------------------------------
   LPDCPRD = FCLP * (1.0 - GREFF) * PRD(K)  
                                               
!------------------------------------------------------------------------------
!  A low Dissolved Oxygen (DO) concentration results in no release of CO2 
!  during diatoms metabolism or predation, but higher DOC concentration. A 
!  low DO concentration is considered anything less than 2 mg/L.
!  Unit conversion: 2 mg/L is equal to 2 * 1.0E-03 kg/m^3
!------------------------------------------------------------------------------
   IF (DO_DO2) THEN

       IF (BIO_N(K,JDO2) <= (2.0 * 1.0E-03)) THEN

           DOCMET = DOCMET + (1.0 - (FCDD + FCLD + FCRD)) * BMD(K) * DIAT

           DOCPRD = DOCPRD + (1.0 - GREFF) * (1.0 - (FCDP + FCLP + FCRP)) * &
                   & PRD(K)

       ENDIF 

   ENDIF

           !SUM_DOCMET(k) = SUM_DOCMET(k) + DOCMET
           !SUM_DOCPRD(k) = SUM_DOCPRD(k) + DOCPRD


!------------------------------------------------------------------------------
!         DOC mass derivative
!------------------------------------------------------------------------------
   DTM(K,JDOC) = DTM(K,JDOC) + DOCMET + DOCPRD 

!------------------------------------------------------------------------------
!         Refractory detritus mass derivative
!------------------------------------------------------------------------------
   DTM(K,JROC) = DTM(K,JROC) + RPOCMET + RPDCPRD 

!------------------------------------------------------------------------------
!         Labile detritus mass derivative 
!------------------------------------------------------------------------------
   DTM(K,JLOC) = DTM(K,JLOC) + LPOCMET + LPDCPRD 

! END DO                                                               


!------------------------------------------------------------------------------
!  Greens
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Green algae carbon concentation
!------------------------------------------------------------------------------
   GRE = BIO_N(K,JGRE)    

!------------------------------------------------------------------------------
!  DOC from green algae metabolism
!------------------------------------------------------------------------------
   DOCMET = FCDG * BMG(K) * GRE           

!------------------------------------------------------------------------------
! LPOC from green algae metabolism
!------------------------------------------------------------------------------
   LPOCMET = FCLG * BMG(K) * GRE          

!------------------------------------------------------------------------------
! RPOC from green algae metabolism
!------------------------------------------------------------------------------
   RPOCMET = FCRG * BMG(K) * GRE          

!------------------------------------------------------------------------------
! DOC from green algae predation
!------------------------------------------------------------------------------
   DOCPRD = FCDP * (1.0 - GREFF) * PRG(K)   

!------------------------------------------------------------------------------
! Refractory particulate from green algae predation
!------------------------------------------------------------------------------
   RPDCPRD = FCRP * (1.0 - GREFF) * PRG(K) 

!------------------------------------------------------------------------------
! Labile particulate from green algae predation 
!------------------------------------------------------------------------------
   LPDCPRD = FCLP * (1.0 - GREFF) * PRG(K)  

!------------------------------------------------------------------------------
!  A low Dissolved Oxygen (DO) concentration results in no release of CO2 
!  during greens metabolism or predation, but higher DOC concentration. A 
!  low DO concentration is considered anything less than 2 mg/L.
!  Unit conversion: 2 mg/L is equal to 2 * 1.0E-03 kg/m^3
!------------------------------------------------------------------------------
   IF (DO_DO2) THEN

       IF (BIO_N(K,JDO2) <= (2.0 * 1.0E-03)) THEN

           DOCMET = DOCMET + (1.0 - (FCDG + FCLG + FCRG)) * BMG(K) * GRE

           DOCPRD = DOCPRD + (1.0 - GREFF) * (1.0 - (FCDP + FCLP + FCRP)) * &
                  & PRG(K)

       ENDIF 

   ENDIF 

           !SUM_DOCMET(k) = SUM_DOCMET(k) + DOCMET
           !SUM_DOCPRD(k) = SUM_DOCPRD(k) + DOCPRD


!------------------------------------------------------------------------------
!         DOC mass derivative
!------------------------------------------------------------------------------
   DTM(K,JDOC) = DTM(K,JDOC) + DOCMET + DOCPRD      

!------------------------------------------------------------------------------
!         Refractory detritus mass derivative
!------------------------------------------------------------------------------
   DTM(K,JROC) = DTM(K,JROC) + RPOCMET + RPDCPRD 

!------------------------------------------------------------------------------
!         Labile detritus mass derivative 
!------------------------------------------------------------------------------
   DTM(K,JLOC) = DTM(K,JLOC) + LPOCMET + LPDCPRD  
         

!------------------------------------------------------------------------------
! Zooplankton mortality
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! ZDEATH --> Zooplankton mortality rate; it was changed from a linear to a
!            quadratic dependence on the zooplankton concentration.  See
!            zoo.F90 subroutine for more details about this change.
!            [ZDTH] = 1/(kg/m^3 s)
!------------------------------------------------------------------------------
   ZDEATH = ZDTH * BIO_N(K,JZOO) * BIO_N(K,JZOO) 

!------------------------------------------------------------------------------
!  DOC component
!------------------------------------------------------------------------------
   DOCZOO  = FCDZ * ZDEATH               

   IF (DO_DO2) THEN

       IF (BIO_N(K,JDO2) <= (2.0 * 1.0E-03)) THEN

           DOCZOO = DOCZOO + (1.0 - (FCDZ + FCLZ + FCRZ)) * ZDEATH

       ENDIF 

   ENDIF

           !SUM_DOCZOO(k) = SUM_DOCZOO(k) + DOCZOO

!------------------------------------------------------------------------------
!  Labile particulate component
!----------!--------------------------------------------------------------------
   LPOCZOO = FCLZ * ZDEATH               

!------------------------------------------------------------------------------
!  Refractory particulate component
!------------------------------------------------------------------------------
   RPOCZOO = FCRZ * ZDEATH               

!------------------------------------------------------------------------------
! DOC time derivative
!------------------------------------------------------------------------------
   DTM(K,JDOC) = DTM(K,JDOC) + DOCZOO   

!------------------------------------------------------------------------------
! Labile particulate time derivative
!------------------------------------------------------------------------------
   DTM(K,JLOC) = DTM(K,JLOC) + LPOCZOO

!------------------------------------------------------------------------------
! Refractory particulate time derivative
!------------------------------------------------------------------------------
   DTM(K,JROC) = DTM(K,JROC) + RPOCZOO  
 

!------------------------------------------------------------------------------
!  Mineralization
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Total phytoplankton concentration
!------------------------------------------------------------------------------
   ALG = BIO_N(K,JGRE) + BIO_N(K,JDIA)

!------------------------------------------------------------------------------
!  Mineralization temperature dependence
!------------------------------------------------------------------------------
   FTMNL(K) = EXP( KTMNL * (TEMP - TRMNL) )   

!------------------------------------------------------------------------------
!  Hydrolysis temperature dependence
!------------------------------------------------------------------------------
   FTHDR(K) = EXP( KTHDR * (TEMP - TRHDR) )   

!------------------------------------------------------------------------------
!  Maximum mineralization rate
!------------------------------------------------------------------------------
   KDOC = KDC + KDCALG * ALG                  

!------------------------------------------------------------------------------
!  Maximum refractory hydrolysis rate
!------------------------------------------------------------------------------
   KRPOC = KRC + KRCALG * ALG           

!------------------------------------------------------------------------------
!  Maximum labile hydrolysis rate 
!------------------------------------------------------------------------------
   KLPOC = KLC + KLCALG * ALG           

!------------------------------------------------------------------------------
!  Mineralization time derivative
!------------------------------------------------------------------------------
   IF (DO_DO2) THEN 
       MNLDOC = KDOC * FTMNL(K) * BIO_N(K,JDOC) *     &
              & ( BIO_N(K,JDO2) / (KHODOC + BIO_N(K,JDO2)) )
   ELSE
       MNLDOC = KDOC * FTMNL(K) * BIO_N(K,JDOC)   
   ENDIF

!------------------------------------------------------------------------------
!  Refractory hydrolysis time derivative
!------------------------------------------------------------------------------
   HDRRPOC = KRPOC * FTHDR(K) * BIO_N(K,JROC)  

!------------------------------------------------------------------------------
!  Labile hydrolysis time derivative 
!------------------------------------------------------------------------------
   HDRLPOC = KLPOC * FTHDR(K) * BIO_N(K,JLOC)  
                      
!------------------------------------------------------------------------------
!  DOC time derivative 
!------------------------------------------------------------------------------
   DTM(K,JDOC)  = DTM(K,JDOC) + HDRRPOC + HDRLPOC - MNLDOC

!------------------------------------------------------------------------------
!  Refractory particulate mass derivative
!------------------------------------------------------------------------------
   DTM(K,JROC) = DTM(K,JROC) - HDRRPOC   

!------------------------------------------------------------------------------
!  Labile particulate mass derivative
!------------------------------------------------------------------------------
   DTM(K,JLOC) = DTM(K,JLOC) - HDRLPOC 

! END DO


!------------------------------------------------------------------------------
!  Denitrification is a microbially mediated process that converts
!  nitrate to molecular nitrogen (N2) under anaerobic conditions using DOC.
!  Assume denitrification happens when D.0. concentration is less than 
!  or equal to 2 mg/L.  A low D.O. concentration is considered anything 
!  less than 2 mg/L.  Denitrification reduces the amount of DOC present in 
!  the system.
!  Unit conversion: 2 mg/L is equal to 2 * 1.0E-03 kg/m^3
!------------------------------------------------------------------------------
 IF (DO_DO2) THEN

       IF (BIO_N(K,JDO2) <= (2.0 * 1.0E-03)) THEN
           DENITRIFICATION_C = ( BIO_N(K,JNO3) / (BIO_N(K,JNO3) + KHDENITR) )   &
                             * KDENITR * BIO_N(K,JDOC)   
       ELSE
           DENITRIFICATION_C = 0.0
       ENDIF

       !SUM_DENITR_C(k) = SUM_DENITR_C(k) + DENITRIFICATION_C

       DTM(K,JDOC) =  DTM(K,JDOC) - DENITRIFICATION_C 

 ENDIF

ENDDO
!
!------------------------------------------------------------------------------
!

RETURN

END SUBROUTINE CARBON

