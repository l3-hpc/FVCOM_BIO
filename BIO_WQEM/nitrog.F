SUBROUTINE NITROG(BIO_N,T_BIO)
!------------------------------------------------------------------------------
!-   Purpose and Methods : Kinetic sources and sinks of nitrogen are
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
!-
!- Revised: 09/02/14  W. Melendez   :   Added parallelization code.
!------------------------------------------------------------------------------
!
USE mod_1D, ONLY: SPP,KBV,KBVM1,NTT
USE STATES
USE EUT
USE FLAGS, ONLY: DO_DO2

IMPLICIT NONE

REAL(SPP), INTENT(IN) :: BIO_N(KBV,NTT),T_BIO(KBV)

INTEGER :: K
!
!------------------------------------------------------------------------------
!
REAL(SPP) :: TEMP
 
REAL(SPP) :: MNDON        ! Dissolved organic N mineralization mass derivative
REAL(SPP) :: HDRLPON      ! Labile detrital N hydrolysis mass mass derivative
REAL(SPP) :: HDRRPON      ! Refractory detrital N hydrolysis mass mass derivative
       
REAL(SPP) :: DON          ! Dissolved organic nitrogen
REAL(SPP) :: NDIS         ! Dissolved nitrogen

REAL(SPP) :: LPON         ! Contributions to labile detritus mass derivative
REAL(SPP) :: RPON         ! Contributions to refractory detritus mass derivative

REAL(SPP) :: NH4D         ! Diatom contribution to ammonia mass derivative
REAL(SPP) :: NH4G         ! Green algae contribution to ammonia mass derivative

REAL(SPP) :: NO3D         ! Diatom contribution to nitrate mass derivative
REAL(SPP) :: NO3G         ! Green algae contribution to nitrate mass derivative
       
REAL(SPP) :: ALGN         ! Algal nitrogen (??)

REAL(SPP) :: NT           ! Nitrification rate
REAL(SPP) :: FTN          ! Nitrification temperature dependence function

REAL(SPP) :: KDON         ! Dissolved organic nitrogen  mineralization rate [s^-1]
REAL(SPP) :: KLPON        ! Labile particulate nitrogen  hydrolysis rate [s^-1]
REAL(SPP) :: KRPON        ! Refractory particulate nitrogen  hydrolysis rate [s^-1]

REAL(SPP) :: ZDEATH              ! Zooplankton concentration from mortality
REAL(SPP) :: NH4ZOO              ! Ammonia released by zooplankton mortality
REAL(SPP) :: DONZOO              ! Dissolved organic N released by zoo. mortality
REAL(SPP) :: LPONZOO             ! Labile organic N released by zoo. mortality
REAL(SPP) :: RPONZOO             ! Refractory organic N released by zoo. mortality

REAL(SPP) :: ALG                 ! Total phytoplankton concentration (dia+gre)
REAL(SPP) :: TIN                 ! Total inorganic nitrogen (no3+nh4)

REAL(SPP) :: FNH4                ! NH4 fraction consumed by production
REAL(SPP) :: FNO3                ! NO3 fraction consumed by production
REAL(SPP) :: FDON                ! DON fraction consumed by production

!------------------------------------------------------------------------------
! Dissolved Oxygen limitation variable
!------------------------------------------------------------------------------
REAL(SPP) :: DO2_LIMITATION 

!------------------------------------------------------------------------------
! Denitrification rate in units of kg/(m^3 * s)
!------------------------------------------------------------------------------
REAL(SPP) :: DENITRIFICATION_N


DO K=1,KBVM1

!
!------------------------------------------------------------------------------
!
TEMP = T_BIO(K)

!------------------------------------------------------------------------------
!  Nitrification 
!------------------------------------------------------------------------------
  IF (TEMP < TMNT) THEN
      FTN = EXP(-KTNT1 * (TEMP - TMNT)**2)
  ELSE
      FTN = EXP(-KTNT2 * (TEMP - TMNT)**2)
  ENDIF

  IF (DO_DO2) THEN
      DO2_LIMITATION = BIO_N(K,JDO2) / ( KHDONT + BIO_N(K,JDO2) )
  ELSE
      DO2_LIMITATION = 1 
  ENDIF

   IF (KHNNT > 0.0) THEN
       NT = FTN * NTM * BIO_N(K,JNH4) / ( KHNNT + BIO_N(K,JNH4) ) &
          * DO2_LIMITATION 
   ELSEIF (BIO_N(K,JNH4) > 0.0) THEN
       NT = FTN * NTM * DO2_LIMITATION
   ELSE
       NT = 0.0
   ENDIF 

  DTM(K,JNH4) = DTM(K,JNH4) - NT 
  DTM(K,JNO3) = DTM(K,JNO3) + NT 
                                                                 
!------------------------------------------------------------------------------
! Algae sources/sinks
!------------------------------------------------------------------------------
!  ALGN --> Algal nitrogen  
!  FNO3 --> NO3 fraction
!  FNH4 --> NH4 fraction
!  NH4D --> Ammonia from diatoms
!  NH4G --> Ammonia from greens
!  NO3D --> Nitrate from diatoms
!  NO3G --> Nitrate from greens
!  DON  --> Dissolved organic nitrogen mass rate
!  LPON --> Labile particulate detrital nitrogen mass rate
!  RPON --> Refactory particulate detrital nitrogen mass rate
!  DTM(k,JNH4) --> Ammonia mass derivative
!  DTM(k,JNO3) --> Nitrate mass derivative
!  DTM(k,JDON) --> Dissolved organic N mass derivative
!  DTM(k,JLON) --> Labile detrital N mass derivative
!  DTM(k,JRON)-->  Refractory detrital N mass derivative
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Diatoms  
!------------------------------------------------------------------------------

   ALGN = ANCP * BIO_N(K,JDIA)                         
   NDIS = BIO_N(K,JNH4) + BIO_N(K,JNO3) + AVFRACDON *  BIO_N(K,JDON)           
   IF (ABS(NDIS) <= TINY(BIO_N(K,JNH4))) THEN
      FNO3 = 0.0
      FNH4 = 0.0
      FDON = 0.0
   ELSE
      FNO3 = BIO_N(K,JNO3) / NDIS        
      FNH4 = BIO_N(K,JNH4) / NDIS
      FDON = (AVFRACDON *  BIO_N(K,JDON)) / NDIS        
   ENDIF

   NH4D = (FNID * BMD(K) - FNH4 * PD(K)) * ALGN  &         
        +  FNIP * (1.0 - GREFF) * PRD(K) * ANCP          
   NO3D = -FNO3 * PD(K) * ALGN                                
   DON = (FNDD * BMD(K) - FDON * PD(K)) * ALGN +  &
         (FNDP * (1.0 - GREFF) * PRD(K) * ANCP)
   LPON = (FNLD * BMD(K) * ALGN) + (FNLP * (1.0 - GREFF) * PRD(K) * ANCP)
   RPON = (FNRD * BMD(K) * ALGN) + (FNRP * (1.0 - GREFF) * PRD(K) * ANCP)
   DTM(K,JNH4) = DTM(K,JNH4) + NH4D              
   DTM(K,JNO3) = DTM(K,JNO3) + NO3D              
   DTM(K,JDON) = DTM(K,JDON) + DON               
   DTM(K,JLON) = DTM(K,JLON) + LPON             
   DTM(K,JRON) = DTM(K,JRON) + RPON             

!------------------------------------------------------------------------------
!  Greens
!------------------------------------------------------------------------------
   ALGN = ANCP * BIO_N(K,JGRE)                                      
   NDIS = BIO_N(K,JNH4) + BIO_N(K,JNO3) + AVFRACDON *  BIO_N(K,JDON)
   IF (ABS(NDIS) <= TINY(BIO_N(K,JNH4))) THEN
      FNO3 = 0.0
      FNH4 = 0.0
      FDON = 0.0
   ELSE
      FNO3 = BIO_N(K,JNO3) / NDIS
      FNH4 = BIO_N(K,JNH4) / NDIS
      FDON = (AVFRACDON *  BIO_N(K,JDON)) / NDIS        
   ENDIF

   NH4G = (FNIG * BMG(K) - FNH4 * PG(K)) * ALGN        &           
        +  FNIP * (1.0 - GREFF) * PRG(K) * ANCP          
   NO3G = -FNO3 * PG(K) * ALGN                                     
   DON = (FNDG * BMG(K) - FDON * PG(K)) * ALGN + &
         FNDP * (1.0 - GREFF) * PRG(K) * ANCP   
   LPON = FNLG * BMG(K) * ALGN + FNLP * (1.0 - GREFF) * PRG(K) * ANCP  
   RPON = FNRG * BMG(K) * ALGN + FNRP * (1.0 - GREFF) * PRG(K) * ANCP  
   DTM(K,JNH4) = DTM(K,JNH4) + NH4G              
   DTM(K,JNO3) = DTM(K,JNO3) + NO3G             
   DTM(K,JDON) = DTM(K,JDON) + DON               
   DTM(K,JLON) = DTM(K,JLON) + LPON              
   DTM(K,JRON) = DTM(K,JRON) + RPON              

!------------------------------------------------------------------------------
!  Zooplankton
!------------------------------------------------------------------------------
! ZDEATH --> Zooplankton mortality rate: it was changed from a linear to a
!            quadratic dependence on the zooplankton concentration.  See
!            zoo.F90 subroutine for more details about this change.
!            [ZDTH] = 1/(kg/m^3 s)
! NH4ZOO --> NH4 component
! DONZOO --> DON component
! LPONZOO --> Labile particulate component
! RPONZOO --> Refractory particulate component
! DTM(k,JNH4) --> NH4 time derivative
! DTM(k,JDON) --> DON time derivative
! DTM(k,JLON) --> Labile part. time derivative
! DTM(k,JRON) --> Refractory part. time derivative
!------------------------------------------------------------------------------
   ZDEATH = ZDTH * BIO_N(K,JZOO) * BIO_N(K,JZOO)                                 
   NH4ZOO  = FNIZ * ANCP * ZDEATH                               
   DONZOO  = FNDZ * ANCP * ZDEATH                               
   LPONZOO = FNLZ * ANCP * ZDEATH                               
   RPONZOO = FNRZ * ANCP * ZDEATH                               
   DTM(K,JNH4) = DTM(K,JNH4) + NH4ZOO           
   DTM(K,JDON) = DTM(K,JDON) + DONZOO 
   DTM(K,JLON) = DTM(K,JLON) + LPONZOO
   DTM(K,JRON) = DTM(K,JRON) + RPONZOO
   

!------------------------------------------------------------------------------
!  Mineralization and hydrolysis
!------------------------------------------------------------------------------
! ALG --> Total phytoplankton concentration (dia+gre)
! TIN --> Total inorganic nitrogen
! KDON --> Maximum DON mineralization rate
! KLPON --> Maximum LON hydrolysis rate
! KRPON --> Maximum RON hydrolysis rate
! MNDON --> N mineralization mass derivative
! HDRLPON --> Labile detrital N hydrolysis mass derivative
! HDRRPON --> Refractory detrital N hydrolysis mass derivative
! DTM(k,JNH4) --> ammonia mass derivative
! DTM(k,JDON) --> dissolved organic N mass derivative
! DTM(k,JLON) --> particulate N mass derivative
! DTM(k,JRON) --> particulate N mass derivative
!------------------------------------------------------------------------------
   ALG = BIO_N(K,JGRE) + BIO_N(K,JDIA)                          
   TIN = BIO_N(K,JNO3) + BIO_N(K,JNH4)                          

   KDON =  KDN + ALG * KDNALG * KHN / (KHN + TIN)                     
   KLPON = KLN + ALG * KLNALG * KHN / (KHN + TIN)                     
   KRPON = KRN + ALG * KRNALG * KHN / (KHN + TIN)                     

   MNDON   =  KDON * FTMNL(K) * BIO_N(K,JDON)                   
   HDRLPON = KLPON * FTHDR(K) * BIO_N(K,JLON)                   
   HDRRPON = KRPON * FTHDR(K) * BIO_N(K,JRON)                   

   DTM(K,JNH4) = DTM(K,JNH4) + MNDON 

   DTM(K,JDON) = DTM(K,JDON) + HDRLPON + HDRRPON - MNDON 
   DTM(K,JLON) = DTM(K,JLON) - HDRLPON
   DTM(K,JRON) = DTM(K,JRON) - HDRRPON


!------------------------------------------------------------------------------
!  Denitrificaiton
!  Denitrification is a microbially mediated process that converts
!  nitrate to molecular nitrogen (N2) under anaerobic conditions using DOC.
!  Assume denitrification happens when D.0. concentration is less than 
!  or equal to 2 mg/L.  A low D.O. concentration is considered anything 
!  less than 2 mg/L.  Denitrification reduces the amount of nitrate present
!  in the system.
!  Unit conversion: 2 mg/L is equal to 2 * 1.0E-03 kg/m^3
!------------------------------------------------------------------------------
IF (DO_DO2) THEN

       IF (BIO_N(K,JDO2) <= (2 * 1.0E-03)) THEN
           DENITRIFICATION_N = BIO_N(K,JNO3) / ( BIO_N(K,JNO3) + KHDENITR ) &
                             * ANCP * KDENITR * BIO_N(K,JDOC)   
       ELSE
           DENITRIFICATION_N = 0
       ENDIF
!       SUM_DENITR(K) = SUM_DENITR(K) + DENITRIFICATION_N

       DTM(K,JNO3) = DTM(K,JNO3) - DENITRIFICATION_N 

ENDIF


ENDDO
!
!------------------------------------------------------------------------------
!

RETURN

END SUBROUTINE NITROG
