    module cgem_utils 
    use schism_glbl, only : rkind
    implicit none

    private :: rnitrate

    contains

!**************************************************************
  FUNCTION thermocons( S, T)     RESULT(RK)

  !------------------------------------------------------------
  ! Compute all the equilibrium constants (RK) for CO2 system.
  ! These are temperature and salinity dependent.
  ! from Whitfield and Turner, The Carbon Dioxide System
  ! in Estuarties--an Inorganic Perspective, The Science of the
  ! total Environment, 49(1986) 235-255, Elsevier.
  !------------------------------------------------------------
    IMPLICIT NONE

    REAL(RKIND),  INTENT(IN) :: S  ! Salinity
    REAL(RKIND),  INTENT(IN) :: T  ! Water Temperature(deg C)

!-----------------------------------------------------------

    REAL(RKIND), DIMENSION(8), PARAMETER :: A0 = (/&
    & 290.9097, 207.6548, 148.0248, 0.0221, &
    & 0.5709  , 0.9805  , 1.4853  , 0.5998 /)

    REAL(RKIND), DIMENSION(8), PARAMETER :: A1 = (/&
    & 14554.21, 11843.79, 8966.9 , 34.02,   &
    & -84.25  , -92.65  , -192.69, -75.25  /)

    REAL(RKIND), DIMENSION(8), PARAMETER :: A2 =   &
    &  (/ 45.0575, 33.6485, 24.4344, 0.0,   &
    &     0.0    , 0.0    , 0.0    , 0.0   /)

    REAL(RKIND),DIMENSION(8),PARAMETER ::  B0 =    &
    & (/ 0.0   , 0.0   , 0.0   , 0.0   ,    &
    &   -1.632 , -3.294, -5.058, -1.767    /)

    REAL(RKIND), DIMENSION(9) :: RK

    REAL(RKIND)    :: SP
    REAL(RKIND)    :: TP
    REAL(RKIND)    :: TP_OVER_100
    REAL(RKIND)    :: RK1S
    INTEGER :: I
    INTEGER :: J
!-----------------------------------------------------------
    SP          = MAX( S, 0.0 )
    TP          = T + 273.15     ! Absolute temperature (deg K).
    TP_OVER_100 = TP * 0.01

    ! Temperature dependence of the thermodynamic stability:
    DO I = 1, 3
      RK(I) = EXP( A0(I) - A1(I) / TP - A2(I) * LOG(TP) ) ! %Ko(i)
    END DO

    ! Sal. and temp. dependence of the stability constant K1:
    I = 1
    DO J = 4, 5
      RK( J ) = EXP( LOG( RK( I ) ) + (A0( J ) + A1( J ) / TP +        &
      &         A2( J ) * LOG( TP ) ) * SQRT( SP )                     &
      &       + B0( J ) * SP / 100.0 )
    END DO

    I = 2
    DO J = 6, 7
      RK( J ) = EXP( LOG( RK( I ) ) + ( A0( J ) + A1( J ) / TP +       &
      &         A2( J ) * LOG( TP ) ) * SQRT( SP )                     &
      &       + B0( J ) * SP / 100.0 )
    END DO

    ! Salinity and temperature dependence of K0:
    RK( 8 ) = EXP( -58.0931 + 90.5069 * ( 100.0 / TP ) +               &
    &           22.2940 * LOG( TP_OVER_100 ) +                         &
    &          (0.027766 - 0.025888 * TP_OVER_100 +                    &
    &           0.0050578 * TP_OVER_100 * TP_OVER_100 ) * SP)

    ! Addition of RK(9) from Sediment Diagenesis Model:
      RK1S = 2.527+1359.96/TP -0.206*SP**(1.0/3.0)
      RK1S = 10.0**(-RK1S)
      RK(9)= RK1S

    RETURN
  END FUNCTION thermocons


!----------------------------------------------------------------------------
  FUNCTION salclosed ( S, TC, T, PH )       RESULT(SPC)
!--------------------------------------------------------------------
! Closed CO2 System Calculations. SPC in units concentration mmol/m3.
!--------------------------------------------------------------------
    IMPLICIT NONE

    REAL(RKIND), INTENT(IN)    :: S     ! Salinity 
    REAL(RKIND), INTENT(IN)    :: TC    ! TC       - Dissolved Inorganic Carbon (DIC) (mmol m-3) in upper
                                 !            model layer
    REAL(RKIND), INTENT(IN)    :: T     ! Temperature (Celsius)
    REAL(RKIND), INTENT(IN)    :: PH    ! dimensionless
!-----------------------------------------------
    REAL(RKIND)                :: RK0        = 0.0
    REAL(RKIND)                :: RK1        = 0.0
    REAL(RKIND)                :: RK2        = 0.0
    REAL(RKIND)                :: RK1RK2     = 0.0
    REAL(RKIND)                :: H1         = 0.0
    REAL(RKIND)                :: H1_SQUARED = 0.0
    REAL(RKIND)                :: ALPHE0     = 0.0
    REAL(RKIND)                :: ALPHE1     = 0.0
    REAL(RKIND)                :: ALPHE2     = 0.0
    REAL(RKIND)                :: HCO3C      = 0.0
    REAL(RKIND)                :: CO3C       = 0.0
    REAL(RKIND)                :: H2CO3C     = 0.0

    REAL(RKIND),DIMENSION(9) :: RK           = 0.0
    REAL(RKIND),DIMENSION(6) :: SPC
!-----------------------------------------------

  !-----------------------------------------------------------
  ! Compute all the equilibrium constants (RK) for CO2 system.
  ! These are temperature and salinity dependent.
  ! from Whitfield and Turner (1986).
  !-----------------------------------------------------------
    RK     = thermocons( S, T)
    RK1    = RK( 4 )
    RK2    = RK( 6 )
    RK0    = RK( 8 )
    RK1RK2 = RK1 * RK2

    !-----------------------------------------------------------------
    ! pH  is given! The expressions below are from Cappellen and Wang,
    !               March 1996, American Journal of Science, "Cycling
    !               of iron and manganese in sfc sediments: A Genl.
    !               Theory for the Coupled Transport and Reaction of
    !               CO2, Oxygen, Nitrogen, Sulfur, Iron, and Manganese
    !-----------------------------------------------------------------
    H1         = 10.0 ** ( -PH )
    H1_SQUARED = H1 * H1
    ALPHE0     = 1.0/( 1.0               + RK1/H1 + RK1RK2/H1_SQUARED )
    ALPHE1     = 1.0/( H1/RK1            + 1.0    + RK2/H1 )
    ALPHE2     = 1.0/( H1_SQUARED/RK1RK2 + H1/RK2 + 1.0 )

    HCO3C      = TC * ALPHE1
    CO3C       = TC * ALPHE2
    H2CO3C     = TC * ALPHE0

! Each element of SPC has concentration units of mmol/m3
    SPC( 1 )   = H2CO3C
    SPC( 2 )   = HCO3C
    SPC( 3 )   = CO3C
    SPC( 4 )   = RK1
    SPC( 5 )   = RK2
    SPC( 6 )   = RK0

    RETURN
  END FUNCTION salclosed
!---------------------------------------------------------------------

      FUNCTION sw_smow( T ) RESULT( RES )

  !------------------------------------------------------------------
  ! Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.
  ! Based on MATLAB code and the following references:
  ! SW_SMOW  $Revision: 1.3 $  $Date: 1994/10/10 05:51:46 $
  ! Copyright (C) CSIRO, Phil Morgan 1992.
  !
  ! INPUT:
  !   T = temperature [degree C (IPTS-68)]
  !
  ! OUTPUT:
  !   dens = density  [kg/m^3]
  !
  ! AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
  !
  ! DISCLAIMER:
  !   This software is provided "as is" without warranty of any kind.
  !   See the file sw_copy.m for conditions of use and licence.
  !
  ! REFERENCES:
  !     Unesco 1983. Algorithms for computation of fundamental properties of
  !     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !     UNESCO 1983 p17  Eqn(14)
  !
  !     Millero, F.J & Poisson, A.
  !     INternational one-atmosphere equation of state for seawater.
  !     Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
  !-------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(RKIND),INTENT(IN) :: T    ! Temperature (deg C)
    REAL(RKIND)            :: RES

    ! Locals:
    REAL(RKIND), PARAMETER :: A0 = 999.842594
    REAL(RKIND), PARAMETER :: A1 =   6.793952E-2
    REAL(RKIND), PARAMETER :: A2 =  -9.095290E-3
    REAL(RKIND), PARAMETER :: A3 =   1.001685E-4
    REAL(RKIND), PARAMETER :: A4 =  -1.120083E-6
    REAL(RKIND), PARAMETER :: A5 =   6.536332E-9

    RES = A0 + ( A1 + ( A2 + ( A3 + ( A4 + A5 * T ) * T ) * T ) * T ) * T

    RETURN
  END FUNCTION sw_smow

      FUNCTION sw_dens0( S, T ) RESULT( RES )

  !--------------------------------------------------------------------------
  ! Density of Sea Water at atmospheric pressure using UNESCO 1983 (EOS 1980).
  ! Based on MATLAB code and the following references:
  ! SW_DENS0  $Revision: 1.3 $  $Date: 1994/10/10 04:54:09 $
  !           Copyright (C) CSIRO, Phil Morgan 1992
  !
  ! INPUT:  (all must have same dimensions)
  !   S = salinity    [psu      (PSS-78)]
  !   T = temperature [degree C (IPTS-68)]
  !
  ! OUTPUT:
  !   dens0 = density  [kg/m^3] of salt water with properties S,T,
  !           P=0 (0 db gauge pressure)
  !
  ! AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
  !
  ! DISCLAIMER:
  !   This software is provided "as is" without warranty of any kind.
  !   See the file sw_copy.m for conditions of use and licence.
  !
  ! REFERENCES:
  !     Unesco 1983. Algorithms for computation of fundamental properties of
  !     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !
  !     Millero, F.J. and  Poisson, A.
  !     International one-atmosphere equation of state of seawater.
  !     Deep-Sea Res. 1981. Vol28A(6) pp625-629.
  !------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(RKIND), INTENT(IN) :: S    ! Salinity (psu)
    REAL(RKIND), INTENT(IN) :: T    ! Temperature (deg C)
    REAL(RKIND) :: RES

    ! Locals from UNESCO 1983 eqn(13) p17.
    REAL(RKIND), PARAMETER:: B0 =  8.24493E-1
    REAL(RKIND), PARAMETER:: B1 = -4.0899E-3
    REAL(RKIND), PARAMETER:: B2 =  7.6438E-5
    REAL(RKIND), PARAMETER:: B3 = -8.2467E-7
    REAL(RKIND), PARAMETER:: B4 =  5.3875E-9

    REAL(RKIND), PARAMETER:: C0 = -5.72466E-3
    REAL(RKIND), PARAMETER:: C1 = +1.0227E-4
    REAL(RKIND), PARAMETER:: C2 = -1.6546E-6

    REAL(RKIND), PARAMETER:: D0 =  4.8314E-4

    RES = sw_smow( T ) +                                                &
    &     ( B0 + ( B1 + ( B2 + ( B3 + B4 * T ) * T ) * T ) * T ) * S +  &
    &     ( C0 + ( C1 + C2 * T ) * T ) * S * SQRT( S ) + D0 * S * S

    RETURN
  END FUNCTION sw_dens0


!------------------------------------------------------------------
  FUNCTION gas_exchange( T, S, TC, H, PH, PCO2 ) RESULT( RES )
!------------------------------------------------------------------
! Gas exchange model from Eldridge and Cifuentes (2000).
! Provides and estimate of the loss or gain of CO2 from the atmosphere
! The atmosphere pCo2 concentration was from Whitfield and Turner (1986)
! and may have changes by now.  D/z is in cm s-1 and depth is in meters.
! Multiplication by 1000 and division by 100 provides units of mmol C m-2
! s-1.  Final units are mmol C m-2 d-1.
!------------------------------------------------------------------
! T        - Temperature (deg C) of upper model layer.
! S        - Salinity (psu) of upper model layer.
! TC       - Dissolved Inorganic Carbon (DIC) (mmol m-3) in upper
!            model layer
! H        -- Thickness of the upper model layer (m.) contacting the
!             atmosphere
! PH       -- pH concentration of seawater.
! PCO2     -- CO2 concentration in atmosphere in ppmv.
!             Table 1 in Huang et al. 2015 for LA shelf, average=380
!------------------------------------------------------------------

    IMPLICIT NONE

    REAL(RKIND)   ,INTENT(IN)  :: T
    REAL(RKIND)   ,INTENT(IN)  :: S
    REAL(RKIND)   ,INTENT(IN)  :: TC
    REAL(RKIND)   ,INTENT(IN)  :: H
    REAL(RKIND)   ,INTENT(IN)  :: PCO2
    REAL(RKIND)   ,INTENT(IN)  :: PH
!--------------------------------------------
    REAL(RKIND), PARAMETER      :: DZ      = 0.005
                                        ! combines diffusion through
                                        ! water-atm layer and
                                        ! thickness(cm s-1).
    REAL(RKIND)                 :: ATM_CO2
    REAL(RKIND)                 :: H2CO3C
    REAL(RKIND)                 :: PW
    REAL(RKIND)                 :: RES
    REAL(RKIND)                 :: RK0
    REAL(RKIND)   ,DIMENSION(6) :: SPC

    ! Use an updated salclosed function to get the distribution of
    ! carbonate alkalinity species and H+.  Whitman and Turner (1986).
    SPC     = salclosed( S, TC, T, PH )  ! SPC units Concentration mmol/m3

    H2CO3C  = SPC( 1 )   ! Concentration of H2CO3 in units of mmol/m3
    RK0     = SPC( 6 )   ! Henry's constant (mol kg-1 atm-1)
    PW      = sw_dens0( S, T )              ! water density  [kg/m3]

    ATM_CO2  = PCO2 * RK0 * (PW * 1.0E-3)

    RES = -DZ  * ( ATM_CO2 - H2CO3C ) / ( H * 100.0) !mmol C m-2 s-1

    RETURN
  END FUNCTION gas_exchange

!***********************************************************************
      FUNCTION SchmidtNumber( S, T, which ) RESULT( RES )

  !-------------------------------------------------------------------------
  ! INPUT:  (all must have same dimensions)
  !   S = salinity    [psu      (PSS-78)]
  !   T = temperature [degree C (IPTS-68)]
  !
  ! OUTPUT:
  !   Schmidt Number for Oxygen (which = 0)
  !or Schmidt Number for CO2    (which = 1)
  !
  ! REFERENCES:
  !     Unesco 1983. Algorithms for computation of fundamental properties of
  !     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(RKIND), INTENT(IN) :: S    ! Salinity (psu)
    REAL(RKIND), INTENT(IN) :: T    ! Temperature (deg C)
    REAL(RKIND) :: RES
    REAL(RKIND) :: muSW,KT
    REAL(RKIND) :: D_02, D_C02
    INTEGER, INTENT(IN) :: which !0 for Oxygen
                                 !1 for CO2

    REAL(RKIND), PARAMETER:: K0 = 1.7910
    REAL(RKIND), PARAMETER:: K1 = -6.144e-2
    REAL(RKIND), PARAMETER:: K2 = 1.4510e-3
    REAL(RKIND), PARAMETER:: K3 = -1.6826e-5
    REAL(RKIND), PARAMETER:: K4 = -1.5290e-4
    REAL(RKIND), PARAMETER:: K5 = 8.3885e-8
    REAL(RKIND), PARAMETER:: K6 = 2.4727e-3
    REAL(RKIND), PARAMETER:: K7 = 6.0574e-6
    REAL(RKIND), PARAMETER:: K8 = -2.6760e-9
    REAL(RKIND), PARAMETER:: K9 = 4.8429e-5
    REAL(RKIND), PARAMETER:: K10 = -4.7172e-6
    REAL(RKIND), PARAMETER:: K11 = 7.5986e-8

    REAL(RKIND), PARAMETER:: D1 = 0.2604e-5
    REAL(RKIND), PARAMETER:: D2 = 0.006383e-5
    REAL(RKIND), PARAMETER:: D3 = 0.1954e-5
    REAL(RKIND), PARAMETER:: D4 = 0.005089e-5

    REAL(RKIND), PARAMETER:: P = 1.01325

    KT = 273.15 + T

    muSW = K0 + K1*T + K2*T*T + K3*T*T*T + K4*P + K5*P*P + K6*S +  & 
  & T*(K7*P + K8*P*P) + S*(K9*T + K10*T*T + K11*T*T*T)

    if(which.eq.0) then
     D_02 = D1 + D2*KT/muSW
     RES = muSW/(sw_dens0(S,T)*D_02)
    endif

    if(which.eq.1) then
     D_C02 = D3 + D4*KT/muSW
     RES = muSW/(sw_dens0(S,T)*D_C02)
    endif


    RETURN
  END FUNCTION SchmidtNumber

!---------------------------------------------------------------------
   FUNCTION   o2sat(S, T)               RESULT(OSAT)

!---------------------------------------------------------------------
! Oxygen concentration at saturation at one atmosphere (umol/kg).
! Source: "The Solubility Of Nitrogen, Oxygen And Argon In Water And
!         Seawater" - Weiss (1970) _Deep Sea Research_ V17(4): 721-735.
! Based on MATLAB code by Edward T Peltzer, MBARI, revised: 1998-01-2
!---------------------------------------------------------------------

    IMPLICIT NONE

    REAL(RKIND)        , INTENT(IN) :: S  ! Salinity
    REAL(RKIND)        , INTENT(IN) :: T  ! Temperature (degrees C)

    ! Declare Locals:
!----------------------------------------------------------------------
    REAL(RKIND), PARAMETER:: CONSTANT1  =  273.15
    REAL(RKIND), PARAMETER:: CONSTANT2  =    0.01
    REAL(RKIND), PARAMETER:: CONSTANT3  = -177.7888
    REAL(RKIND), PARAMETER:: CONSTANT4  =  255.5907
    REAL(RKIND), PARAMETER:: CONSTANT5  =  146.4813
    REAL(RKIND), PARAMETER:: CONSTANT6  =   22.204
    REAL(RKIND), PARAMETER:: CONSTANT7  =   -0.037362
    REAL(RKIND), PARAMETER:: CONSTANT8  =    0.016504
    REAL(RKIND), PARAMETER:: CONSTANT9  =    0.0020564
    REAL(RKIND), PARAMETER:: ML_PER_UM  =   44.658
    REAL(RKIND)          :: T1
    REAL(RKIND)          :: OSAT ! Oxygen (umol/kg) --function RESULT
!----------------------------------------------------------------------

    T1 = ( T + CONSTANT1 ) * CONSTANT2

    OSAT = CONSTANT3 + CONSTANT4 /                                     &
    &      T1 + CONSTANT5 * LOG( T1 ) - CONSTANT6 * T1

    OSAT = OSAT + S * (CONSTANT7 + T1 * (CONSTANT8 - CONSTANT9 * T1))
    OSAT = EXP( OSAT )
    OSAT = OSAT * ML_PER_UM ! Oxygen (umol/kg)

    RETURN
  END FUNCTION o2sat

  
!****************************************************************      
!---------------------------------------------------------------------
   FUNCTION   Q10_T(T,K)               RESULT(Tadj)
!--------------------------------------------------------------------- 
! Q10 rate adjustment, taken from reaction subroutine
!---------------------------------------------------------------------    
   
    IMPLICIT NONE

!----------------------------------------------------------------------
    REAL(RKIND), INTENT(IN) :: T  ! Temperature (degrees C) 
    REAL(RKIND), INTENT(IN) :: K  ! Constant     
!----------------------------------------------------------------------    
    REAL(RKIND), PARAMETER  :: Tk = 25.
    REAL(RKIND)             :: FACTOR    
    REAL(RKIND)             :: Tadj ! Temperature adjusted rate function 
!----------------------------------------------------------------------       
    FACTOR = LOG10( 2.0 ) * 0.1 * ( Tk - T )
    Tadj = LOG10( K ) - FACTOR
    Tadj = 10.0 ** Tadj 

    RETURN

  END FUNCTION Q10_T

  
!--------------------------------------------------------------------
  SUBROUTINE Nitrification( O2, NH4, KO2, KNH4, nitmax, TEMP, R_11 )     
  !-------------------------------------------------------------------  
  ! Provides Nitrification rate for the model.
  ! Note that the time units are in days, not seconds. 
  !-------------------------------------------------------------------   
    IMPLICIT NONE
    
    REAL(RKIND), INTENT(IN) :: O2
    REAL(RKIND), INTENT(IN) :: NH4
    REAL(RKIND), INTENT(IN) :: KO2
    REAL(RKIND), INTENT(IN) :: KNH4
    REAL(RKIND), INTENT(IN) :: nitmax 
    REAL(RKIND), INTENT(IN) :: TEMP
    REAL(RKIND), INTENT(OUT):: R_11  

    REAL(RKIND) :: FACTOR
    REAL(RKIND) :: TQ1
    REAL(RKIND) :: TQ2
    REAL(RKIND) :: RQ1
    REAL(RKIND) :: RQ11
!-------------------------------------------------
    ! Use the Q10 relationship to determine the rates.

    TQ1  = 25.0
    TQ2  = TEMP
    FACTOR = LOG10( 2.0 ) * 0.1 * ( TQ1 - TQ2 )

    RQ1  = nitmax
    RQ11 = LOG10( RQ1 ) - FACTOR
    RQ11 = 10.0 ** RQ11

    R_11 = RQ11 * O2/(KO2+O2) * NH4/(KNH4+NH4)

    RETURN
  END SUBROUTINE Nitrification 

!***********************************************************************
  SUBROUTINE reaction ( OM1, OM2, O2, NO3, KG1, KG2, KO2, KstarO2, &
  &   KNO3, X1, Y1, Z1, X2, Y2, Z2, TEMP, RC )

  !----------------------------------------------------------
  ! Provides all reaction rates for the model.
  ! Note that the time units are in years, not days.
  !---------------------------------------------------------- 
    IMPLICIT NONE
   
    REAL(RKIND), INTENT(IN) :: OM1
    REAL(RKIND), INTENT(IN) :: OM2
    REAL(RKIND), INTENT(IN) :: O2
    REAL(RKIND), INTENT(IN) :: NO3
    REAL(RKIND), INTENT(IN) :: KG1
    REAL(RKIND), INTENT(IN) :: KG2
    REAL(RKIND), INTENT(IN) :: KO2
    REAL(RKIND), INTENT(IN) :: KstarO2
    REAL(RKIND), INTENT(IN) :: KNO3
    REAL(RKIND), INTENT(IN) :: X1
    REAL(RKIND), INTENT(IN) :: Y1
    REAL(RKIND), INTENT(IN) :: Z1
    REAL(RKIND), INTENT(IN) :: X2
    REAL(RKIND), INTENT(IN) :: Y2
    REAL(RKIND), INTENT(IN) :: Z2
    REAL(RKIND), INTENT(IN) :: TEMP
    REAL(RKIND), DIMENSION(10),INTENT(OUT):: RC
    
    ! Begin Locals:
!---------------------------------------------    
    REAL(RKIND):: FACTOR         
    REAL(RKIND):: RQ1           
    REAL(RKIND):: RQ2          
    REAL(RKIND):: TQ1         
    REAL(RKIND):: TQ2        
    REAL(RKIND):: RQ21      
    REAL(RKIND):: RQ22     
    REAL(RKIND):: RCT1           
    REAL(RKIND):: RCT2          
    REAL(RKIND):: R11          
    REAL(RKIND):: R12         
    REAL(RKIND):: R21        
    REAL(RKIND):: R22       
    REAL(RKIND):: R1_SUM    
    REAL(RKIND):: R2_SUM        
    REAL(RKIND):: FBNO3         
    REAL(RKIND):: GAM14      
    REAL(RKIND):: GAM24
    REAL(RKIND):: BET14
    REAL(RKIND):: BET24
    REAL(RKIND):: A14
    REAL(RKIND):: A24     
    REAL(RKIND):: RCH2O1         
    REAL(RKIND):: RCH2O2        
    REAL(RKIND):: RO2          
    REAL(RKIND):: RNO3        
    REAL(RKIND):: RPO4       
    REAL(RKIND):: RTC       
    REAL(RKIND):: RNH4
    REAL(RKIND):: RSi
    REAL(RKIND):: RALK
    REAL(RKIND):: RN2
    REAL(RKIND),DIMENSION(2) :: R   
!---------------------------------------------     
    ! End Locals   

    ! Use the Q10 relationship to determine the rates.
    ! Assume that TEMP is the maximum temperature.
    RQ1  = KG1
    RQ2  = KG2
    TQ1  = 25.0
    TQ2  = TEMP
    FACTOR = LOG10( 2.0 ) * 0.1 * ( TQ1 - TQ2 )
    RQ21 = LOG10( RQ1 ) - FACTOR
    RQ21 = 10.0 ** RQ21
    RQ22 = LOG10( RQ2 ) - FACTOR
    RQ22 = 10.0 ** RQ22

    ! Lets oxidants determine the rate of organic degradation using
    ! the full Monod relationship.
    ! Calculate the concentration of OMs from Flux.
    ! Calculate the fluxes using the concentration of OMs and 
    !    reaction rates.

    ! Provide Monod Ks, i.e. KNO3,  for denitrification--see
    ! parameter statement
    RCT1     = RQ21 * OM1
    RCT2     = RQ22 * OM2
    R( 1 ) = O2 / ( KO2 + O2 )
    R( 1 ) = MAX( 0.0, R( 1 ) )

    FBNO3 = rnitrate( O2, KstarO2 ) ! Feedback on denitrification by O2.    

    R( 2 ) = NO3 / ( KNO3 + NO3 ) * FBNO3
    R( 2 ) = MAX( 0.0, R( 2 ) )

    R11 = RCT1 * R( 1 )
    R12 = RCT1 * R( 2 )
    R1_SUM = R11 + R12
    
    R21 = RCT2 * R( 1 )
    R22 = RCT2 * R( 2 )
    R2_SUM = R21 + R22


    ! Reaction rates from TABLE 5 of Van Cappellen and Wang (1996).
    ! Use J. Lehrter's updated equations, where OM is first converted
    ! to NH3 (instead of converting directly to NO3)
    ! Use the stoichiometry from OM to determine O2 and NO3 use.
    GAM14 = ( 4.0 * X1 + 3.0 * Y1 ) / 5.0 / X1 ! denitrification NO3-/CH20
    GAM24 = ( 4.0 * X2 + 3.0 * Y2 ) / 5.0 / X2 ! denitrification NO3-/CH20

    BET14 = ( 2.0 * X1 + 4.0 * Y1 ) / 5.0 / X1 ! denitrification NO3-/CH20
    BET24 = ( 2.0 * X2 + 4.0 * Y2 ) / 5.0 / X2 ! denitrification NO3-/CH20
    RN2 = ( BET14 * R12 + BET24 * R22 )

    A14 = ( 4.0 * X1 + 3.0 * Y1 - 10.0 * Z1 ) / 5.0 / X1
    A24 = ( 4.0 * X2 + 3.0 * Y2 - 10.0 * Z2 ) / 5.0 / X2

    RO2 = -( R11 + R21 )

    RNO3 = - ( GAM14 * R12 + GAM24 * R22 )
 
    RALK = (A14 * R12 + A24 * R22)

    RNH4 = ( Y1 / X1 * R11 + Y2 / X2 * R21 )

    RSi    =  R1_SUM * Y1/X1 + R2_SUM * Y2/X2 

    RPO4   =  R1_SUM * Z1/X1 + R2_SUM * Z2/X2
    RCH2O1 = -R1_SUM   
    RCH2O2 = -R2_SUM   
    RTC    =  R1_SUM + R2_SUM   

    RC( 1 ) = RCH2O1
    RC( 2 ) = RCH2O2
    RC( 3 ) = RO2
    RC( 4 ) = RNO3
    RC( 5 ) = RPO4
    RC( 6 ) = RTC
    RC( 7 ) = RNH4
    RC( 8 ) = RSi 
    RC( 9 ) = RALK
    RC( 10 ) = RN2

    RETURN

  END SUBROUTINE reaction 
!***********************************************************************

!***********************************************************************
  FUNCTION rnitrate( O2, KO2 ) RESULT( RES ) 

  !--------------------------------------------------------------------  
  ! Monod type function and feedbacks for NO3 - as an electron acceptor.
  !--------------------------------------------------------------------   
    IMPLICIT NONE
    
    REAL(RKIND),INTENT(IN) :: O2  ! Oxygen concentration
    REAL(RKIND),INTENT(IN) :: KO2 ! Inhibition constant
    REAL(RKIND)            :: RES
    REAL(RKIND):: PO2

    PO2  = MAX( 0.0, O2 )

    RES  = KO2 / ( KO2 + PO2 )  

    RETURN
  END FUNCTION rnitrate  
!***********************************************************************

  end module cgem_utils

