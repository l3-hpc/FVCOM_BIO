SUBROUTINE GREENS(BIO_N,T_BIO,PAR)
!------------------------------------------------------------------------------
!-   Purpose and Methods : Greens kinetics calculations
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
USE FLAGS, ONLY: SILIM

IMPLICIT NONE

!
!------------------------------------------------------------------------------
!
REAL(SPP), INTENT(IN) :: PAR
REAL(SPP), INTENT(IN) :: BIO_N(KBV,NTT),T_BIO(KBV)

INTEGER :: K
REAL(SPP)  :: TEMP
REAL(SPP) :: FN               ! Nutrient limitation factor
REAL(SPP) :: PO4AVL           ! Available phosphate
REAL(SPP) :: NDIS             ! Dissolved nitrogen

DO K=1,KBV

!
!------------------------------------------------------------------------------
!
TEMP=T_BIO(K)

!------------------------------------------------------------------------------
!  Greens metabolism rate
!------------------------------------------------------------------------------
BMG(k) = BMRG * EXP(KTBG * (TEMP - TRG))


!------------------------------------------------------------------------------
!                             Light limitation
!------------------------------------------------------------------------------
! Reference: Jassby, Alan D. and Trevor Platt. July 1976. 
!            Mathematical formulation of the relationship between
!            photosynthesis and light for phytoplankton.
!            Limnology and Oceanography. Vol. 21(4). 
!
! ALPHA_GRE = slope of light saturation curve at low light levels, i.e. in 
!             the linear range.  Light saturation curve is also known as a
!             Photosynthesis-Light (P-I) curve.
!             Its units are g C [g Chl a]^-1 h^-1 (umol quanta)^-1 m^-2 s^-1.
!
! PBMAX_GRE = specific photosynthetic rate at optimal illumination, i.e the
!             assimilation number. It's the maximum assimilation number 
!             normalized to phytoplankton biomass; biomass is abbreviated as
!             B. Its units are g C [g Chl a]^-1 h^-1.
!
! The light limitation equation is parametrized in terms of PAR.
! PAR = photosynthetic active radiation = 400nm - 700nm.  The part
! of the electromagnetic spectrum that can be utilized for photosynthesis
! is between 400 and 700 nm and the energy in this part of the spectrum
! is called photosynthetic active radiation or PAR. Short-wave radiation 
! or solar radiation is the energy of wavelengths in the solar spectrum.
! The solar spectrum is between 300 and 3000 nm. Approximately 48% of the
! energy in the solar spectrum is PAR: see Brad Penta et al., 
! "An underwater light attenuation scheme for marine ecosystem models"  
! Optics Express, vol. 16, no. 21 (2008). 
! The units of PAR are umol/m2/s. 
!
! IAVG = average solar radiation in Watts/m^2. We need to convert Watts/m^2
! to umol/(m^2 s).  This is done by using the factor 4.57 which converts 
! W/m^2 to umol/(m^2 s).  
! See website "http://www.egc.com/useful_info_lighting.php" for more
! details regarding the conversion factor of 4.57.
!
!------------------------------------------------------------------------------
   IFG(K)  = TANH (ALPHA_GRE * PAR / PBMAX_GRE)

!------------------------------------------------------------------------------
! Nutrient limitations
!------------------------------------------------------------------------------

!-------------------------------------
! Dissolved phosphorous
!-------------------------------------
   PO4AVL = BIO_N(K,JSRP) + AVFRAC * BIO_N(K,JDOP)   

!-------------------------------------
! Dissolved nitrogen
!-------------------------------------
   NDIS = BIO_N(K,JNH4) + BIO_N(K,JNO3) + AVFRACDON *  BIO_N(K,JDON) 

!-------------------------------------
! Phosphorous limitation
!-------------------------------------
   PFG(K) = PO4AVL / (KHPG + PO4AVL)             

!-------------------------------------
! Nitrogen limitation
!-------------------------------------
   IF (NDIS > 1.0E-20) THEN
       NFG(K) = NDIS / (KHNG + NDIS)            
   ELSE
       NFG(K) = 0.0
   ENDIF

!------------------------------------------------------------------------------
!  Temperature and nutrient limitations
!------------------------------------------------------------------------------
   IF (TEMP < TMG) THEN
       TFG(K) = EXP(-KTGG1 * (TEMP - TMG)**2)                       
   ELSE
       TFG(K) = EXP(-KTGG2 * (TMG - TEMP)**2)
   ENDIF
   IF (SILIM == 3) THEN 
       FN = PFG(K)                      ! Phosphorus limitation only
   ELSE
       FN = MIN(NFG(K),PFG(K))       ! Nutrient limitation
   ENDIF
   PG(K) = PMG * FN * IFG(K) * TFG(K)      ! Production
   IF (PG(K) < 0.0) PG(K) = 0.0


!------------------------------------------------------------------------------
!  Update green algae time derivative
!------------------------------------------------------------------------------
   DTM(K,JGRE) = DTM(K,JGRE) + (PG(K) - BMG(K)) * BIO_N(K,JGRE) - PRG(K)                  

!   PG_AVG(k)  = PG_AVG(k)  + PG(k) *  f(JGRE) * Vol * real(dT,4)

!
!------------------------------------------------------------------------------
!
ENDDO

RETURN

END SUBROUTINE GREENS

