SUBROUTINE DIATOMS(BIO_N,T_BIO,PAR)
!------------------------------------------------------------------------------
!-   Purpose and Methods : Diatoms kinetics calculations
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
USE mod_1D, ONLY: SPP,KBV,KBVM1,NTT
USE STATES
USE EUT
USE FLAGS, ONLY: SILIM

IMPLICIT NONE

!
!------------------------------------------------------------------------------
!
REAL(SPP), INTENT(IN) :: PAR(KBV)
REAL(SPP), INTENT(IN) :: BIO_N(KBV,NTT),T_BIO(KBV)

INTEGER :: K
REAL(SPP)  :: TEMP
REAL(SPP) :: FN  ! Nutrient limitation factor
REAL(SPP) :: PO4AVL            ! Dissolved phosphorous
REAL(SPP) :: SA                ! Available silica
REAL(SPP) :: NDIS              ! Dissolved nitrogen

DO K=1,KBVM1

!
!------------------------------------------------------------------------------
!
TEMP = T_BIO(K)
!write(6,*) "TEMP",TEMP
!------------------------------------------------------------------------------
!  Diatom metabolism
!------------------------------------------------------------------------------
BMD(K) = BMRD * EXP(KTBD * (TEMP - TRD))                           
!write(6,*) "BMD,BMRD,KTBD,TRD",BMD(K),BMRD,KTBD,TRD
!------------------------------------------------------------------------------
!                             Light limitation
!------------------------------------------------------------------------------
! Reference: Jassby, Alan D. and Trevor Platt. July 1976. 
!            Mathematical formulation of the relationship between
!            photosynthesis and light for phytoplankton.
!            Limnology and Oceanography. Vol. 21(4). 
!
! ALPHA_DIA = slope of light saturation curve at low light levels, i.e. in 
!             the linear range. Light saturation curve is also known as a
!             Photosynthesis-Light (P-I) curve.
!             Its units are g C [g Chl a]^-1 h^-1 (umol quanta)^-1 m^-2 s^-1.
!
! PBMAX_DIA = specific photosynthetic rate at optimal illumination, i.e the
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

!   OUTPUT_PAR(k) = PAR  ! Store instantaneous PAR in this array.

   IFD(K) = TANH (ALPHA_DIA * PAR(K) / PBMAX_DIA)
   !write(6,*) "IFD,ALPHA_DIA,PAR,PBMAX_DIA",IFD(k),ALPHA_DIA,PAR(k),PBMAX_DIA
!------------------------------------------------------------------------------
!  Nutrient limitations
!------------------------------------------------------------------------------

!-------------------------------------
! Total phosphate
!-------------------------------------
   PO4AVL = BIO_N(K,JSRP) + AVFRAC * BIO_N(K,JDOP)   
!write(6,*) "PO4AVL,JSRP,AVFRAC,JDOP",PO4AVL,BIO_N(K,JSRP),AVFRAC,BIO_N(K,JDOP)
!-------------------------------------
! Available silica
!-------------------------------------
   SA = BIO_N(K,JSA)                              
!write(6,*) "SA",SA
!-------------------------------------
! Dissolved mitrogen
!-------------------------------------
   NDIS = BIO_N(K,JNH4) + BIO_N(K,JNO3) + AVFRACDON *  BIO_N(K,JDON)
!write(6,*) "NDIS,JNH4,JNO3,AVFRACDON,JDON",NDIS,BIO_N(K,JNH4),BIO_N(K,JNO3),AVFRACDON,BIO_N(K,JDON)
!-------------------------------------
! Phosphorous limitation
!-------------------------------------
   PFD(K) = PO4AVL / (KHPD + PO4AVL)            

!-------------------------------------
! Silica limitation
!-------------------------------------
   SFD(K) = SA / (KHSD + SA)              

!-------------------------------------
! Nitrogen limitation
!-------------------------------------
   IF (NDIS > 1.0E-20_SPP) THEN
       NFD(K) = NDIS / (KHND + NDIS)            
   ELSE
       NFD(K) = 0.0_SPP
   ENDIF


!------------------------------------------------------------------------------
!  Temperature limitation
!------------------------------------------------------------------------------
   IF (TEMP < TMD) THEN
      TFD(K) = EXP(-KTGD1 * (TEMP - TMD)**2)     
   ELSE
      TFD(K) = EXP(-KTGD2 * (TMD - TEMP)**2)
   ENDIF

!------------------------------------------------------------------------------
! Nutrient limitation 
!------------------------------------------------------------------------------
IF (SILIM == 1) FN = SFD(K) * MIN(NFD(K),PFD(K))   
IF (SILIM == 2) FN = MIN(NFD(K),SFD(K),PFD(K))   
IF (SILIM == 3) FN = SFD(K) * PFD(K)

!------------------------------------------------------------------------------
!  Production
!------------------------------------------------------------------------------
PD(K) = PMD * FN * IFD(K) * TFD(K)               

!------------------------------------------------------------------------------
!  Set production to zero if it's negative
!------------------------------------------------------------------------------
   IF (PD(K) < 0.0_SPP) PD(K) = 0.0_SPP


!------------------------------------------------------------------------------
!  Update diatom time derivatives 
!------------------------------------------------------------------------------
   DTM(K,JDIA) = DTM(K,JDIA) + (PD(K) - BMD(K)) * BIO_N(K,JDIA) - PRD(K)     


!    PD_AVG(k)  = PD_AVG(k)  + PD(k) * f(JDIA) * Vol * real(dT,4)

!
!------------------------------------------------------------------------------
!
ENDDO

!write(6,*) "DTM",DTM(:,JDIA)
!write(6,*) "PD",PD
!write(6,*) "BIO_N",BIO_N(:,JDIA)
!write(6,*) "PRD",PRD
!write(6,*) "BMD",BMD

RETURN

END SUBROUTINE DIATOMS
