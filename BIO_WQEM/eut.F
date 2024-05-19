MODULE EUT
!------------------------------------------------------------------------------
!-
!-  $Id: eut.F90,v 1.0.6.1 2014/08/26 22:54:04 wmelende Exp $
!-  $Locker:  $
!-
!-  Purpose and Methods : Include file for eutrophication variables.
!-
!-  Created: 12/09/97  M. Settles
!------------------------------------------------------------------------------
USE MOD_1D, ONLY: SPP,KBV,NTT

IMPLICIT NONE

SAVE

!DTM, need to allocate
REAL(SPP),ALLOCATABLE :: DTM(:,:)          ! Base metabolic rate for diatoms

!------------------------------------------------------------------------------
!  Miscellaneous arrays/variables; not parameters
!------------------------------------------------------------------------------

!calculated in diatoms, used in carbon, dissolved oxygen, nitrog, phosph, silica 
REAL(SPP),ALLOCATABLE :: BMD(:)          ! Base metabolic rate for diatoms

!calculated in greens, used in carbon, dissolved oxygen, nitrog, phosph, silica
REAL(SPP),ALLOCATABLE :: BMG(:)          ! Base metabolic rate for greens

!Calculated in carbon, used in carbon, nitrog, and phosph
REAL(SPP),ALLOCATABLE :: FTHDR(:)        ! Temperature correction for hydrolysis

!Calculated in carbon, used in carbon, dissolved_oxygen, nitrog, phosph
REAL(SPP),ALLOCATABLE :: FTMNL(:)        ! Temperature correction for mineralization

!Used for output, calculated in diatoms
REAL(SPP),ALLOCATABLE :: OUTPUT_PAR(:)    ! Daily cumulative PAR

!calculated in diatoms, used in dissolved oxygen, nitrog, phosph, silica
REAL(SPP),ALLOCATABLE :: PD(:)           ! Production rates for diatoms

!Used for output, calculated in diatoms
REAL(SPP),ALLOCATABLE :: PD_AVG(:)       ! Average production rates for 
                                              ! diatoms

!calculated in diatoms, used in dissolved oxygen, nitrog, phosph
REAL(SPP),ALLOCATABLE :: PG(:)           ! Production rates for greens

!Used for output, calculated in greens
REAL(SPP),ALLOCATABLE :: PG_AVG(:)       ! Average production rates for 
                                              ! greens

!Calculated in zoo, used in carbon, dissolved oxygen, nitrog, phosph, silica
REAL(SPP),ALLOCATABLE :: PRD(:)          ! Predation time derivative for diatoms

!Calculated in zoo, used in carbon, dissolved oxygen, nitrog, phosph
!FV REAL(SPP),ALLOCATABLE :: PRG(:)          ! Predation time derivative for greens
REAL(SPP),ALLOCATABLE :: PRG(:)          ! Predation time derivative for greens


!------------------------------------------------------------------------------
!  Fine segmentation growth limitation factors
!------------------------------------------------------------------------------
!Calculated and used in diatoms, used in output
REAL(SPP),ALLOCATABLE :: PFD(:)       ! Phosphorous limitation for diatoms

!Calculated and used in diatoms, used in output
REAL(SPP),ALLOCATABLE :: SFD(:)       ! Silica limitation for diatoms

!Calculated and used in diatoms, used in output
REAL(SPP),ALLOCATABLE :: NFD(:)       ! Nitrogen limitation for diatoms

!Calculated and used in diatoms, used in output
REAL(SPP),ALLOCATABLE :: IFD(:)       ! Light limitation for diatoms

!Calculated and used in diatoms, used in output
REAL(SPP),ALLOCATABLE :: TFD(:)       ! Temperature limitation for diatoms

!Calculated and used in greens, used in output
REAL(SPP),ALLOCATABLE :: PFG(:)       ! Phosphorous limitation for greens

!Calculated and used in greens, used in output
REAL(SPP),ALLOCATABLE :: NFG(:)       ! Nitrogen limitation for greens

!Calculated and used in greens, used in output
REAL(SPP),ALLOCATABLE :: IFG(:)       ! Light limitation for greens

!Calculated and used in greens, used in output
REAL(SPP),ALLOCATABLE :: TFG(:)       ! Temperature limitation for greens

!New for Droop
!diatoms
!Calculated in diatoms, used in nitrog
REAL(SPP),ALLOCATABLE :: QND(:)       !
!Calculated in diatoms, used in phosph
REAL(SPP),ALLOCATABLE :: QPD(:)       !
REAL(SPP) :: QminND = 0.003_SPP
REAL(SPP) :: QminPD = 0.003_SPP
!greens
!Calculated in greens, used in nitrog
REAL(SPP),ALLOCATABLE :: QNG(:)       !
!Calculated in greens, used in phosph
REAL(SPP),ALLOCATABLE :: QPG(:)       !
REAL(SPP) :: QminNG = 0.003_SPP
REAL(SPP) :: QminPG = 0.003_SPP
!nitrog
REAL(SPP) :: QmaxND = 1.0_SPP
REAL(SPP) :: QmaxNG = 1.0_SPP
REAL(SPP) :: QmaxPD = 1.0_SPP
REAL(SPP) :: QmaxPG = 1.0_SPP
REAL(SPP) :: KhIntND = 2.0E-05
REAL(SPP) :: KhIntNG = 2.0E-05
REAL(SPP) :: KhIntPD = 2.0E-05
REAL(SPP) :: KhIntPG = 2.0E-05
!eut.F90
REAL(SPP) :: FIntNID = 0.0_SPP    ! Fraction of inorganic nitrogen produced by algal metabolism
REAL(SPP) :: FIntNIG = 0.0_SPP    ! Fraction of inorganic nitrogen produced by algal metabolism
REAL(SPP) :: FIntNDD = 0.0_SPP    ! Fraction of dissolved nitrogen produced by algal metabolism
REAL(SPP) :: FIntNDG = 0.0_SPP    ! Fraction of dissolved nitrogen produced by algal metabolism
REAL(SPP) :: FIntNLD = 0.0_SPP    ! Fraction of labile nitrogen produced by algal metabolism
REAL(SPP) :: FIntNLG = 0.0_SPP    ! Fraction of labile nitrogen produced by algal metabolism
REAL(SPP) :: FIntNRD = 0.0_SPP    ! Fraction of refractory nitrogen produced by algal metabolism
REAL(SPP) :: FIntNRG = 0.0_SPP    ! Fraction of refractory nitrogen produced by algal metabolism
REAL(SPP) :: FIntPID = 0.0_SPP    ! Fraction of inorganic phosphorus produced by algal metabolism
REAL(SPP) :: FIntPIG = 0.0_SPP    ! Fraction of inorganic phosphorus produced by algal metabolism
REAL(SPP) :: FIntPDD = 0.0_SPP    ! Fraction of dissolved phosphorus produced by algal metabolism
REAL(SPP) :: FIntPDG = 0.0_SPP    ! Fraction of dissolved phosphorus produced by algal metabolism
REAL(SPP) :: FIntPLD = 0.0_SPP    ! Fraction of labile phosphorus produced by algal metabolism
REAL(SPP) :: FIntPLG = 0.0_SPP    ! Fraction of labile phosphorus produced by algal metabolism
REAL(SPP) :: FIntPRD = 0.0_SPP    ! Fraction of refractory phosphorus produced by algal metabolism
REAL(SPP) :: FIntPRG = 0.0_SPP    ! Fraction of refractory phosphorus produced by algal metabolism
REAL(SPP) :: FIntLuxNIP = 0.0_SPP ! Fraction of luxury inorganic nitrogen produced by predation
REAL(SPP) :: FIntLuxNDP = 0.0_SPP ! Fraction of luxury dissolved organic nitrogen produced by predation
REAL(SPP) :: FIntLuxNLP = 0.0_SPP ! Fraction of luxury labile nitrogen produced by predation
REAL(SPP) :: FIntLuxNRP = 0.0_SPP ! Fraction of luxury refractory nitrogen produced by predation
REAL(SPP) :: FIntStrNIP = 0.0_SPP ! Fraction of structural inorganic nitrogen produced by predation
REAL(SPP) :: FIntStrNDP = 0.0_SPP ! Fraction of structural dissolved organic nitrogen produced by predation
REAL(SPP) :: FIntStrNLP = 0.0_SPP ! Fraction of structural labile nitrogen produced by predation
REAL(SPP) :: FIntStrNRP = 0.0_SPP ! Fraction of structural refractory nitrogen produced by predation
REAL(SPP) :: FIntLuxPIP = 0.0_SPP ! Fraction of luxury inorganic phosphorus produced by predation
REAL(SPP) :: FIntLuxPDP = 0.0_SPP ! Fraction of luxury dissolved organic phosphorus produced by predation
REAL(SPP) :: FIntLuxPLP = 0.0_SPP ! Fraction of luxury labile phosphorus produced by predation
REAL(SPP) :: FIntLuxPRP = 0.0_SPP ! Fraction of luxury refractory phosphorus produced by predation
REAL(SPP) :: FIntStrPIP = 0.0_SPP ! Fraction of structural inorganic phosphorus produced by predation
REAL(SPP) :: FIntStrPDP = 0.0_SPP ! Fraction of structural dissolved organic phosphorus produced by predation
REAL(SPP) :: FIntStrPLP = 0.0_SPP ! Fraction of structural labile phosphorus produced by predation
REAL(SPP) :: FIntStrPRP = 0.0_SPP ! Fraction of structural refractory phosphorus produced by predation
!------------------------------------------------------------------------------
! Maximum uptake rate of nitrogen by diatoms.
! Units: g N / (g algal C) / s
!------------------------------------------------------------------------------
REAL(SPP) :: UpNMaxD = 2.315E-06_SPP

!------------------------------------------------------------------------------
! Maximum uptake rate of nitrogen by greens.
! Units: g N / (g algal C) / s
!------------------------------------------------------------------------------
REAL(SPP) :: UpNMaxG = 2.315E-06_SPP

!------------------------------------------------------------------------------
! Maximum uptake rate of phosphorus by diatoms.
! Units: g P / (g algal C) / s
!------------------------------------------------------------------------------
REAL(SPP) :: UpPMaxD = 2.315E-06_SPP

!------------------------------------------------------------------------------
! Maximum uptake rate of phosphorus by greens.
! Units: g P / (g algal C) / s
!------------------------------------------------------------------------------
REAL(SPP) :: UpPMaxG = 2.315E-06_SPP



!------------------------------------------------------------------------------
!  Eutro parameter declarations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! ALPHA_DIA = slope of diatoms light saturation curve at low light levels.
!             Its units are gC [gChl-a]^-1 h^-1 [umol quanta/(m^2 s)]^-1.
!------------------------------------------------------------------------------
REAL(SPP) :: ALPHA_DIA = 0.035_SPP

!------------------------------------------------------------------------------
! ALPHA_GRE = slope of greens light saturation curve at low light levels.
!             Its units are gC [gChl-a]^-1 h^-1 [umol quanta/(m^2 s)]^-1.
!------------------------------------------------------------------------------
REAL(SPP) :: ALPHA_GRE = 0.035_SPP

REAL(SPP) :: ANCP                         ! N:C ratio for plankton [-]
                                              
REAL(SPP) :: APCP                         ! P:C ratio for plankton [-]
                                              
REAL(SPP) :: ASCD                         ! Si:C ratio for diatoms [-]
                                              
REAL(SPP) :: AVFRAC                       ! Fraction of DOP available
                                              ! to phytoplankton.

!------------------------------------------------------------------------------
! Fraction of DON available to phytoplankton.
!------------------------------------------------------------------------------
REAL(SPP) :: AVFRACDON = 0.05_SPP
                                              

REAL(SPP) :: BMRD                        ! Diatom base metabolic rates [s^-1]
                                              
REAL(SPP) :: BMRG                        ! Green base metabolic rates [s^-1]
                                              
REAL(SPP) :: CCHLD                       ! Carbon:chlorophyll ratio for diatoms [-]
                                              
REAL(SPP) :: CCHLG                       ! Carbon:chlorophyll ratio for green algae [-]
                                              
REAL(SPP) :: CGZ                         ! Zooplankton grazing coefficient [sec^-1]

REAL(SPP) :: FCDD                        ! Fraction of basal metabolism exuded
                                              ! as dissolved organic carbon by diatoms [-] 

REAL(SPP) :: FCDG                        ! Fraction of basal metabolism exuded
                                              ! as dissolved organic carbon by greens [-]

REAL(SPP) :: FCDP                        ! Fraction of dissolved organic carbon
                                              ! produced by predation [-]

REAL(SPP) :: FCDZ                        ! Fraction of DOC from 
                                              ! zooplankton mortality [-]

REAL(SPP) :: FCLD                        ! Fraction of labile particulate organic
                                              ! carbon produced by diatoms metabolism.

REAL(SPP) :: FCLG                        ! Fraction of labile particulate organic
                                              ! carbon produced by greens metabolism.

REAL(SPP) :: FCLP                        ! Fraction of labile pdc from predation [-]

REAL(SPP) :: FCLZ                        ! Fraction of labile pdc from 
                                              ! zooplankton mortality [-]

REAL(SPP) :: FCRD                        ! Fraction of refractory particulate organic
                                              ! carbon produced by diatoms metabolism.

REAL(SPP) :: FCRG                        ! Fraction of refractory particulate organic
                                              ! carbon produced by greens metabolism.

REAL(SPP) :: FCRP                        ! Fraction of refractory pdc from
                                              ! predation [-]

REAL(SPP) :: FCRZ                        ! Fraction of refractory pdc from 
                                              ! zooplankton mortality [-]

REAL(SPP) :: FNDD                        ! Fraction dissolved organic nitrogen
                                              ! produced by diatom metabolism [-]

REAL(SPP) :: FNDG                        ! Fraction dissolved organic nitrogen
                                              ! produced by green algae metabolism [-]

REAL(SPP) :: FNDP                        ! Fraction of dissolved organic
                                              ! nitrogen produced by predation [-]              

REAL(SPP) :: FNDZ                        ! Fraction of dissolved organic
                                              ! nitrogen produced by zoo. mortality [-]              

REAL(SPP) :: FNID                        ! Fraction of inorganic nitrogren
                                              ! produced by diatom metabolism [-]

REAL(SPP) :: FNIG                        ! Fraction of inorganic nitrogren
                                              ! produced by green algae metabolism [-]

REAL(SPP) :: FNIP                        ! Fraction of inorganic nitrogen
                                              ! produced by predation [-]

REAL(SPP) :: FNIZ                        ! Fraction of inorganic nitrogen
                                              ! produced by zoo. mortality [-]

REAL(SPP) :: FNLD                        ! Fraction labile particulate nitrogen
                                              ! produced by diatom metabolism [-]

REAL(SPP) :: FNLG                        ! Fraction labile particulate nitrogen
                                              ! produced by green algae metabolism [-]

REAL(SPP) :: FNLP                        ! Fraction of labile organic
                                              ! nitrogen produced by predation [-]             

REAL(SPP) :: FNLZ                        ! Fraction of labile organic
                                              ! nitrogen produced by zoo. mortality [-]             

REAL(SPP) :: FNRD                        ! Fraction refractory particulate nitrogen
                                              ! produced by diatom metabolism [-]

REAL(SPP) :: FNRG                        ! Fraction refractory particulate nitrogen
                                              ! produced by green algae metabolism [-]

REAL(SPP) :: FNRP                        ! Fraction of refractory organic
                                              ! nitrogen produced by predation [-]            

REAL(SPP) :: FNRZ                        ! Fraction of refractory organic nitrogen
                                              ! produced by zooplankton mortality [-]

REAL(SPP) :: FPDD                        ! Fraction dissolved organic phosphorous
                                              ! produced by diatom metabolism [-]

REAL(SPP) :: FPDG                        ! Fraction dissolved organic phosphorous
                                              ! produced by green metabolism [-]

REAL(SPP) :: FPDP                        ! Fraction of dissolved organic phosphorous
                                              ! produced by predation [-]

REAL(SPP) :: FPDZ                        ! Fraction of dissolved organic
                                              ! phosphorous produced by zoo. mortality [-]

REAL(SPP) :: FPID                        ! Fraction of inorganic phosphorous
                                              ! produced by diatom metabolism [-]

REAL(SPP) :: FPIG                        ! Fraction of inorganic phosphorous
                                              ! produced by green metabolism [-]

REAL(SPP) :: FPIP                        ! Fraction of inorganic phosphorous
                                              ! produced by diatom metabolism [-]

REAL(SPP) :: FPIZ                        ! Fraction of inorganic phosphorous
                                              ! produced by zooplankton mortality [-]

REAL(SPP) :: FPLD                        ! Fraction of labile organic
                                              ! phosphorous produced by diatom 
                                              ! metabolism [-]

REAL(SPP) :: FPLG                        ! Fraction of labile organic
                                              ! phosphorous produced by green 
                                              ! algae metabolism [-]

REAL(SPP) :: FPLP                        ! Fraction of labile organic
                                              ! phosphorous produced by predation [-]

REAL(SPP) :: FPLZ                        ! Fraction of labile organic phosphorous
                                              ! produced by zoo. mortality [-]

REAL(SPP) :: FPRD                        ! Fraction of refractory organic
                                              ! phosphorous produced by diatom metabolism [-]

REAL(SPP) :: FPRG                        ! Fraction of refractory organic
                                              ! phosphorous produced by greeen 
                                              ! algae metabolism [-]

REAL(SPP) :: FPRP                        ! Fraction of refractory organic
                                              ! phosphorous produced by predation [-]

REAL(SPP) :: FPRZ                        ! Fraction of refractory organic
                                              ! phosphorous produced by zoo. mortality [-]

REAL(SPP) :: FSAP                        ! Fraction of silica made available
                                              ! through predation [-]

REAL(SPP) :: GREFF                       ! Zooplankton grazing efficiency [-]

REAL(SPP) :: ILMUL                       ! Scaling for surface illumination [-]

REAL(SPP) :: KDC                         ! Minimum DOC mineralization rate [1/sec]

REAL(SPP) :: KDCALG                      ! Mineralization algal dependence factor [m^3-kg^-1-s^1]

REAL(SPP) :: KDN                         ! Minimum dissolved organic nitrogen
                                              ! mineralization rate [s^-1]

REAL(SPP) :: KDNALG                      ! Dissolved organic nitrogen mineralization
                                              ! algal proportionality dependence [m^3-kg^-1-s^-1]

REAL(SPP) :: KDP                         ! Minimum dissolved organic phosphorous
                                              ! mineralization rate [s^-1]

REAL(SPP) :: KDPALG                      ! Dissolved phosphorous mineralization
                                              ! algal proportionality dependence [m^3-kg^-1-s^-1] 

REAL(SPP) :: KE                          ! Background light attenuation

REAL(SPP) :: KECHL                       ! Light attenuation factor for chlorophyll-a [m^2-kg^-1]

REAL(SPP) :: KHN                         ! Organic nitrogen decay
                                              ! half saturation constant [kg-m^-3]

REAL(SPP) :: KHND                        ! Mean half-saturation constant
                                              ! for nitrogen uptake by diatoms [kg-m^-3]

REAL(SPP) :: KHNG                        ! Mean half-saturation constant
                                              ! for nitrogen uptake by greens [kg-m^-3]

REAL(SPP) :: KHNNT                       ! Half saturation concentration of NH4
                                              ! required for nitrification [kg-m^-3]

REAL(SPP) :: KHP                         ! Organic phosphorous decay
                                              ! half saturation constant [kg-m^-3]

REAL(SPP) :: KHPD                        ! Mean half-saturation constant
                                              ! for algal phosphorous uptake [kg-m^-3]

REAL(SPP) :: KHPG                        ! Mean half-saturation constant
                                              ! for algal phosphorous uptake [kg-m^-3] 

REAL(SPP) :: KHSD                        ! Mean half-saturation constant
                                              ! for diatom silica uptake [kg-m^-3]

REAL(SPP) :: KLC                         ! Min. hydrolysis rate of labile PDC [s^-1]

REAL(SPP) :: KLCALG                      ! Labile hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

REAL(SPP) :: KLN                         ! Minimum labile particulate nitrogen
                                              ! hydrolysis rate [s^-1]

REAL(SPP) :: KLNALG                      ! Labile nitrogen hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

REAL(SPP) :: KLP                         ! Minimum labile organic phosphorous
                                              ! hydrolysis rate [s^-1]

REAL(SPP) :: KLPALG                      ! Labile phosphorous hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

REAL(SPP) :: KRC                         ! Min. hydrolysis rate of refractory PDC [s^-1]

REAL(SPP) :: KRCALG                      ! Refractory hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

REAL(SPP) :: KRN                         ! Minimum refractory particulate nitrogen
                                              ! hydrolysis rate [s^-1]

REAL(SPP) :: KRNALG                      ! Refractory nitrogen hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

REAL(SPP) :: KRP                         ! Minimum refractory organic phosphorous
                                              ! hydrolysis rate [s^-1]

REAL(SPP) :: KRPALG                      ! Refractory phosphorous hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

REAL(SPP) :: KSUA                        ! Particulate silica dissolution
                                              ! rate constant [s^-1]

REAL(SPP) :: KSZ                         ! half saturation coeff. of zooplankton  
                                              ! for phytoplankton [kg-m^-3]

REAL(SPP) :: KTBD                        ! Metabolism temperature dependence
                                              ! factor for diatoms [deg-C]

REAL(SPP) :: KTBG                        ! Metabolism temperature dependence
                                              ! factor for greens [deg-C]

REAL(SPP) :: KTGD1                       ! Effect of temperature below optimal
                                              ! temperature for diatoms [deg-C^2]

REAL(SPP) :: KTGD2                       ! Effect of temperature above optimal
                                              ! temperature for diatoms [deg-C^2]

REAL(SPP) :: KTGG1                       ! Effect of temperature below optimal
                                              ! temperature for greens [deg-C^2]

REAL(SPP) :: KTGG2                       ! Effect of temperature above optimal
                                              ! temperature for greens [deg-C^2]

REAL(SPP) :: KTHDR                       ! Hydrolysis temperature dependence 
                                              ! [deg-C^-1]

REAL(SPP) :: KTMNL                       ! Mineralization temperature 
                                              ! dependence [deg-C^-1]

REAL(SPP) :: KTNT1                       ! Effect of temperature below optimal T
                                              ! for nitrification [deg-C^-2]

REAL(SPP) :: KTNT2                       ! Effect of temperature above optimal T
                                              ! for nitrification [deg-C^-2]

REAL(SPP) :: KTSUA                       ! Silica dissolution temperature
                                              ! rate constant [deg-C^-1

REAL(SPP) :: NTM                         ! Nitrification rate coefficient at
                                              ! optimal temperature [kg-m^-3-s-1]


!------------------------------------------------------------------------------
! PBMAX_DIA = specific diatoms photosynthetic rate at optimal 
!             illumination, i.e the assimilation number. It's the maximum
!             assimilation number normalized to diatoms biomass.
!             Its units are gC [gChl-a]^-1 h^-1.
!------------------------------------------------------------------------------
REAL(SPP) :: PBMAX_DIA = 10.5_SPP

!------------------------------------------------------------------------------
! PBMAX_GRE = specific greens photosynthetic rate at optimal 
!             illumination, i.e the assimilation number. It's the maximum
!             assimilation number normalized to greens biomass.
!             Its units are gC [gChl-a]^-1 h^-1.
!------------------------------------------------------------------------------
REAL(SPP) :: PBMAX_GRE = 10.5_SPP


REAL(SPP) :: PMD                         ! Diatom production under
                                              ! optimal conditions [s^-1]

REAL(SPP) :: PMG                         ! Green algae production under
                                              ! optimal conditions [s^-1]

REAL(SPP) :: TMD                         ! Temperature of optimum growth
                                              ! for diatoms [deg-C] 

REAL(SPP) :: TMG                         ! Temperature of optimum growth
                                              ! for greens [deg-C]

REAL(SPP) :: TMNT                        ! Optimal temperature for
                                              ! nitrification [deg-C]

REAL(SPP) :: TRD                         ! Metabolism reference temperature
                                              ! for diatoms [deg-C]

REAL(SPP) :: TRG                         ! Metabolism reference temperature
                                              ! for greens [deg-C]

REAL(SPP) :: TRHDR                       ! Reference temperature for
                                              ! hydrolysis [deg-C]

REAL(SPP) :: TRMNL                       ! Reference temperature for
                                              ! mineralization [deg-C]

REAL(SPP) :: TRSUA                       ! Silica dissolution reference
                                              ! temperature [deg-C]

REAL(SPP) :: TZREF                       ! Reference predation temperature

REAL(SPP) :: ZDTH                        ! Zooplankton death/die-off coefficient

REAL(SPP) :: ZTHET                       ! temperature coefficient 
                                              ! for predation [-]




!------------------------------------------------------------------------------
!  Default values of the eutro parameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Eutrophication parameters
!------------------------------------------------------------------------------

DATA ANCP/0.2_SPP/                       ! N:C ratio for plankton [-]
  
DATA APCP/0.01_SPP/                      ! P:C ratio for plankton [-] 
  
DATA ASCD/0.8_SPP/                       ! Si:C ratio for diatoms [-]
  
DATA AVFRAC/0.05_SPP/                   ! Fraction of DOP available to
                                              ! phytoplankton.

DATA BMRD/8.6E-07_SPP/                   ! Diatom base metabolic rates [s^-1]
  
DATA BMRG/8.6E-07_SPP/                   ! Green base metabolic rates [s^-1]
  
DATA CCHLD/50.0_SPP/                    ! Carbon:chlorophyll ratio for diatoms [-]
  
DATA CCHLG/50.0_SPP/                    ! Carbon:chlorophyll ratio for green algae [-]
  
DATA CGZ/4.6E-06_SPP/                     ! Zooplankton grazing coefficient [sec^-1]

DATA FCDD/0.1_SPP/                       ! Fraction of basal metabolism exuded
                                              ! as dissolved organic carbon by diatoms [-] 
  
DATA FCDG/0.1_SPP/                       ! Fraction of basal metabolism exuded
                                              ! as dissolved organic carbon by greens [-]
  
DATA FCDP/0.1_SPP/                       ! Fraction of dissolved organic carbon
                                              ! produced by predation [-]

DATA FCDZ/0.1_SPP/                       ! Fraction of DOC from
                                              ! zooplankton mortality [-]

DATA FCLD/0.10_SPP/                      ! Fraction of labile particulate
                                              ! organic carbon produced from
                                              ! diatoms metabolic and/or
                                              ! mortality rate.

DATA FCLG/0.10_SPP/                      ! Fraction of labile particulate
                                              ! organic carbon produced from
                                              ! greens metabolic and/or
                                              ! mortality rate.

DATA FCLP/0.55_SPP/                      ! Fraction of labile pdc 
                                              ! from predation [-]

DATA FCLZ/0.55_SPP/                      ! Fraction of labile pdc from
                                              ! zooplankton mortality [-]

DATA FCRD/0.10_SPP/                      ! Fraction of refractory 
                                              ! particulate organic carbon
                                              ! produced from diatoms   
                                              ! metabolic and/or mortality
                                              ! rate.

DATA FCRG/0.10_SPP/                      ! Fraction of refractory 
                                              ! particulate organic carbon
                                              ! produced from greens   
                                              ! metabolic and/or mortality
                                              ! rate.

DATA FCRP/0.35_SPP/                      ! Fraction of refractory pdc from
                                              ! predation [-]

DATA FCRZ/0.35_SPP/                      ! Fraction of refractory pdc from 
                                              ! zooplankton mortality [-]

DATA FNDD/0.8_SPP/                       ! Fraction dissolved organic nitrogen
                                              ! produced by diatom metabolism [-]

DATA FNDG/0.8_SPP/                       ! Fraction dissolved organic nitrogen
                                              ! produced by green algae metabolism [-]

DATA FNDP/0.1_SPP/                       ! Fraction of dissolved organic
                                              ! nitrogen produced by predation [-]              

DATA FNDZ/0.4_SPP/                       ! Fraction of dissolved organic
                                              ! nitrogen produced by zoo. mortality [-]              

DATA FNID/0.1_SPP/                       ! Fraction of inorganic nitrogren
                                              ! produced by diatom metabolism [-]

DATA FNIG/0.1_SPP/                       ! Fraction of inorganic nitrogren
                                              ! produced by green algae metabolism [-]

DATA FNIP/0.0_SPP/                       ! Fraction of inorganic nitrogen
                                              ! produced by predation [-]

DATA FNIZ/0.1_SPP/                       ! Fraction of inorganic nitrogen
                                              ! produced by zoo. mortality [-]

DATA FNLD/0.1_SPP/                       ! Fraction labile particulate nitrogen
                                              ! produced by diatom metabolism [-]

DATA FNLG/0.1_SPP/                       ! Fraction labile particulate nitrogen
                                              ! produced by green algae metabolism [-]

DATA FNLP/0.55_SPP/                      ! Fraction of labile organic
                                              ! nitrogen produced by predation [-]             

DATA FNLZ/0.35_SPP/                      ! Fraction of labile organic
                                              ! nitrogen produced by zoo. mortality [-]             

DATA FNRD/0.0_SPP/                       ! Fraction refractory particulate nitrogen
                                              ! produced by diatom metabolism [-]

DATA FNRG/0.0_SPP/                       ! Fraction refractory particulate nitrogen
                                              ! produced by green algae metabolism [-]

DATA FNRP/0.35_SPP/                      ! Fraction of refractory organic
                                              ! nitrogen produced by predation [-]            

DATA FNRZ/0.15_SPP/                      ! Fraction of refractory organic nitrogen
                                              ! produced by zooplankton mortality [-]

DATA FPDD/0.8_SPP/                       ! Fraction dissolved organic phosphorous
                                              ! produced by diatom metabolism [-]

DATA FPDG/0.8_SPP/                       ! Fraction dissolved organic phosphorous
                                              ! produced by green metabolism [-]

DATA FPDP/0.1_SPP/                       ! Fraction of dissolved organic phosphorous
                                              ! produced by predation [-]

DATA FPDZ/0.4_SPP/                       ! Fraction of dissolved organic 
                                              ! phosphorous produced by zoo. mortality [-]

DATA FPID/0.1_SPP/                       ! Fraction inorganic phosphorous
                                              ! produced by diatom metabolism [-]

DATA FPIG/0.1_SPP/                       ! Fraction inorganic phosphorous
                                              ! produced by green metabolism [-]

DATA FPIP/0.0_SPP/                       ! Fraction inorganic phosphorous
                                              ! produced by predation  [-]

DATA FPIZ/0.1_SPP/                       ! Fraction inorganic phosphorous
                                              ! produced by zooplankton mortality [-]

DATA FPLD/0.1_SPP/                       ! Fraction of labile organic
                                              ! phosphorous produced by diatom 
                                              ! metabolism [-]

DATA FPLG/0.1_SPP/                       ! Fraction of labile organic
                                              ! phosphorous produced by green 
                                              ! algae metabolism [-]

DATA FPLP/0.55_SPP/                      ! Fraction of labile organic
                                              ! phosphorous produced by predation [-]

DATA FPLZ/0.35_SPP/                      ! Fraction of labile organic phosphorous
                                              ! produced by zoo. mortality [-]

DATA FPRD/0.0_SPP/                       ! Fraction of refractory organic
                                              ! phosphorous produced by diatom 
                                              ! metabolism [-]

DATA FPRG/0.0_SPP/                       ! Fraction of refractory organic
                                              ! phosphorous produced by greeen 
                                              ! algae metabolism [-]

DATA FPRP/0.35_SPP/                      ! Fraction of refractory organic
                                              ! phosphorous produced by predation [-]

DATA FPRZ/0.15_SPP/                      ! Fraction of refractory organic
                                              ! phosphorous produced by zoo. mortality [-]

DATA FSAP/0.0_SPP/                       ! Fraction silica made available
                                              ! through predation [-]


DATA GREFF/0.6_SPP/                     ! Zooplankton grazing efficiency [-]

DATA ILMUL/1.0_SPP/                     ! Scaling for surface illumination [-]


DATA KDC/1.2E-07/                     ! Respiration rate of DOC [sec^-1]

DATA KDCALG/0.0_SPP/                   ! Mineralization algal dependence factor [m^3-kg^-1-s^1]

DATA KDN/2.3E-07/                     ! Minimum dissolved organic nitrogen
                                              ! mineralization rate [sec^-1]

DATA KDNALG/0.0_SPP/                   ! Minimum dissolved organic nitrogen
                                              ! mineralization rate [m^3-kg^-1-s^1]

DATA KDP/5.0E-07_SPP/                     ! Minimum dissolved organic phosphorous
                                              ! mineralization rate [s^-1]

DATA KDPALG/2.3E-03_SPP/               ! Dissolved phosphorous mineralization
                                              ! algal proportionality dependence [m^3-kg^-1-s^-1] 

DATA KE/0.15_SPP/                      ! Background light attenuation [m^-1]

DATA KECHL/1.7E+04_SPP/                ! Light attenuation factor for chlorophyll-a [m^2-kg^-1]

DATA KHN/2.5E-05_SPP/                     ! Organic nitrogen decay
                                              ! half saturation constant [kg-m^-3]

DATA KHND/2.5E-05_SPP/                   ! Mean half-saturation constant
                                              ! for nitrogen uptake by diatoms [kg-m^-3]

DATA KHNG/2.5E-05_SPP/                   ! Mean half-saturation constant
                                              ! for nitrogen uptake by greens [kg-m^-3]

DATA KHNNT/1.0E-04_SPP/                 ! Half saturation concentration of NH4
                                              ! required for nitrification [kg-m^-3]

DATA KHP/1.0E-06_SPP/                     ! Organic phosphorous decay
                                              ! half saturation constant [kg-m^-3]

DATA KHPD/1.0E-06_SPP/                   ! Mean half-saturation constant
                                              ! for algal phosphorous uptake [kg-m^-3]

DATA KHPG/1.0E-06_SPP/                   ! Mean half-saturation constant
                                              ! for algal phosphorous uptake [kg-m^-3] 

DATA KHSD/3.0E-05_SPP/                   ! Mean half-saturation constant
                                              ! for diatom silica uptake [kg-m^-3]

DATA KLC/8.6E-07_SPP/                     ! Hydrolysis rate of labile PDC [s^-1]

DATA KLCALG/0.0_SPP/                   ! Labile hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

DATA KLN/2.3E-07_SPP/                     ! Minimum labile particulate nitrogen
                                              ! hydrolysis rate [s^-1]

DATA KLNALG/0.0_SPP/                   ! Labile nitrogen hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

DATA KLP/2.3E-07_SPP/                    ! Minimum labile organic phosphorous
                                              ! hydrolysis rate [s^-1]

DATA KLPALG/2.3E-03_SPP/                   ! Labile phosphorous hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

DATA KRC/5.0E-08_SPP/                     ! Min. hydrolysis rate of refractory PDC 
                                              ! [s^-1]

DATA KRCALG/0.0_SPP/                   ! Refractory hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

DATA KRN/5.7E-08_SPP/                     ! Minimum refractory particulate nitrogen
                                              ! hydrolysis rate [s^-1]

DATA KRNALG/0.0_SPP/                   ! Refractory nitrogen hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

DATA KRP/5.7E-08_SPP/                     ! Minimum labile organic phosphorous
                                              ! hydrolysis rate [s^-1]

DATA KRPALG/2.3E-03_SPP/               ! Refractory phosphorous hydrol. algal proportionality
                                              ! dependence [m^3-kg^-1-s^-1]

DATA KSUA/5.0E-08_SPP/                   ! Particulate silica dissolution
                                              ! rate constant [s^-1]

DATA KSZ/7.5E-04_SPP/                     ! half saturation coeff. of zooplankton 
                                              ! for phytoplankton [kg-m^-3]
 
DATA KTBD/6.9E-02_SPP/                   ! Metabolism temperature dependence
                                              ! factor for diatoms [deg-C]

DATA KTBG/6.9E-02_SPP/                   ! Metabolism temperature dependence
                                              ! factor for greens [deg-C]

DATA KTGD1/2.5E-03_SPP/                 ! Effect of temperature below optimal
                                              ! temperature for diatoms [deg-C^2]

DATA KTGD2/6.0E-03_SPP/                 ! Effect of temperature above optimal
                                              ! temperature for diatoms [deg-C^2]

DATA KTGG1/2.5E-03_SPP/                 ! Effect of temperature below optimal
                                              ! temperature for greens [deg-C^2]

DATA KTGG2/1.0E-02_SPP/                 ! Effect of temperature above optimal
                                              ! temperature for greens [deg-C^2]

DATA KTHDR/6.9E-02_SPP/                 ! Hydrolysis temperature dependence 
                                              ! [deg-C^-1]

DATA KTMNL/6.9E-02_SPP/                 ! Mineralization temperature 
                                              ! dependence [deg-C^-1]

DATA KTNT1/9.0E-02_SPP/                 ! Effect of temperature below optimal T
                                              ! for nitrification [deg-C^-2]

DATA KTNT2/9.0E-02_SPP/                 ! Effect of temperature above optimal T
                                              ! for nitrification [deg-C^-2]

DATA KTSUA/1.1E-06_SPP/                 ! Silica dissolution temperature
                                              ! rate constant [deg-C^-1]

DATA NTM/1.2E-12_SPP/                     ! Nitrification rate coefficient at
                                              ! optimal temperature [kg-m^-3-s-1]

DATA PMD/2.4E-05_SPP/                     ! Diatom production under
                                              ! optimal conditions [s^-1]

DATA PMG/2.3E-05_SPP/                     ! Green algae production under
                                              ! optimal conditions [s^-1]

DATA TMD/28.0_SPP/                        ! Temperature of optimum growth
                                              ! for diatoms [deg-C] 

DATA TMG/28.0_SPP/                        ! Temperature of optimum growth
                                              ! for greens [deg-C]

DATA TMNT/25.0_SPP/                      ! Optimal temperature for
                                              ! nitrification [deg-C]

DATA TRD/25.0_SPP/                        ! Metabolism reference temperature
                                              ! for diatoms [deg-C]

DATA TRG/25.0_SPP/                        ! Metabolism reference temperature
                                              ! for greens [deg-C]

DATA TRHDR/25.0_SPP/                    ! Reference temperature for
                                              ! hydrolysis [deg-C]

DATA TRMNL/25.0_SPP/                    ! Reference temperature for
                                              ! mineralization [deg-C]

DATA TRSUA/25.0_SPP/                    ! Silica dissolution reference
                                              ! temperature [deg-C]

DATA TZREF/25.0_SPP/                    ! Predation reference 
                                              ! temperature [deg-C]

DATA ZDTH/3.5E-07_SPP/                   ! Zooplankton death/die-off 
                                              ! coefficient [sec^-1]

DATA ZTHET/1.05_SPP/                    ! temperature coefficient 
                                              ! for predation [-]


!------------------------------------------------------------------------------
!  Dissolved Oxygen parameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Maximum observed carbon sod mass flux.
!  It has units of (kg of O2)/(m^2 s).
!------------------------------------------------------------------------------
REAL(SPP) :: C_SOD_FLUX_MAX = 9.259E-09_SPP

!------------------------------------------------------------------------------
!  Denitrification carbon-to-nitrogen ratio
!------------------------------------------------------------------------------
REAL(SPP) :: DENIT_CN_RATIO = 1.07_SPP

!------------------------------------------------------------------------------
!  COD rate coefficient
!------------------------------------------------------------------------------
REAL(SPP) :: KCOD = 1.0_SPP

!------------------------------------------------------------------------------
!  Maximum denitrification rate coefficient [s^-1]
!------------------------------------------------------------------------------
REAL(SPP) :: KDENITR = 1.0_SPP

!------------------------------------------------------------------------------
!  Half-Saturation concentration of NO3 required for 
!  denitrification [kg NO3 / m^-3]
!------------------------------------------------------------------------------
REAL(SPP) :: KHDENITR = 1.0_SPP

!------------------------------------------------------------------------------
!  Half-Saturation concentration of NO3 required for
!  denitrification in the sediment bed [kg NO3 / m^-3].
!------------------------------------------------------------------------------
REAL(SPP) :: KHDENITR_SED = 1.0E-03_SPP

!------------------------------------------------------------------------------
!  Half-Saturation concentration of dissolved oxygen required for 
!  nitrification [kg O2 / m^-3]
!------------------------------------------------------------------------------
REAL(SPP) :: KHDONT = 1.0_SPP

!------------------------------------------------------------------------------
!  Half-Saturation concentration of dissolved oxygen required for
!  nitrification [kg O2 / m^-3] in the sediment bed.
!------------------------------------------------------------------------------
REAL(SPP) :: KHDONT_SED = 2.0E-3_SPP

!------------------------------------------------------------------------------
!  Half-Saturation concentration of dissolved oxygen required for oxic
!  respiration [kg O2 / m^-3]
!------------------------------------------------------------------------------
REAL(SPP) :: KHODOC = 1.0_SPP

!------------------------------------------------------------------------------
!  Half-Saturation concentration of dissolved oxygen required for oxic
!  respiration [kg O2 / m^-3] in the sediment bed.
!------------------------------------------------------------------------------
REAL(SPP) :: KHODOC_SED = 2.0E-03_SPP

!------------------------------------------------------------------------------
!  Half-Saturation concentration of dissolved oxygen required for
!  exertion of chemical oxygen demand [kg O2 / m^-3]
!------------------------------------------------------------------------------
REAL(SPP) :: KHOCOD = 1.0_SPP

!------------------------------------------------------------------------------
! D.O. reaeration velocity (m/s)
!------------------------------------------------------------------------------
REAL(SPP) :: KREA_VEL_DO = 1.1574E-05_SPP


!------------------------------------------------------------------------------
!  Minimum observed NO3 mass flux across the water-sediment interface.
!  The number is negative because it represents a mass flux from the water
!  to the sediment bed.
!  It has units of (kg of N2)/(m^2 s)
!------------------------------------------------------------------------------
REAL(SPP) :: SED_NO3_FLUX_MIN = -1.62118E-10_SPP


!------------------------------------------------------------------------------
!  Dissolved Oxygen-to-Carbon ratio in respiration = 2.67 grams O2/grams C
!------------------------------------------------------------------------------
REAL(SPP) :: RCDO = 2.67_SPP

!------------------------------------------------------------------------------
!  Nitrogen to Dissolved Oxygen ratio
!------------------------------------------------------------------------------
REAL(SPP) :: RNTO = 1.0_SPP

!------------------------------------------------------------------------------
! D.O. reaeration coefficient (m/s)
!------------------------------------------------------------------------------
REAL(SPP) :: KRDO = 2.77E-05_SPP

!------------------------------------------------------------------------------
! SOD to nitrification ratio?? 
!------------------------------------------------------------------------------
REAL(SPP) :: RSODNTR = 0.24_SPP

!------------------------------------------------------------------------------
!  SOD reference temperature (degrees Celsius)
!------------------------------------------------------------------------------
REAL(SPP) :: TREF = 25.0_SPP

!------------------------------------------------------------------------------
!  Dissolved Oxygen arrays used to store kinetic processes.
!------------------------------------------------------------------------------
!Calculated in dissolved oxygen for output purposes
!REAL(SPP), ALLOCATABLE :: NITDO2(:)
!Calculated in dissolved oxygen for output purposes
!REAL(SPP), ALLOCATABLE :: DOMETD_ARR(:)
!Calculated in dissolved oxygen for output purposes
!REAL(SPP), ALLOCATABLE :: DOMETG_ARR(:)
!Calculated in dissolved oxygen for output purposes
!REAL(SPP), ALLOCATABLE :: DOPREDD_ARR(:)
!Calculated in dissolved oxygen for output purposes
!REAL(SPP), ALLOCATABLE :: DOPREDG_ARR(:)
!Calculated in dissolved oxygen for output purposes
!REAL(SPP), ALLOCATABLE :: DOZOO_ARR(:)
!Calculated in dissolved oxygen for output purposes
!REAL(SPP), ALLOCATABLE :: DOMNLDOC_ARR(:)
 

!
!------------------------------------------------------------------------------
!
contains

SUBROUTINE EUT_allocate
!------------------------------------------------------------------------------
!-
!-  $Id: eut.F90,v 1.0.6.1 2014/08/26 22:54:04 wmelende Exp $
!-  $Locker:  $
!-
!-  Purpose and Methods : Include file for eutrophication variables.
!-
!-  Created: 12/09/97  M. Settles
!------------------------------------------------------------------------------

IMPLICIT NONE


!------------------------------------------------------------------------------
!  Miscellaneous arrays/variables; not parameters
!------------------------------------------------------------------------------
ALLOCATE(DTM(KBV,NTT))         ;  DTM      =  0.0_SPP        ! DTM

  
ALLOCATE(BMD(KBV))        ! Base metabolic rate for diatoms

ALLOCATE(BMG(KBV))        ! Base metabolic rate for greens

ALLOCATE(FTHDR(KBV))        ! Temperature correction for hydrolysis

ALLOCATE(FTMNL(KBV))        ! Temperature correction for mineralization

ALLOCATE(OUTPUT_PAR(KBV))     ! Daily cumulative PAR
OUTPUT_PAR = 0._SPP


ALLOCATE(PD(KBV))        ! Production rates for diatoms
PD = 0._SPP
ALLOCATE(PD_AVG(KBV))        ! Average production rates for 
PD_AVG = 0._SPP                                              ! diatoms

ALLOCATE(PG(KBV))        ! Production rates for greens
PG = 0._SPP
ALLOCATE(PG_AVG(KBV))        ! Average production rates for 
PG_AVG = 0._SPP                                              ! greens

ALLOCATE(PRD(KBV))        ! Predation time derivative for diatoms

ALLOCATE(PRG(KBV))        ! Predation time derivative for diatoms


!------------------------------------------------------------------------------
!  Fine segmentation growth limitation factors
!------------------------------------------------------------------------------

ALLOCATE(PFD(KBV))     ! Phosphorous limitation for diatoms

ALLOCATE(SFD(KBV))     ! Silica limitation for diatoms

ALLOCATE(NFD(KBV))     ! Nitrogen limitation for diatoms

ALLOCATE(IFD(KBV))     ! Light limitation for diatoms

ALLOCATE(TFD(KBV))     ! Temperature limitation for diatoms

ALLOCATE(PFG(KBV))     ! Phosphorous limitation for greens

ALLOCATE(NFG(KBV))     ! Nitrogen limitation for greens

ALLOCATE(IFG(KBV))     ! Light limitation for greens

ALLOCATE(TFG(KBV))     ! Temperature limitation for greens

!------------------------------------------------------------------------------
! Initialize arrays to zero.
!------------------------------------------------------------------------------
PRD = 0.0_SPP
PRG = 0.0_SPP
PFD = 0.0_SPP
SFD = 0.0_SPP
NFD = 0.0_SPP
IFD = 0.0_SPP
TFD = 0.0_SPP
PFG = 0.0_SPP
NFG = 0.0_SPP
IFG = 0.0_SPP
TFG = 0.0_SPP

!Droop
ALLOCATE(QND(KBV))  
ALLOCATE(QPD(KBV))  
ALLOCATE(QNG(KBV))  
ALLOCATE(QPG(KBV))  

!------------------------------------------------------------------------------
!  Dissolved Oxygen arrays used to store kinetic processes.
!------------------------------------------------------------------------------
!ALLOCATE(NITDO2(KBV)) 
!ALLOCATE(DOMETD_ARR(KBV)) 
!ALLOCATE(DOMETG_ARR(KBV))
!ALLOCATE(DOPREDD_ARR(KBV)) 
!ALLOCATE(DOPREDG_ARR(KBV))
!ALLOCATE(DOZOO_ARR(KBV)) 
!ALLOCATE(DOMNLDOC_ARR(KBV)) 
 
!DOMETD_ARR = 0.
!DOMETG_ARR = 0.
!DOPREDD_ARR = 0.
!DOPREDG_ARR = 0.
!DOZOO_ARR = 0.
!DOMNLDOC_ARR = 0.
!NITDO2 = 0.

!
!------------------------------------------------------------------------------
!
return
END SUBROUTINE EUT_allocate

END Module EUT
