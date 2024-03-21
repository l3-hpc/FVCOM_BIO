MODULE FLAGS
!
!------------------------------------------------------------------------------
!-
!- $Id: flags.F90,v 1.0.6.1 2014/08/26 22:54:04 wmelende Exp $
!- $Locker:  $
!-
!- Purpose and Methods:  Miscellaneous program flags
!-
!- Created:   08/17/98  M. Settles
!-
!-
!- Revised: 02/03/14  W. Melendez :  Added SOC_LEHRTER and NO3SED_LEHRTER
!-                                   flags. Added NO3SED_ML and SOCSED_ML 
!-                                   flags. 
!- Revised: 04/30/13  W. Melendez :  Added OUTPUT_ZOOMORT flag.
!- Revised: 04/30/13  W. Melendez :  Added OUTPUT_DOCMINRL flag.
!- Revised: 04/25/13  W. Melendez :  Added WRITE_NITROGEN_PROCESSES flag.
!- Revised: 03/07/13  W. Melendez :  Added WRITE_MASSCHECK flag.
!- Revised: 02/21/13  W. Melendez :  Added WRITE_MASSBD flag.
!- Revised: 01/17/13  W. Melendez :  Added WRITE_DO2_KINETICS flag.
!- Revised: 06/28/12  W. Melendez :  Added DO_PHYTO_PROCESSES flag.
!- Revised: 06/27/12  W. Melendez :  Added DO_PRINT_HYDRO flag.
!- Revised: 05/24/12  W. Melendez :  Added KDWD and DO_KD_WIND_DEPTH flags.
!- Revised: 10/04/10  W. Melendez :  Added DO_WIND_SPEED flag.
!- Revised: 07/27/10  W. Melendez :  Initialized TITLE to empty strings.
!- Revised: 06/29/10  W. Melendez :  Added declaration and initialization of 
!-                                   OUTPUT_FILE_NAME character string.
!- Revised: 06/23/10  W. Melendez :  Added TITLE character array. This array
!-                                   is used to store user-input descriptions
!-                                   of what is in the model output files.
!- Revised: 06/10/10  W. Melendez :  Added GENERATE_GEOGRAPHICAL_OUTPUT 
!-                                   flag.  Added DO_TRACER flag.
!- Revised: 06/22/09  W. Melendez : ----> GOMDOM. This application is based
!-                                   on LM3-Eutro (LM3 version 3.2.13) and 
!-                                   its purpose is to serve as a dissolved
!-                                   oxygen model for the Gulf of Mexico
!-                                   hypoxia project. NetCDF is the library 
!-                                   used to handle input/output binary
!-                                   files.
!- Revised: 05/16/05  W. Melendez :  Added new flag  ---> 
!-                                   DORESCALE_PLOADS_ATM_TRIBS.
!- Revised: 12/03/04  W. Melendez :  Added new flag DOSED_LOADS_INCREASE.
!- Revised: 10/21/03  W. Melendez :  Added three new flags :
!-                                   DOTRIB_LOADS_ELIMINATION,
!-                                   DOATM_LOADS_ELIMINATION, and
!-                                   DOLOADS_INCREASE.
!- Revised: 10/09/03  W. Melendez :  Added new flag : DOEXP_DECAY_LOAD that 
!-                                   enables exponential decay on the
!-                                   phosphorus loads.  Made some cosmetic
!-                                   changes.
!- Revised: 07/30/03  W. Melendez :  Added new flag, LOADS_1982_83, that
!-                                   enables loads for 1982-83.
!- Revised: 09/23/02  W. Melendez :  Added boundary flows/dispersion
!-                                   flag, BD_FLOW_DISP.
!- Revised: 10/24/00  W. Melendez :  Add chloride flag
!- Revised: 06/08/00  M. Settles  :  Add DOL2PROD flag for level 2
!-                                   phytoplankton productivity load
!- Revised: 10/14/99  M. Settles  :  Add settling and sediment
!-                                   calculation flags
!- Revised: 10/01/99  W. Melendez :  Add resuspension flag
!- Revised: 08/13/99  M. Settles  :  Growth limitation and shear flags
!-                                   made logical
!- Revised: 07/19/99  M. Settles  :  Add silica limitation calc. flag
!- Revised: 07/14/99  M. Settles  :  Add eutrophication calculation flag
!- Revised: 06/18/99  M. Settles  :  Add shear calculation flag
!- Revised: 06/15/99  M. Settles  :  Add growth limitation flag
!- Revised: 04/13/99  M. Settles  :  Default parameter values
!- Revised: 04/13/99  M. Settles  :  Add former ICM control parameters
!- Revised: 03/01/99  M. Settles  :  Add DOAGGREG aggregation flag
!-
!------------------------------------------------------------------------------


IMPLICIT NONE

INTEGER :: DoDroop = 1

INTEGER :: OVERWRITE                     ! Overwrite output
LOGICAL :: END_RUN                       ! Stop the run
LOGICAL :: DOAGGREG                      ! Perform aggregation

INTEGER :: BCC,BCC__                     ! Boundary conditions flag
INTEGER :: BD_FLOW_DISP                  ! Boundary flows/disp. flag
INTEGER :: DIAC                          ! Diagnostics flag
INTEGER :: FLC                           ! Flow flag
INTEGER :: VBC                           ! Volume balance flag
INTEGER :: XYDFC                         ! Horizontal diffusion flag
INTEGER :: ZDFC                          ! Vertical diffusion flag

INTEGER :: CONST_TSTEP   = 0             ! Constant time-step flag

!------------------------------------------------------------------------------
! Wind-depth light attenuation coefficient flags
!------------------------------------------------------------------------------
INTEGER :: KDWD = 0

LOGICAL :: DO_KD_WIND_DEPTH = .FALSE.

!------------------------------------------------------------------------------
! Flag enabling write of growth limitation factors
!------------------------------------------------------------------------------
LOGICAL :: DOGRLM                        

LOGICAL :: DOSHEAR                     ! Flag enabling shear calculation
LOGICAL :: DOSETTLING                  ! Flag enabling settling
LOGICAL :: DOSEDIMENT                  ! Flag enabling sediments
LOGICAL :: DOTRANS                     ! Flag enabling water column transport
LOGICAL :: DOEUTRO                     ! Flag enabling eutrophication kinetics

!------------------------------------------------------------------------------
! Flag that enables the writing of phytoplankton related processes to
! output.
!------------------------------------------------------------------------------
LOGICAL :: DO_PHYTO_PROCESSES = .FALSE.

!------------------------------------------------------------------------------
! Flag that enables the writing of dissolved oxygen's sources and sinks to
! output.
!------------------------------------------------------------------------------
LOGICAL :: WRITE_DO2_KINETICS = .FALSE.

LOGICAL :: DO_SALINITY                 ! Flag enabling simulation of salinity
LOGICAL :: DO_DO2 = .TRUE.
LOGICAL :: DO_WOD = .FALSE.

!------------------------------------------------------------------------------
! Flag that enables the writing of mass exported/imported across the
! model grid boundary to output.
!------------------------------------------------------------------------------
LOGICAL :: WRITE_MASSBD = .FALSE.

!------------------------------------------------------------------------------
! Flag that enables the writing of mass created by the negative 
! concentration fix.
!------------------------------------------------------------------------------
LOGICAL :: WRITE_MASSCHECK = .FALSE.

!------------------------------------------------------------------------------
! Flag that enables the writing of denitrification and the sediment-water
! nitrate flux kinetic processes to output.
!------------------------------------------------------------------------------
LOGICAL :: WRITE_NITROGEN_PROCESSES = .FALSE.

!------------------------------------------------------------------------------
! Flag that enables the output of the DOC mass lost due to mineralization.
!------------------------------------------------------------------------------
LOGICAL :: OUTPUT_DOCMINRL = .FALSE.

!------------------------------------------------------------------------------
! Flag that enables the output of the amount of carbon lost from 
! zooplankton due to mortality.
!------------------------------------------------------------------------------
LOGICAL :: OUTPUT_ZOOMORT = .FALSE.

!------------------------------------------------------------------------------
! Flag enabling the use of wind speed in the program.
!------------------------------------------------------------------------------
LOGICAL :: DO_WIND_SPEED = .FALSE.

!------------------------------------------------------------------------------
! Flag enabling the simulation of a tracer (conservative substance)
!------------------------------------------------------------------------------
LOGICAL :: DO_TRACER = .FALSE.

!------------------------------------------------------------------------------
! Flag enabling mass balance calculations
!------------------------------------------------------------------------------
LOGICAL :: DOMASSBALANCE                 

!------------------------------------------------------------------------------
! Flag enabling generation of model output using geographical coordinates.
!------------------------------------------------------------------------------
LOGICAL :: GENERATE_GEOGRAPHICAL_OUTPUT = .FALSE.

!------------------------------------------------------------------------------
! Flag enabling writing of hydrodynamic data to output.
!------------------------------------------------------------------------------
LOGICAL :: DO_PRINT_HYDRO = .FALSE.

!------------------------------------------------------------------------------
! Flags enabling Murrell and Lehrter's empirical formulas used to 
! calculate nitrate and SOC fluxes across the water-sediment interface.
!------------------------------------------------------------------------------
INTEGER :: NO3SED_ML = 0
INTEGER :: SOCSED_ML = 0

LOGICAL :: NO3SED_LEHRTER = .FALSE.
LOGICAL :: SOC_LEHRTER = .FALSE.

!------------------------------------------------------------------------------
!  Silica limitation calculation flag
!------------------------------------------------------------------------------
INTEGER :: SILIM                 

DATA BCC/0/      
DATA BD_FLOW_DISP/1/                 
DATA DIAC/0/
DATA FLC/1/
DATA OVERWRITE/0/        
DATA VBC/1/
DATA XYDFC/1/
DATA ZDFC/1/

DATA DOGRLM/.FALSE./
DATA DOSHEAR/.FALSE./
DATA DOSETTLING/.FALSE./
DATA DOSEDIMENT/.FALSE./
DATA DOTRANS/.FALSE./
DATA DOEUTRO/.TRUE./
DATA DO_SALINITY/.FALSE./
DATA DOMASSBALANCE/.FALSE./

DATA SILIM/1/

!------------------------------------------------------------------------------
!  This character array stores user-input descriptions of what is in 
!  the model output files.
!------------------------------------------------------------------------------
CHARACTER (LEN = 90), DIMENSION (10) :: TITLE = ''


!------------------------------------------------------------------------------
!  This character string stores the user-input name of the model output file.
!------------------------------------------------------------------------------
CHARACTER (LEN = 90) :: OUTPUT_FILE_NAME = 'Model_Output.nc'


END MODULE FLAGS

