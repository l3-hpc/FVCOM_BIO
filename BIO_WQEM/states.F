MODULE STATES

!------------------------------------------------------------------------------
!-
!- $Id: states.F90,v 1.0.6.1 2014/08/26 22:54:04 wmelende Exp $
!- $Locker:  $
!-
!- Purpose and Methods:  Indices of the states for the fixed states
!-
!- Created:   08/11/98  M. Settles
!-
!- Revised: 06/10/10  W. Melendez :  Added tracer state index JTR.
!- Revised: 06/22/09  W. Melendez : ----> GOMDOM. This application is based
!-                                   on LM3-Eutro (LM3 version 3.2.13) and 
!-                                   its purpose is to serve as a dissolved
!-                                   oxygen model for the Gulf of Mexico
!-                                   hypoxia project. NetCDF is the library 
!-                                   used to handle input/output binary
!-                                   files.
!- Revised: 10/14/99  M. Settles  :  Add sediment/burial state
!-                                   lookup tables
!- Revised: 08/13/99  M. Settles  :  Add sand state variable
!- Revised: 07/15/99  M. Settles  :  Remove PCB/partitioning states
!- Revised: 04/07/99  M. Settles  :  New state variable index convention
!- Revised: 01/06/99  M. Settles  :  Add state index of last POC: JLASTPOC
!- Revised: 12/16/98  M. Settles  :  Add settling and resuspending states
!- Revised: 12/15/98  M. Settles  :  cpp away states when solids disabled
!- Revised: 12/15/98  M. Settles  :  Add atrazine state
!- Revised: 11/24/98  M. Settles  :  Remove algal phosphate states
!-
!------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER,PARAMETER  ::  JDOC  =    1        ! Dissolved organic Carbon
INTEGER,PARAMETER  ::  JDIA  =    2        ! Diatom
INTEGER,PARAMETER  ::  JGRE  =    3        ! Green
INTEGER,PARAMETER  ::  JZOO  =    4        ! Zooplankton
INTEGER,PARAMETER  ::  JLOC  =    5        ! Labile detritus
INTEGER,PARAMETER  ::  JROC  =    6        ! Refractory detritus

INTEGER,PARAMETER ::   JLASTPOC = JROC     ! State index of last particulate

INTEGER,PARAMETER  ::  JSRP   =  7         ! Inorganic phosphorous
                                    
INTEGER,PARAMETER  ::  JDOP   =  8         ! Dissolved organic phosphorous
INTEGER,PARAMETER  ::  JLOP   =  9         ! Labile detrital organic phosphorous
INTEGER,PARAMETER  ::  JROP   =  10        ! Refractory detrital organic phosphorous
                                    
INTEGER,PARAMETER :: JNH4     = 11         ! Ammonia
INTEGER,PARAMETER :: JNO3     = 12         ! Nitrate
INTEGER,PARAMETER :: JDON     = 13         ! Organic nitrogen
INTEGER,PARAMETER :: JLON     = 14         ! Labile detrital organic nitrogen
INTEGER,PARAMETER :: JRON     = 15         ! Refractory detrital organic nitrogen
                                    
INTEGER,PARAMETER  ::  JSA    = 16         ! Available silica
INTEGER,PARAMETER  ::  JSU    = 17         ! Unavailable silica

INTEGER,PARAMETER ::  JDO2   = 18          ! Dissolved oxygen

INTEGER,PARAMETER ::  JTR    = 19          ! Tracer 

INTEGER,PARAMETER  :: JDIAN = 20           ! Diatoms's Internal Nitrogen
INTEGER,PARAMETER  :: JDIAP = 21           ! Diatoms's Internal Phosphorus
INTEGER,PARAMETER  :: JGREN = 22           ! Greens's Internal Nitrogen
INTEGER,PARAMETER  :: JGREP = 23           ! Greens's Internal Phosphorus


!
!--------------------------------------------------------------
!

END MODULE STATES
