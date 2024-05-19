   MODULE MOD_PARAMETER
!** This routine reads in the control parameters for the biological model.  
!      USE mod_fvcom
      USE mod_1D

!      USE LIMS, ONLY: KB
      IMPLICIT NONE
      SAVE
      CHARACTER(LEN=50)              ::  BIONAME,            PARMINDEX
      CHARACTER(LEN=50), ALLOCATABLE ::  NUTRIENT_NAME(:),   NUTRIENT_UNIT(:)
!      CHARACTER(LEN=50)              ::  MODEL_STRUCTURE     
!      REAL(SPP)                      ::  TOTALL
      INTEGER                        ::  INN = 1


!*****************************  CONTAINS   ********************************
      CONTAINS
      SUBROUTINE GET_PARAMETER_NEW(fname)
      character(len=*)    :: fname
      logical             :: fexist
      INTEGER                        ::  I,J,K
!FV
      character(len=250)             :: inputfile 
      integer                        :: idum,jdum
      REAL(SPP) :: rdum
!************************    BEGINNING EXECUTABLE    **********************
      KBVM1=KBV-1
      call grid_setup(KBVM1)

      !ensure biological parameter file exists    
      inquire(file=trim(fname),exist=fexist)
      if(.not.fexist)then
        write(*,*)'biological parameter file: ',trim(fname),' does not exist'
        write(*,*)'stopping'
        stop
      end if

      OPEN (1, FILE=trim(fname),status='old')

!Read 'Nutrient' parameters
write(6,*) "read bio input, number of state variables"
      READ(1,*) NTT 
            INN=1
write(6,*) "nf",NTT
      CLOSE(1)

      call cgem_setup(NTT)
write(6,*) "after cgem_setup",NTT

      ALLOCATE (NUTRIENT_NAME(NTT), NUTRIENT_UNIT(NTT))

#ifdef DEBUG
write(6,*) "NTT",NTT
#endif

!      DO J=1,NTT
!       READ(1,*) NUTRIENT_NAME(J),NUTRIENT_UNIT(J)
!      END DO

   NUTRIENT_NAME(1) = "A1"
   NUTRIENT_UNIT(1) = "cells/m3" 
   NUTRIENT_NAME(2) = "Qn1"
   NUTRIENT_UNIT(2) = "mmol-N/cell"
   NUTRIENT_NAME(3) = "Qp1"
   NUTRIENT_UNIT(3) = "mmol-P/cell"
   NUTRIENT_NAME(4) = "Zoo1"
   NUTRIENT_UNIT(4) = "individuals/m3"
   NUTRIENT_NAME(5) = "Zoo2"
   NUTRIENT_UNIT(5) = "individuals/m3"
   NUTRIENT_NAME(6) = "NO3"
   NUTRIENT_UNIT(6) = "mmol-N/m3"
   NUTRIENT_NAME(7) = "NH4"
   NUTRIENT_UNIT(7) = "mmol-N/m3"
   NUTRIENT_NAME(8) = "PO4"
   NUTRIENT_UNIT(8) = "mmol-P/m3"
   NUTRIENT_NAME(9) = "DIC"
   NUTRIENT_UNIT(9) = "mmol-C/m3"
   NUTRIENT_NAME(10) = "O2"
   NUTRIENT_UNIT(10) = "mmol-O2/m3"
   NUTRIENT_NAME(11) = "OM1CA"
   NUTRIENT_UNIT(11) = "mmol-C/m3"
   NUTRIENT_NAME(12) = "OM1NA"
   NUTRIENT_UNIT(12) = "mmol-N/m3"
   NUTRIENT_NAME(13) = "OM1PA"
   NUTRIENT_UNIT(13) = "mmol-P/m3"
   NUTRIENT_NAME(14) = "OM2CA"
   NUTRIENT_UNIT(14) = "mmol-C/m3"
   NUTRIENT_NAME(15) = "OM2NA"
   NUTRIENT_UNIT(15) = "mmol-N/m3"
   NUTRIENT_NAME(16) = "OM2PA"
   NUTRIENT_UNIT(16) = "mmol-P/m3"
   NUTRIENT_NAME(17) = "OM1CZ"
   NUTRIENT_UNIT(17) = "mmol-C/m3"
   NUTRIENT_NAME(18) = "OM1NZ"
   NUTRIENT_UNIT(18) = "mmol-N/m3"
   NUTRIENT_NAME(19) = "OM1PZ"
   NUTRIENT_UNIT(19) = "mmol-P/m3"
   NUTRIENT_NAME(20) = "OM2CZ"
   NUTRIENT_UNIT(20) = "mmol-C/m3"
   NUTRIENT_NAME(21) = "OM2NZ"
   NUTRIENT_UNIT(21) = "mmol-N/m3"
   NUTRIENT_NAME(22) = "OM2PZ"
   NUTRIENT_UNIT(22) = "mmol-P/m3"
   NUTRIENT_NAME(23) = "OM1R"
   NUTRIENT_UNIT(23) = "mmol-C/m3"
   NUTRIENT_NAME(24) = "OM2R"
   NUTRIENT_UNIT(24) = "mmol-C/m3"
   NUTRIENT_NAME(25) = "CDOM"
   NUTRIENT_UNIT(25) = "ppb"
   NUTRIENT_NAME(26) = "Si"
   NUTRIENT_UNIT(26) = "mmol-Si/m3"
   NUTRIENT_NAME(27) = "OM1BC"
   NUTRIENT_UNIT(27) = "mmol-C/m3"
   NUTRIENT_NAME(28) = "OM2BC"
   NUTRIENT_UNIT(28) = "mmol-C/m3"
   NUTRIENT_NAME(29) = "Alk"
   NUTRIENT_UNIT(29) = "mmol-HCO3/m3"
   NUTRIENT_NAME(30) = "Tr"
   NUTRIENT_UNIT(30) = "none"

!#ifdef DEBUG
      DO J=1,NTT
       write(6,*) NUTRIENT_NAME(J),NUTRIENT_UNIT(J)
      END DO
!#endif


!*******************    BIO_NAME ASSIGNMENT FOR NETCDF OUT PUT   ***************
            ALLOCATE (BIO_NAME(NTT,4))
            BIO_NAME(:,3)  = 'Concentration'
            BIO_NAME(:,4)  = 'Array'
            DO J=1,NTT
              BIO_NAME(J,1)=NUTRIENT_NAME(J)
              BIO_NAME(J,2)=NUTRIENT_UNIT(J)
              write(6,*) BIO_NAME(J,1),BIO_NAME(J,2)
            END DO

!************* ALLOCATING                   GLOBAL PARAMETERS      AND INITIALIZAITON
      ALLOCATE (BIO_VAR(KBV,NTT))        ;  BIO_VAR     =  0.0_SPP  !ENSEMBLE VARIABLES FOR MIXING
      ALLOCATE (DELTA_D(KBV))            ;  DELTA_D     =  1.0_SPP  !LAYER THICKNESS
      ALLOCATE (DELTA_Z(KBV))            ;  DELTA_Z     =  1.0_SPP  !THICKNESS BETWEEN LAYER CENTRE
      ALLOCATE (DEPTH_Z(KBV))            ;  DEPTH_Z     =  1.0_SPP  !DEPTH AT EACH LAYER CENTER
      ALLOCATE (KM_BIO(KBV))             ;  KM_BIO      =  0.0_SPP  !DIFFUSIVITY 
      ALLOCATE (T_BIO(KBV))              ;  T_BIO       =  0.0_SPP  !TEMPERATURE
      ALLOCATE (S_BIO(KBV))              ;  S_BIO       =  0.0_SPP  !SALINITY


       RETURN
END SUBROUTINE GET_PARAMETER_NEW


  SUBROUTINE BIO_PARAMETER_REPORT
    IMPLICIT NONE

  END SUBROUTINE BIO_PARAMETER_REPORT

END MODULE MOD_PARAMETER
