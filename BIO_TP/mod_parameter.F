   MODULE MOD_PARAMETER
!** This routine reads in the control parameters for the biological model.  
      USE MOD_1D
!      USE LIMS, ONLY: KB
      IMPLICIT NONE
      SAVE
      CHARACTER(LEN=50)              ::  BIONAME,            PARMINDEX
      CHARACTER(LEN=50), ALLOCATABLE ::  NUTRIENT_NAME(:),   NUTRIENT_UNIT(:)
      CHARACTER(LEN=50)              ::  MODEL_STRUCTURE     
      REAL(SPP)                      ::  TOTALL


!*****************************  CONTAINS   ********************************
      CONTAINS
      SUBROUTINE GET_PARAMETER_NEW(fname)
      character(len=*)    :: fname
      logical             :: fexist
      INTEGER                        ::  I,J,K

!************************    BEGINNING EXECUTABLE    **********************
      KBVM1=KBV-1
      INN = 1

      !ensure biological parameter file exists    
      inquire(file=trim(fname),exist=fexist)
      if(.not.fexist)then
        write(*,*)'biological parameter file: ',trim(fname),' does not exist'
        write(*,*)'stopping'
        stop
      end if

      OPEN (1, FILE=trim(fname),status='old')

!Read parameters
      READ(1,*) NTT
            !INN=1

      ALLOCATE (NUTRIENT_NAME(NTT), NUTRIENT_UNIT(NTT))
!*******************    BIO_NAME ASSIGNMENT FOR NETCDF OUT PUT   ***************
            ALLOCATE (BIO_NAME(NTT,4))
            BIO_NAME(:,3)  = 'Concentration'
            BIO_NAME(:,4)  = 'Array'
            DO J=1,NTT
              BIO_NAME(J,1)=NUTRIENT_NAME(J)
              BIO_NAME(J,2)=NUTRIENT_UNIT(J)
            END DO

!************* ALLOCATING                   GLOBAL PARAMETERS      AND INITIALIZAITON
      ALLOCATE (BIO_VAR(KBV,NTT))        ;  BIO_VAR     =  0.0_SPP  !ENSEMBLE VARIABLES FOR MIXING
      ALLOCATE (DELTA_D(KBV))            ;  DELTA_D     =  1.0_SPP  !LAYER THICKNESS
      ALLOCATE (DELTA_Z(KBV))            ;  DELTA_Z     =  1.0_SPP  !THICKNESS BETWEEN LAYER CENTRE
      ALLOCATE (DEPTH_Z(KBV))            ;  DEPTH_Z     =  1.0_SPP  !DEPTH AT EACH LAYER CENTER
      ALLOCATE (KM_BIO(KBV))             ;  KM_BIO      =  0.0_SPP  !DIFFUSIVITY 
      ALLOCATE (T_BIO(KBV))              ;  T_BIO       =  0.0_SPP  !TEMPERATURE
      ALLOCATE (S_BIO(KBV))              ;  S_BIO       =  0.0_SPP  !SALINITY
      ALLOCATE (ISINK_OUT(NTT))         ;  ISINK_OUT  =  0        !0=no sink out, 1=sinks out
!******************  ALLOCATING             NUTRIENT       AND INITIALIZATION   *************************
   !Calculated
      !ALLOCATE (NUTRIENT(KBV,NTT))        ;  NUTRIENT     =  0.0_SPP  !NUTRIENT CONCENTRATION 
   !Input
      ALLOCATE (N_0(NTT))                ;  N_0         =  0.0_SPP  !NUTRIENT THRESHOLD
      ALLOCATE (W_N(NTT))                ;  W_N         =  0.0_SPP  !NUTRIENT SINKING VELOCITY
    
   !Read inputs
      DO J=1,NTT
        READ(1,*)  NUTRIENT_NAME(J),NUTRIENT_UNIT(J)
        READ(1,*)  N_0(J)       != 0.00 Nutrient threshold, (units of nutrient)
        READ(1,*)  W_N(J)       ! = TP (as NUTRIENT) sinking rate, (m s-1)
        READ(1,*)  ISINK_OUT(J) !1=sink out, 0=no sink out
      END DO

      write(6,*) "name",NUTRIENT_NAME
      write(6,*) "unit",NUTRIENT_UNIT
      write(6,*) "N_0", N_0 
      write(6,*) "W_N",W_N
      write(6,*) "ISINK_OUT",ISINK_OUT

      CLOSE(1)

       RETURN
END SUBROUTINE GET_PARAMETER_NEW


  SUBROUTINE BIO_PARAMETER_REPORT
    IMPLICIT NONE

  END SUBROUTINE BIO_PARAMETER_REPORT

END MODULE MOD_PARAMETER
