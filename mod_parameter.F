   MODULE MOD_PARAMETER
!** This routine reads in the control parameters for the biological model.  
      USE MOD_1D
      USE MOD_DETRITUS
!      USE LIMS, ONLY: KB
      IMPLICIT NONE
      SAVE
      CHARACTER(LEN=50)              ::  BIONAME,            PARMINDEX
      CHARACTER(LEN=50), ALLOCATABLE ::  DETRITUS_NAME(:),   DETRITUS_UNIT(:)
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

      !ensure biological parameter file exists    
      inquire(file=trim(fname),exist=fexist)
      if(.not.fexist)then
        write(*,*)'biological parameter file: ',trim(fname),' does not exist'
        write(*,*)'stopping'
        stop
      end if

      OPEN (1, FILE=trim(fname),status='old')

!Read Detritus parameters
      READ(1,*) NND
            NTT=NND   !TOTAL VARIABLE NUMBER FOR VERTICAL MIXING
            IND=1

      ALLOCATE (DETRITUS_NAME(NND), DETRITUS_UNIT(NND))
      DO J=1,NTT
        READ(1,*)DETRITUS_NAME(J),DETRITUS_UNIT(J)
      END DO

       write(6,*) "Det name unit", DETRITUS_NAME(1), DETRITUS_UNIT(1)


!      DETRITUS_NAME(1) = "TP"
!      DETRITUS_UNIT(1) = "mg/L"

!*******************    BIO_NAME ASSIGNMENT FOR NETCDF OUT PUT   ***************
            ALLOCATE (BIO_NAME(NTT,4))
            BIO_NAME(:,3)  = 'Concentration'
            BIO_NAME(:,4)  = 'Array'
            DO J=1,NTT
              BIO_NAME(J,1)=DETRITUS_NAME(J)
              BIO_NAME(J,2)=DETRITUS_UNIT(J)
            END DO

!************* ALLOCATING                   GLOBAL PARAMETERS      AND INITIALIZAITON
      ALLOCATE (BIO_VAR(KBV,NTT))        ;  BIO_VAR     =  0.0_SPP  !ENSEMBLE VARIABLES FOR MIXING
      ALLOCATE (DELTA_D(KBV))            ;  DELTA_D     =  1.0_SPP  !LAYER THICKNESS
      ALLOCATE (DELTA_Z(KBV))            ;  DELTA_Z     =  1.0_SPP  !THICKNESS BETWEEN LAYER CENTRE
      ALLOCATE (DEPTH_Z(KBV))            ;  DEPTH_Z     =  1.0_SPP  !DEPTH AT EACH LAYER CENTER
      ALLOCATE (KM_BIO(KBV))             ;  KM_BIO      =  0.0_SPP  !DIFFUSIVITY 
      ALLOCATE (T_BIO(KBV))              ;  T_BIO       =  0.0_SPP  !TEMPERATURE

!******************  ALLOCATING             DETRITUS       AND INITIALIZATION   *************************
   !Calculated
      ALLOCATE (BIO_D(KBV,NND))          ;  BIO_D       =  0.0_SPP  !DETRITUS ABUNDANCE
      ALLOCATE (DETRITE(KBV,NND))        ;  DETRITE     =  0.0_SPP  !DETRITUS CONCENTRATION 
      ALLOCATE (D_SINK(KBV,NND))         ;  D_SINK      =  0.0_SPP  !DETRITUS SINKING FLUX
      ALLOCATE (D_SINK1(KBV))            ;  D_SINK1     =  0.0_SPP  !SINKING FLUX WITHIN 1 ITERATION
      ALLOCATE (W_D_TOTAL(NND))          ;  W_D_TOTAL   =  0.0_SPP  !SINKING DISTANCE WITIN 1 T_STEP
   !Input
      ALLOCATE (D_0(NND))                ;  D_0         =  0.0_SPP  !DETRITUS THRESHOLD
      ALLOCATE (W_D(NND))                ;  W_D         =  0.0_SPP  !DETRITUS SINKING VELOCITY
    
   !Read inputs
      READ(1,*)  D_0 != 0.00
      READ(1,*)  W_D ! = 9.2e-7
       write(6,*) "D_0", D_0 
       write(6,*) "W_D",W_D

!      'D_0'         (NND)          mmol C m-3        DETRITUS THRESHOLLD
!      0.00
!      'W_D'         (NND)          m s-1             TP (as DETRITUS) SINKING
!      VELOCITY
!      9.2e-7                                         !(29 m/y as used by Chapra
!      and Dolan (2012)
      CLOSE(1)

       RETURN
END SUBROUTINE GET_PARAMETER_NEW


  SUBROUTINE BIO_PARAMETER_REPORT
    IMPLICIT NONE

  END SUBROUTINE BIO_PARAMETER_REPORT

END MODULE MOD_PARAMETER
