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


!*****************************  CONTAINS   ********************************
      CONTAINS
      SUBROUTINE GET_PARAMETER_NEW(fname)
      character(len=*)    :: fname
      logical             :: fexist
      INTEGER                        ::  I,J,K,INN
!FV
      character(len=250)             :: inputfile 
      integer                        :: idum,jdum
      REAL(SPP) :: rdum
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
write(6,*) "before read"
      READ(1,*) NTT 
            INN=1

write(6,*) "NTT",NTT



      ALLOCATE (NUTRIENT_NAME(NTT), NUTRIENT_UNIT(NTT))

#ifdef DEBUG
write(6,*) "NTT",NTT
#endif

!      DO J=1,NTT
!       READ(1,*) NUTRIENT_NAME(J),NUTRIENT_UNIT(J)
!      END DO

   NUTRIENT_NAME(1) = "DOC"
   NUTRIENT_UNIT(1) = "C kg/m3" 
   NUTRIENT_NAME(2) = "DIA"
   NUTRIENT_UNIT(2) = "C kg/m3"
   NUTRIENT_NAME(3) = "GRE"
   NUTRIENT_UNIT(3) = "C kg/m3"
   NUTRIENT_NAME(4) = "ZOO"
   NUTRIENT_UNIT(4) = "C kg/m3"
   NUTRIENT_NAME(5) = "LOC"
   NUTRIENT_UNIT(5) = "C kg/m3"
   NUTRIENT_NAME(6) = "ROC"
   NUTRIENT_UNIT(6) = "C kg/m3"
   NUTRIENT_NAME(7) = "SRP"
   NUTRIENT_UNIT(7) = "P kg/m3"
   NUTRIENT_NAME(8) = "DOP"
   NUTRIENT_UNIT(8) = "P kg/m3"
   NUTRIENT_NAME(9) = "LOP"
   NUTRIENT_UNIT(9) = "P kg/m3"
   NUTRIENT_NAME(10) = "ROP"
   NUTRIENT_UNIT(10) = "P kg/m3"
   NUTRIENT_NAME(11) = "NH4"
   NUTRIENT_UNIT(11) = "N kg/m3"
   NUTRIENT_NAME(12) = "NO3"
   NUTRIENT_UNIT(12) = "N kg/m3"
   NUTRIENT_NAME(13) = "DON"
   NUTRIENT_UNIT(13) = "N kg/m3"
   NUTRIENT_NAME(14) = "LON"
   NUTRIENT_UNIT(14) = "N kg/m3"
   NUTRIENT_NAME(15) = "RON"
   NUTRIENT_UNIT(15) = "N kg/m3"
   NUTRIENT_NAME(16) = "SA"
   NUTRIENT_UNIT(16) = "Si kg/m3"
   NUTRIENT_NAME(17) = "SU"
   NUTRIENT_UNIT(17) = "Si kg/m3"
   NUTRIENT_NAME(18) = "DO2"
   NUTRIENT_UNIT(18) = "O2 kg/m3"
   NUTRIENT_NAME(19) = "TR"
   NUTRIENT_UNIT(19) = "kg/m3"
!   NUTRIENT_NAME(20) = "JDIAN"
!   NUTRIENT_UNIT(20) = "N"
!   NUTRIENT_NAME(21) = "JDIAP"
!   NUTRIENT_UNIT(21) = "P"
!   NUTRIENT_NAME(22) = "JGREN"
!   NUTRIENT_UNIT(22) = "N"
!   NUTRIENT_NAME(23) = "JGREP"
!   NUTRIENT_UNIT(23) = "P"

#ifdef DEBUG
      DO J=1,NTT
       write(6,*) NUTRIENT_NAME(J),NUTRIENT_UNIT(J)
      END DO
#endif


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

!FV
      read(1,*) inputfile
      write(6,*) "inputfile",inputfile

      CLOSE(1)

!FV GoMDOM, Allocated vars and read inputfile
      call dummy() !Contains Allocate_Input, with ws(:)
      call Read_InputFile_GD(inputfile,idum,jdum)

       RETURN
END SUBROUTINE GET_PARAMETER_NEW


  SUBROUTINE BIO_PARAMETER_REPORT
    IMPLICIT NONE

  END SUBROUTINE BIO_PARAMETER_REPORT

END MODULE MOD_PARAMETER
