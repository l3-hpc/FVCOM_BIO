SUBROUTINE DISSOLVED_OXYGEN(BIO_N,T_BIO)
!------------------------------------------------------------------------------
!-  Description : This subroutine performs kinetics calculations of 
!-                dissolved oxygen. 
!-  Created: 09/03/2009  James Pauer
!------------------------------------------------------------------------------
USE mod_1D, ONLY: SPP,KBV,KBVM1,NTT
USE STATES
USE EUT
USE FLAGS, ONLY: DO_DO2

IMPLICIT NONE

!
!------------------------------------------------------------------------------
!
REAL(SPP), INTENT(IN) :: BIO_N(KBV,NTT),T_BIO(KBV)
INTEGER :: K
REAL(SPP) TEMP

REAL(SPP) :: ALG               ! Phytoplankton concentration (dia+gre)
REAL(SPP) :: DOPRODD
REAL(SPP) :: DOMETD
REAL(SPP) :: DOPRODG
REAL(SPP) :: DOMETG
REAL(SPP) :: DOPREDD           ! DO consumption due to diatoms predation
REAL(SPP) :: DOPREDG           ! DO consumption due to greens predation
REAL(SPP) :: DOZOO             ! DO consumption due to zooplankton mortality
REAL(SPP) :: DOMNLDOC
REAL(SPP) :: DOCOD
REAL(SPP) :: FTN
REAL(SPP) :: KDOC              ! Mineralization rate of DOC [sec^-1]
REAL(SPP) :: NT
REAL(SPP) :: ZDEATH            ! Zooplankton mortality rate


!------------------------------------------------------------------------------
! Skip if not performing dissolved oxygen kinetics.
!------------------------------------------------------------------------------
IF (.NOT. DO_DO2) RETURN

DO K=1,KBVM1

TEMP=T_BIO(K)

!------------------------------------------------------------------------------
!- Nitrification based on LM3-Eutro/CE-Qual-ICM models
!------------------------------------------------------------------------------
   IF (TEMP < TMNT) THEN
       FTN = EXP(-KTNT1 * (TEMP - TMNT)**2)
   ELSE
       FTN = EXP(-KTNT2 * (TEMP - TMNT)**2)
   ENDIF


   IF (KHNNT > 0.0) THEN
       NT = FTN * NTM *                                    &
            ( BIO_N(K,JNH4) / (KHNNT  + BIO_N(K,JNH4)) ) *   &
            ( BIO_N(K,JDO2) / (KHDONT + BIO_N(K,JDO2)) )
   ELSEIF (BIO_N(K,JNH4) > 0.0) THEN
       NT = FTN * NTM * BIO_N(K,JDO2)/(KHDONT + BIO_N(K,JDO2))
   ELSE
       NT = 0.0
   ENDIF

   DTM(K,JDO2) = DTM(K,JDO2) - RNTO * NT 


!------------------------------------------------------------------------------
!- Oxygen sources and sinks from diatoms and greens. We
!- assume a threshold oxygen value below which all carbon 
!- is excreted as DOC, thus no DO consumed and CO2 generated
!- however, changes are necessary in other subroutines to 
!- change DOC fraction, FCDD when DO is below threshold value.
!- Also a oxygen:carbon conversion factor must be defined and
!- a value assigned. 
!------------------------------------------------------------------------------

    DOPRODD = PD(K) * BIO_N(K,JDIA) * RCDO
    DOPRODG = PG(K) * BIO_N(K,JGRE) * RCDO

    IF (BIO_N(K,JDO2) > (2.0 * 1.0E-03)) THEN
        DOMETD = (1.0 - (FCDD + FCLD + FCRD)) * BMD(K) * BIO_N(K,JDIA) * RCDO
        DOMETG = (1.0 - (FCDG + FCLG + FCRG)) * BMG(K) * BIO_N(K,JGRE) * RCDO
        DOPREDD = (1.0 - GREFF) * (1. - (FCDP + FCLP + FCRP)) * &
                 &  PRD(K) * RCDO
        DOPREDG = (1.0 - GREFF) * (1. - (FCDP + FCLP + FCRP)) * &
                 &  PRG(K) * RCDO
        ZDEATH = ZDTH * BIO_N(K,JZOO) * BIO_N(K,JZOO)
        DOZOO   =  (1.0 - (FCDZ + FCLZ + FCRZ)) * ZDEATH * RCDO
    ELSE
        DOMETD  = 0.0
        DOMETG  = 0.0
        DOPREDD = 0.0
        DOPREDG = 0.0
        DOZOO   = 0.0
    ENDIF


    DTM(K,JDO2) = DTM(K,JDO2) + DOPRODD - DOMETD + DOPRODG - DOMETG -  & 
                   &  DOPREDD - DOPREDG - DOZOO


!------------------------------------------------------------------------------
!- DOC and COD subroutine
!------------------------------------------------------------------------------

    ALG      = BIO_N(K,JGRE) + BIO_N(K,JDIA)
    KDOC     = (KDC + (KDCALG * ALG)) * FTMNL(K) 
    DOMNLDOC = ( BIO_N(K,JDO2) / (KHODOC + BIO_N(K,JDO2)) )      & 
             *  KDOC * BIO_N(K,JDOC) * RCDO
    
    !(WOD undefined)  IF (DO_WOD) THEN
    !    DOCOD  = ( BIO_N(K,JDO2) / (KHOCOD + BIO_N(K,JDO2)) )      &
    !           *  KCOD * WOD(K) * RCDO 
    !ELSE
        DOCOD  = 0.0
    !ENDIF

    DTM(K,JDO2) = DTM(K,JDO2) -(DOMNLDOC + DOCOD) 

    !DOMETD_ARR(K) = DOMETD_ARR(K) + DOMETD * Vol * real(dT,4)
    !DOMETG_ARR(K) = DOMETG_ARR(K) + DOMETG * Vol * real(dT,4)
    !DOPREDD_ARR(K) = DOPREDD_ARR(K) + DOPREDD * Vol * real(dT,4)
    !DOPREDG_ARR(K) = DOPREDG_ARR(K) + DOPREDG * Vol * real(dT,4)
    !DOZOO_ARR(K) = DOZOO_ARR(K) + DOZOO * Vol * real(dT,4)
    !DOMNLDOC_ARR(K) = DOMNLDOC_ARR(K) + DOMNLDOC * Vol * real(dT,4)
    !NITDO2(K) = NITDO2(K) + RNTO * NT * Vol * real(dT,4)

ENDDO

!------------------------------------------------------------------------------
END SUBROUTINE DISSOLVED_OXYGEN
