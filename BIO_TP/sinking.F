subroutine sinking(BIO_N)
!==========================================================================
!  This routine computes source and sink terms of NUTRIENT (using NUTRIENT) in 1D              ===
!  without vertical mixing                                              ===
!==========================================================================
!
!  Define global data.
   USE MOD_1D, only: SPP,KBV,KBVM1,NTT,DELTA_D,T_STEP,W_N,N_0,ISINK_OUT
!  USE MOD_NUTRIENT
   
   IMPLICIT NONE
   REAL(SPP), intent(inout) :: BIO_N(KBV,NTT)
   REAL(SPP) :: W_N_TOTAL(NTT)
   REAL(SPP) :: NUTRIENT(KBV,NTT),N_SINK(KBV,NTT)
   REAL(SPP) :: N_SINK1(KBV)
   INTEGER :: I,K,J,I1,I2, I_SINK
   REAL(SPP) :: DELTA_DMIN,W_N1

!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================
!Definitions from other subroutines
!mod_bio_3D.F:        DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS
!*****************   SINKING   ************************
      DO J=1,NTT
         DO K=1,KBVM1
           N_SINK(K,J)=0.0_SPP
           ! NUTRIENT AFTER THRESHOLD SUBSTRACTION
           NUTRIENT(K,J)=BIO_N(K,J)-N_0(J)
         END DO
      END DO

      DELTA_DMIN=1.0E+5_SPP    !CFL REQUIRES W*DTI<DELTA_DMIN, THE THINNEST LAYER
      DO K=1,KBVM1      
        DELTA_DMIN=MIN(DELTA_DMIN,DELTA_D(K))  
      END DO

      DO 10 J=1,NTT
        W_N_TOTAL(J)=W_N(J)*T_STEP
        IF (W_N_TOTAL(J) > DELTA_DMIN) THEN
          I_SINK=INT(W_N_TOTAL(J)/DELTA_DMIN)+1
        ELSE
          I_SINK=1
        END IF
        W_N1=W_N_TOTAL(J)/REAL(I_SINK)

        DO I1=1,I_SINK
          N_SINK1(1)=-1.*W_N1*(NUTRIENT(1,J))/DELTA_D(1)    !SURFACE         
          DO K=2,KBV-2
            N_SINK1(K)=1.*W_N1*(NUTRIENT(K-1,J)-NUTRIENT(K,J))/DELTA_D(K)
          END DO

!WM/L3 modification, add flag for sinking
          if(isink_out(j).eq.1) then
           N_SINK1(KBVM1)=W_N1*(NUTRIENT(KBVM1-1,J)-NUTRIENT(KBVM1,J))/DELTA_D(KBVM1) !Sink out
          else
           N_SINK1(KBVM1)=W_N1*(NUTRIENT(KBVM1-1,J))/DELTA_D(KBVM1) !No sink out
          endif

          DO K=1,KBVM1
            NUTRIENT(K,J)=NUTRIENT(K,J)+N_SINK1(K)
            N_SINK(K,J)=N_SINK(K,J)+N_SINK1(K)
          END DO
        END DO !I1

 10    CONTINUE

!*****************   NUTRIENT SOURCES AND SINKS   *******************
      DO 20 J=1,NTT
         DO K=1,KBVM1
            BIO_N(K,J)=BIO_N(K,J)+N_SINK(K,J) ! SINKING
         END DO !K=1,KBVM1
 20    CONTINUE
      RETURN
END SUBROUTINE sinking 
