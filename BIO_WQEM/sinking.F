subroutine sinking(BIO_N) 
!==========================================================================
!  This routine computes source and sink terms 
!==========================================================================
!
!  Define global data.
   USE MOD_1D, only:SPP,NTT,KBV,KBVM1,DELTA_D,T_STEP,ws,BIO_0,isink_out
   
   IMPLICIT NONE
   INTEGER :: I,K,J,I1,I2, I_SINK
   REAL(SPP) :: DELTA_DMIN,W_N1,W_N_TOTAL(NTT)
   REAL(SPP),intent(inout) :: BIO_N(KBV,NTT)
   REAL(SPP) :: BIO_SINK(KBV,NTT),N_SINK(KBV,NTT),N_SINK1(KBV)

!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================
!Definitions from other subroutines
!mod_bio_3D.F:        DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS
!*****************   SINKING   ************************
      DO J=1,NTT
         DO K=1,KBVM1
           N_SINK(K,J)=0.0_SPP
           !FVD(K,J)=BIO_D(K,J)-D_0(J)
           ! DETRITUS AFTER THRESHOLD SUBSTRACTION
           BIO_SINK(K,J)=BIO_N(K,J)-BIO_0(J)
         END DO
      END DO

      DELTA_DMIN=1.0E+5_SPP    !CFL REQUIRES W*DTI<DELTA_DMIN, THE THINNEST LAYER
      DO K=1,KBVM1      
        DELTA_DMIN=MIN(DELTA_DMIN,DELTA_D(K))  
      END DO

      DO 10 J=1,NTT
        W_N_TOTAL(J)=ws(J)*T_STEP
        IF (W_N_TOTAL(J) > DELTA_DMIN) THEN
          I_SINK=INT(W_N_TOTAL(J)/DELTA_DMIN)+1
        ELSE
          I_SINK=1
        END IF
        W_N1=W_N_TOTAL(J)/REAL(I_SINK)

        DO I1=1,I_SINK
          N_SINK1(1)=-1.*W_N1*(BIO_SINK(1,J))/DELTA_D(1)    !SURFACE         
          DO K=2,KBV-2
            N_SINK1(K)=1.*W_N1*(BIO_SINK(K-1,J)-BIO_SINK(K,J))/DELTA_D(K)
          END DO

!WM/L3 modification, add flag for sinking
          if(isink_out.eq.1) then
           N_SINK1(KBVM1)=W_N1*(BIO_SINK(KBVM1-1,J)-BIO_SINK(KBVM1,J))/DELTA_D(KBVM1) !Sink out
          else
           N_SINK1(KBVM1)=W_N1*(BIO_SINK(KBVM1-1,J))/DELTA_D(KBVM1) !No sink out
          endif

          DO K=1,KBVM1
            BIO_SINK(K,J)=BIO_SINK(K,J)+N_SINK1(K)
            N_SINK(K,J)=N_SINK(K,J)+N_SINK1(K)
          END DO
        END DO

 10    CONTINUE

!*****************   NUTRIENT SOURCES AND SINKS   *******************
      DO 20 J=1,NTT
         DO K=1,KBVM1
            BIO_N(K,J)=BIO_N(K,J)+N_SINK(K,J) ! SINKING
         END DO !K=1,KBVM1
 20    CONTINUE
      RETURN
END SUBROUTINE sinking 
