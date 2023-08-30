subroutine DETRITUS 
!==========================================================================
!  This routine computes source and sink terms of DETRITE (using DETRITE) in 1D              ===
!  without vertical mixing                                              ===
!==========================================================================
!
!  Define global data.
   USE MOD_1D
   USE MOD_DETRITUS
   
   IMPLICIT NONE
   INTEGER :: I,K,J,I1,I2, N_SINK
   REAL(SPP) :: DELTA_DMIN,W_D1,SOURCE
!=======================================================================
!  BEGINNING EXECUTABLE
!=======================================================================
!Definitions from other subroutines
!mod_bio_3D.F:        DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS
!*****************   SINKING   ************************
      DO J=1,NND
         DO K=1,KBVM1
           D_SINK(K,J)=0.0_SPP
           !FVD(K,J)=BIO_D(K,J)-D_0(J)
           ! DETRITUS AFTER THRESHOLD SUBSTRACTION
           DETRITE(K,J)=BIO_D(K,J)-D_0(J)
         END DO
      END DO

      DELTA_DMIN=1.0E+5_SPP    !CFL REQUIRES W*DTI<DELTA_DMIN, THE THINNEST LAYER
      DO K=1,KBVM1      
        DELTA_DMIN=MIN(DELTA_DMIN,DELTA_D(K))  
      END DO

      DO 10 J=1,NND
        W_D_TOTAL(J)=W_D(J)*T_STEP
        IF (W_D_TOTAL(J) > DELTA_DMIN) THEN
          N_SINK=INT(W_D_TOTAL(J)/DELTA_DMIN)+1
        ELSE
          N_SINK=1
        END IF
        W_D1=W_D_TOTAL(J)/REAL(N_SINK)

        DO I1=1,N_SINK
          D_SINK1(1)=-1.*W_D1*(DETRITE(1,J))/DELTA_D(1)    !SURFACE         
          DO K=2,KBV-2
            D_SINK1(K)=1.*W_D1*(DETRITE(K-1,J)-DETRITE(K,J))/DELTA_D(K)
          END DO

!WM/L3 modification, add flag for sinking
          if(isink_out.eq.1) then
           D_SINK1(KBVM1)=W_D1*(DETRITE(KBVM1-1,J)-DETRITE(KBVM1,J))/DELTA_D(KBVM1) !Sink out
          else
           D_SINK1(KBVM1)=W_D1*(DETRITE(KBVM1-1,J))/DELTA_D(KBVM1) !No sink out
          endif

          DO K=1,KBVM1
            DETRITE(K,J)=DETRITE(K,J)+D_SINK1(K)
            D_SINK(K,J)=D_SINK(K,J)+D_SINK1(K)
          END DO
        END DO
        DO K=1,KBVM1
          D_SINK(K,J)=D_SINK(K,J)/T_STEP
        END DO
 10    CONTINUE

!*****************   NUTRIENT SOURCES AND SINKS   *******************
      DO 20 J=1,NND
         DO K=1,KBVM1
           SOURCE=0.
           SOURCE=SOURCE                                   &
                     +D_SINK(K,J)                          ! SINKING
    
          BIO_D(K,J)=BIO_D(K,J)+SOURCE*T_STEP

         END DO !K=1,KBV
 20    CONTINUE
      RETURN
END SUBROUTINE DETRITUS 
