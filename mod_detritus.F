MODULE MOD_DETRITUS
!******************************************************************
!***        DEFINE PARAMETERS RELATED TO DETRITUS               ***
!***   PARAMETERS ARE LISTED ALPHABETICALLY FOR EACH CATEGORY   ***
!******************************************************************
     USE MOD_1D
     IMPLICIT NONE
     SAVE

!Needed for mod_bio_3D
     INTEGER                :: NND                  ! TOTAL NUMBER OF BACTERIA
     INTEGER                :: IND                  ! DETRITUS INDEX USED IN VERTICAL MIXING
     REAL(SPP), ALLOCATABLE :: BIO_D(:,:)           ! DETRITUS ABUNDANCE

!Input parameters
     REAL(SPP), ALLOCATABLE :: D_0(:)                ! DETRITUS THRESHOLD
     REAL(SPP), ALLOCATABLE :: W_D(:)                ! DETRITUS SINKING VELOCITY
     REAL(SPP)  N2CD                      ! Element ratio in detritus 
     INTEGER ISINK_OUT

!Used for calculations
     REAL(SPP), ALLOCATABLE :: DETRITE(:,:)          ! DETRITUS CONCENTRATION 
     REAL(SPP), ALLOCATABLE :: D_SINK(:,:)           ! DETRITUS SINKING FLUX
     REAL(SPP), ALLOCATABLE :: D_SINK1(:)            ! SINKING FLUX WITHIN 1 ITERATION
     REAL(SPP), ALLOCATABLE :: W_D_TOTAL(:)          ! SINKING DISTANCE WITIN 1 T_STEP

END MODULE MOD_DETRITUS
