MODULE MOD_NUTRIENT
!******************************************************************
!***        DEFINE PARAMETERS RELATED TO each state var         ***
!***   PARAMETERS ARE LISTED ALPHABETICALLY FOR EACH CATEGORY   ***
!******************************************************************
     USE MOD_1D
     IMPLICIT NONE
     SAVE

!Needed for mod_bio_3D
     INTEGER                :: NNN                  ! TOTAL NUMBER OF BACTERIA
     INTEGER                :: INN                  ! NUTRIENT INDEX USED IN VERTICAL MIXING
     REAL(SPP), ALLOCATABLE :: BIO_N(:,:)           ! NUTRIENT ABUNDANCE

!Input parameters
     REAL(SPP), ALLOCATABLE :: N_0(:)                ! NUTRIENT THRESHOLD
     REAL(SPP), ALLOCATABLE :: W_N(:)                ! NUTRIENT SINKING VELOCITY
     INTEGER ISINK_OUT

!Used for calculations
     REAL(SPP), ALLOCATABLE :: NUTRIENT(:,:)          ! NUTRIENT CONCENTRATION 
     REAL(SPP), ALLOCATABLE :: N_SINK(:,:)           ! NUTRIENT SINKING FLUX
     REAL(SPP), ALLOCATABLE :: N_SINK1(:)            ! SINKING FLUX WITHIN 1 ITERATION
     REAL(SPP), ALLOCATABLE :: W_N_TOTAL(:)          ! SINKING DISTANCE WITIN 1 T_STEP

END MODULE MOD_NUTRIENT
