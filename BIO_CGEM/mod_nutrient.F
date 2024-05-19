MODULE MOD_NUTRIENT
!******************************************************************
!***        DEFINE PARAMETERS RELATED TO each state var         ***
!***   PARAMETERS ARE LISTED ALPHABETICALLY FOR EACH CATEGORY   ***
!******************************************************************
     USE MOD_1D
     IMPLICIT NONE
     SAVE

!Needed for mod_bio_3D
     REAL(SPP), ALLOCATABLE :: BIO_N(:,:)           ! NUTRIENT ABUNDANCE
     INTEGER                :: INN                  ! NUTRIENT INDEX USED IN VERTICAL MIXING

END MODULE MOD_NUTRIENT
