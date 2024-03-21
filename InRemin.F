Module InRemin 
! =========================================================
! Fluxes for Instant Remineralization 
! =========================================================
       IMPLICIT NONE
   
       SAVE

       REAL,ALLOCATABLE :: SED_CARBON_RATE(:,:) 
       REAL,ALLOCATABLE :: TSOD(:,:) 
       REAL,ALLOCATABLE :: SED_NO3_RATE(:,:) 
       REAL,ALLOCATABLE :: SED_NH3_RATE(:,:) 

contains

SUBROUTINE InRemin_allocate 
! =========================================================
! Fluxes for Instant Remineralization 
! =========================================================
       USE Model_dim

       IMPLICIT NONE

       ALLOCATE(SED_CARBON_RATE(im,jm)) 
       ALLOCATE(TSOD(im,jm)) 
       ALLOCATE(SED_NO3_RATE(im,jm)) 
       ALLOCATE(SED_NH3_RATE(im,jm)) 

       SED_CARBON_RATE=0.
       TSOD=0.
       SED_NO3_RATE=0.
       SED_NH3_RATE=0.

return

END SUBROUTINE InRemin_allocate

END Module InRemin
 
