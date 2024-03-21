Module MASS_BALANCE_GD 

IMPLICIT NONE

SAVE

REAL,ALLOCATABLE :: SUM_DENITR(:,:,:)    ! N loss to denitrification 
REAL,ALLOCATABLE :: SUM_DENITR_C(:,:,:)  ! C loss to denitrification 
REAL,ALLOCATABLE :: SUM_DOCPRD(:,:,:)    ! C loss to predation 
REAL,ALLOCATABLE :: SUM_DOCMET(:,:,:)    ! C loss to metabolism 
REAL,ALLOCATABLE :: SUM_DOCZOO(:,:,:)    ! C loss to zooplankton mortality 

contains

Subroutine MASS_BALANCE_GD_allocate

USE Model_dim

IMPLICIT NONE

SAVE

ALLOCATE(SUM_DENITR(im,jm,km))    ! N loss to denitrification 
ALLOCATE(SUM_DENITR_C(im,jm,km))  ! C loss to denitrification 
ALLOCATE(SUM_DOCPRD(im,jm,km))    ! C loss to predation 
ALLOCATE(SUM_DOCMET(im,jm,km))    ! C loss to metabolism 
ALLOCATE(SUM_DOCZOO(im,jm,km))    ! C loss to zooplankton mortality 

SUM_DENITR = 0.
SUM_DENITR_C = 0.
SUM_DOCPRD = 0.
SUM_DOCMET = 0.
SUM_DOCZOO = 0.

return

END Subroutine MASS_BALANCE_GD_allocate

END Module MASS_BALANCE_GD
