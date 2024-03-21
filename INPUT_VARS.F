Module INPUT_VARS

IMPLICIT NONE

SAVE

!Calculated:
integer(kind=8) :: START_SECONDS ! Seconds from iYr0 to start of run.
integer(kind=8) :: END_SECONDS   ! Seconds from iYr0 to end of run.
integer :: nstep  !number of timesteps in a run
integer :: iout   !output time-interval in timesteps

!--Code Run Identifier-----------
character(150) code_ID
!--Run Specifics---------------
integer iYrS,iMonS,iDayS,iHrS,iMinS,iSecS
integer iYrE,iMonE,iDayE,iHrE,iMinE,iSecE
integer dT, dT_out, dT_restart
integer icent,jcent

!--Switches ---------
integer Read_Solar,Read_Wind,Read_T,Read_Sal 
integer InitializeHow
integer Which_hydro
integer Which_Output
integer Which_VMix
integer Which_Adv
!----Sinking Terms----------------------------------
real,allocatable :: ws(:)

contains

Subroutine Allocate_Input

USE Model_dim, ONLY : nf

IMPLICIT NONE

integer ierr

!----Sinking Terms----------------------------------
ALLOCATE( ws(nf),stat=ierr ) 
!FV if(ierr.ne.0) call error("ws",ierr)
if(ierr.ne.0) write(6,*) "Error when allocating ws"


ws= 0.  !Default is no sinking

#ifdef DEBUG
write(6,*) "---INPUT_VARS---"
write(6,*) " general allocate, just ws, default is no sinking"
write(6,*)
#endif

return

END Subroutine Allocate_Input

END MODULE INPUT_VARS

