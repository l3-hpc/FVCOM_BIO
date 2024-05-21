subroutine grid_setup(nvrts)

  use grid
  use DATE_TIME 
  implicit none

  integer, intent(in) :: nvrts

  km=nvrts-1

  ! Compute starting time of run in seconds since Model_dim::iYrS:
  START_SECONDS = &
  TOTAL_SECONDS( iYrS, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )

#if defined(DEBUG)
write(6,*) "km",km
write(6,*) "START_SECONDS",START_SECONDS
write(6,*) "iYr,etc,",iYrS, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS
#endif


return
end subroutine

