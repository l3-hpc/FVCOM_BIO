subroutine grid_setup(nvrts)

  use grid
  use DATE_TIME 
  implicit none

  integer, intent(in) :: nvrts

  km=nvrts-1

  ! Compute starting time of run in seconds since Model_dim::iYrS:
  START_SECONDS = &
  TOTAL_SECONDS( iYrS, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )

return
end subroutine

