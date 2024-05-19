!Dummy subroutine, defining rkind(used in schism)
MODULE schism_glbl
  USE mod_1d, only: SPP
  implicit none
  public 

  ! Get real datatype from FVCOM
  integer,parameter :: rkind = SPP

END MODULE schism_glbl
