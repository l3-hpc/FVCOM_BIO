subroutine cgem_sink(ff,ff_new,ws,nz,dz,dT,i,istep,myrank)
  !This is called after cgem_step which returns ff
  use grid, only:km
  !use cgem, only:nf,adjust_ws,adjust_fac,debug
  use cgem, only:nf,debug
  use schism_glbl, only : rkind

  implicit none

  integer, intent(in) :: i,istep,myrank,nz
  real(rkind), intent(in) :: dT
  real(rkind), intent(in) :: ws 
  real(rkind), intent(in) :: dz(km)
  real(rkind), intent(in) :: ff(km)
  real(rkind), intent(out) :: ff_new(km)
  integer :: k, km1
  real(rkind) :: x
  real(rkind) :: cmin
  real(rkind) :: wsc
  real(rkind), dimension(km) :: mass_in,mass_out,d_mass

  km1 = nz-1
  ! !If sinking rate violates courant condition, change it
  !if(adjust_ws) then
  !  cmin = adjust_fac*MINVAL(dz)/dT
  !  wsc = DMIN1(cmin,ws)
  !else  !or not
    wsc = ws
  !endif

  !Nothing sinking in
  mass_in(1) = 0.d0
  !area is constant, cancels
  mass_out(1:km1) = ff(1:km1)*wsc*dz(1:km1)
  !Don't sink out
  mass_out(nz) = 0.d0
  !mass in should be mass out
  mass_in(2:nz) = mass_out(1:km1)
  d_mass = mass_in - mass_out
  ff_new(1:nz) = ff(1:nz) + d_mass(1:nz)/dz(1:nz)*dt

  if(debug.eq.1.and.i.eq.10) then
    write(6,*) "sinking"
    write(6,'(*(g0,:,", "))') istep,i,SUM(ff(1:nz)),SUM(ff_new(1:nz))
    do k=1,nz
      write(6,'(*(g0,:,", "))') ff(k),ff_new(k),mass_in(k),mass_out(k),d_mass(k),dz(k)
    enddo
 endif


  do k=1,nz !Check for zeros...if this is zero, fix something.
    if(ff_new(k).le.0) then
       write(6,'(*(g0,:,", "))') "k,fold,fnew,wsc,dz=",k,ff(k),ff_new(k),wsc,dz(k)
        ff_new(k) = 0. 
      stop 
     endif
  enddo

  return
end subroutine cgem_sink
