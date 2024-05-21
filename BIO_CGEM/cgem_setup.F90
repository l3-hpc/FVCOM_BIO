subroutine cgem_setup(tracers)

  use cgem

  integer, intent(in) :: tracers !value of ntrs(3)
  integer check_tracers

#if defined(DEBUG)
write(6,*) "Begin cgem_setup"
#endif

  call cgem_dim  !Read nospA and nospZ

#if defined(DEBUG)
write(6,*) "Begin cgem_init"
#endif

!SCHISM 
!We have to figure out what to do about tracers
!check that nospA and nospZ are such that tracers was set correctly in param.nml
  check_tracers = nospA*3+nospZ+25
#if defined(DEBUG)
write(6,*) "In cgem_setup,nospA,nospZ,tracers,check_tracers",nospA,nospZ,tracers,check_tracers
#endif


  if(tracers.ne.check_tracers) then
    write(6,*) "Schism tracers",tracers,"not equal cgem tracers",check_tracers
    write(6,*) "Either change cgem.nml or param.nml"
    stop 
  endif

  call cgem_allocate
  call cgem_read
  call cgem_init

#if defined(DEBUG)
write(6,*) "End cgem_setup,tracers,check_tracers",tracers,check_tracers
#endif


return
end subroutine cgem_setup

