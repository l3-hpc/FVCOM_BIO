!=======================================================================
!  CGEM-FVCOM bio-geo-chem subroutine
!=======================================================================
SUBROUTINE bio_geo_chem(B_N,T_B,S_B,IRRAD,Wind)
!------------------------------------------------------------------------------
!Input: bio variables(column), temperature(column),
!           surface irradiance(cell), and windspeed(cell)
!Output: new bio variables
!mod_main.F:    INTEGER(itime) :: IINT     !!INTERNAL TIME STEP ITERATION NUMBER (ISTART => IEND)
!USE mod_main, only: IINT
!FVCOM  Define global data.
USE mod_1d, only: T_STEP,DELTA_D,DEPTH_Z,ITS,INEA,KBV
use schism_glbl, only: rkind
use cgem, only: ws,skipcgem,checkwindrad,sinkwcgem,nf,debug,Which_rad,Which_wind,&
                  PARfac,cgemcoords,cgemlat,cgemlon
use grid, only: iYrS,START_SECONDS,km

!Definitions from other subroutines
!mod_bio_3D.F:        DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS
!
IMPLICIT NONE

!GOMDOM
!INTEGER*8, INTENT(IN) :: TC_8
REAL(rkind), INTENT(IN) :: IRRAD,Wind
REAL(rkind), INTENT(IN) :: T_B(KBV),S_B(KBV)
REAL(rkind), INTENT(INOUT) :: B_N(KBV,nf)
real(rkind) :: ff(nf),ff_new(nf),ff_in(km,nf),ff_out(km,nf),dz(km),PAR(km),d_sfc(km)
real(rkind) :: cgem_wind,cgem_dz
integer k, i, myrank,nz,cgem_inea
logical :: is_surface,is_bottom
integer :: TC_8
real :: x
real(rkind) :: Rad
real(rkind), parameter :: cv        = 2.77e14_rkind ! multiplicative factor used
                                             ! to convert from watts/m2 
                                             ! to photons/cm2/sec
                                             ! Morel and Smith (1974)
i = 1
myrank = 1
nz=km

  !Time in integer seconds since start of run
  TC_8 = START_SECONDS + ITS*int(T_STEP)

    if(Which_rad.eq.1) then
      !Convert Rad from W/m2 to photons/cm2/sec
      Rad = IRRAD*cv
      !if(debug.eq.2.and.ITS.eq.1) write(6,*) "RAD,IRRAD,Wind,INEA",Rad,IRRAD,Wind,INEA
      !setting to nonzero for testing
      !Rad = cv
      !f(debug.eq.2) write(6,*) "For testing,RAD",Rad,IRRAD,Wind
    else
      call getSolar( iYrS, TC_8, cgemlon, cgemlat, Rad)
      ! Rad is just above sea surface.
      ! Rad was short wave generated by NRL, multiplied by SWtoPAR: ratio of PAR to
      ! shortwave radiation (hardcoded 4/30/14 to 0.43).
      ! Hardcoded to 0.47 on 2/11/16, Re: Tsubo and Walker, 2005
      ! PARfac is a multiplication factor for testing
      Rad = (0.47*Rad) * PARfac
    endif

    if(Which_wind.eq.1) then
      cgem_wind=Wind
    else
      !Set as constant if no wind fluxes available
      cgem_wind = 5.d0
    endif
!    if(debug.eq.2) write(6,*) "cgem_wind",cgem_wind

!L3 Check these!!!
!         DELTA_D(K)=DZ(I,K)*D(I)                   !LAYER THICKNESS
!         DELTA_Z(K)=DZZ(I,K)*D(I)                  !DISTANCE BETWEEN LAYERS
!         DEPTH_Z(K)=Z(I,K)*D(I)                    !LAYER CENTER DEPTH
ff_in = B_N(1:km,:)
d_sfc = -1.*DEPTH_Z + 0.5*DELTA_D
#if defined(DEBUG)
write(6,*) "d_sfc", d_sfc
write(6,*) "delta_d",DELTA_D
write(6,*) "delta_z",DEPTH_Z
write(6,*) "TC_8",TC_8
write(6,*) "Rad,lat,lon",Rad,cgemlat,cgemlon
#endif

call call_iop_par(ff_in,nz,TC_8,Rad,DELTA_D,cgemlat,cgemlon,PAR)
#if defined(DEBUG)
write(6,*) "Rad,PAR",Rad,PAR
#endif

     do k=1,km
       ff = ff_in(k,1:nf)
       ff_new = 0
       is_surface = .FALSE.
       is_bottom = .FALSE.
       !check these...does FVCOM always have equal layers?
       if(k.eq.1) is_surface = .TRUE.
       if(k.eq.nz) is_bottom = .TRUE.
       !write(6,*) "k,ff",k,ff
       call cgem_step(ff,ff_new,T_STEP, S_B(k), T_B(k), PAR(k), cgem_wind, cgemlat, DELTA_D(k), d_sfc(k), is_surface, is_bottom, Rad,INEA) !write(6,*) "ff_new",ff_new
       do i=1,nf
         if(ws(i).gt.tiny(x)) then
          !sink
          !call cgem_sink(ff(i),ff_new(i),ws(i),km,DELTA_D,T_STEP,i,ITS,myrank)
         endif
       enddo
       B_N(k,:) = ff_new
     enddo

     !write(6,*) "B_N",B_N

     !stop

!-----------------------------------------------------------------------------
RETURN

END SUBROUTINE bio_geo_chem
