subroutine surface_flux(ff_new,dT,dz,T,S,Wind,pH)

use cgem, only:Which_fluxes,iO2surf,iDICsurf,iO2,iDIC,iALK,iSi,iPO4, &
            & iA,nospA,nf,fmin,pCO2,debug
use cgem_utils 
use schism_glbl, only : rkind

implicit none

real(rkind), intent(in) :: dT,dz,T,S,Wind
real, intent(in) :: pH
real(rkind) :: pHrk
real(rkind), intent(inout) :: ff_new(nf)
real(rkind) :: O2_sfc, Sc, Op_umole, rhow, Op, OsDOp
real(rkind) :: Vtrans, alpha_O2, O2_atF,DIC_sfc, CO2_atF
real(rkind) :: SDay = 86400.d0
!------------------------------------------------------------------
if(debug.eq.3) write(6,*) "PH IN - SURFACE FLUX",pH

if(Which_fluxes(iO2surf).eq.1) then
!--------------------------------------------------------------
! Calc  O2_atF, the sea surface vertical flux of O2
!--------------------------------------------------------------
   O2_sfc   = ff_new(iO2) ! O2 (mmol-O2/m3) in sfc layer, k=1

   Sc       = SchmidtNumber(S,T,0)  ! Schmidt number,
                                                          !   0 (zero)
                                                          !   for O2

   Op_umole = o2sat(S,T)     ! O2 saturation,
                                                   !    (umol-O2/kg)

   rhow     = sw_dens0(S,T)  ! water density [kg/m3]

   Op       = rhow * Op_umole * 1.0E-3 ! O2 saturation,
                                                   !    (mmol-O2/m3)
   OsDOp    = O2_sfc/Op

!--------------------------------------------------------------
!  Vtrans below is the O2 transfer vel (m/s)
!
!  Vtrans   = (5.9*(kw)*(OsDOp*OsDOp))*(Sc)**X
!    where kw and Sc are dependent on Wind Speed.
!  Values kw and X are from Liss and Merlivat, 1986.
!  Factor of OsDOp**2 is from Justic, et. al 2002 
!   for when saturation levels are above 125%.
!---------------------------------------------------------------
 if(Wind.lt.3.6) then
   Vtrans        = DMAX1((5.9 * (0.17*Wind)         &
   &              * Sc**(-2./3.) / SDay), 0.)
 else if(Wind.le.13.) then
   Vtrans        = DMAX1((5.9 *(2.85*Wind - 9.65 )    &
   &              / SQRT(Sc) / SDay), 0.)
 else
   Vtrans        = DMAX1((5.9 *(5.9*Wind - 49.3 )    &
   &              / SQRT(Sc) / SDay), 0.)
 endif
 if(OsDOp.gt.1.25) Vtrans = Vtrans * (OsDOp*OsDOp)

   alpha_O2       = 1.025

   O2_atF         = Vtrans*(O2_sfc - alpha_O2*Op)
                                           ! flux of O2 thru
                                           ! the
                                           ! sea sfc
                                           ! ((mmol-O2/m2/sec)
                                           ! negative means
                                           ! into
   ff_new(iO2) = DMAX1(ff_new(iO2) - O2_atF/dz*dT,fmin(iO2))

endif

if(Which_fluxes(iDICsurf).eq.1) then
!--------------------------------------------------------------
! Calc  SFLUX_CO2, the sea surface vertical flux of CO2
!--------------------------------------------------------------
               DIC_sfc = ff_new(iDIC) ! Dissolved Inorganic Carbon
                                         !    (mmol m-3) in sfc layer, k=1

             !----------------------------------------------------------
             ! Units of gas_exchange are mmol CO2 m-2 s-1 
             !----------------------------------------------------------
!use mocsy instead but calculate to compare
             pHrk = real(pH,rkind)
             CO2_atF = gas_exchange(T,S,DIC_sfc,dz,pHrk,pCO2)
             ff_new(iDIC) = DMAX1(ff_new(iDIC) - CO2_atF/dz*dT,fmin(iDIC))
        if(debug.eq.3) write(6,*) "DIC_sfc,CO2_atF,pH,pHrk,pCO2,T,S,dz",DIC_sfc,CO2_atF,pH,pHrk,pCO2,T,S,dz

else

        write(6,*) "Invalid DIC surface flux",Which_fluxes(iDICsurf)
        stop
endif

RETURN

end subroutine surface_flux 
