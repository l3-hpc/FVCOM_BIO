  subroutine call_iop_par(ff,nz,TC_8,PARsurf,dz,lat,lon,PARdepth)

  use cgem, only: nf,nospA,Kw,Kchla,Kspm,Kcdom,       &
   &              aw490,astar490,astarOMA,astarOMZ,astarOMR,astarOMBC,&
   &              iOM1CA,iOM1CZ,iOM1R,iOM1BC,iCDOM,CChla,iA,Qc,debug,CF_SPM
  use grid, only: iYrS,km
  use schism_glbl, only : rkind
  use date_time

  implicit none

  integer, intent(in)    :: TC_8         ! Model time in seconds from iYrS
  integer, intent(in)    :: nz           ! Number of layers
  real(rkind) :: lat,lon
  real(rkind), intent(in)    :: PARsurf      ! Irradiance just below sea surface
  real(rkind), intent(in), dimension(km,nf) :: ff      ! state variables
  real(rkind), intent(in), dimension(km)    :: dz      ! depth of cell
  real(rkind), intent(out), dimension(km)    :: PARdepth! PAR, visible irradiance at the middle 
                                                   !  of layer k (quanta/cm**2/sec)
  real(rkind), parameter :: C_cf  = 12.0E-3    ! C conversion factor (mmol-C/m3 to g-C/m3)
!----------------------------------------------------------------
! Calculate absorption (490 nm) components: seawater, chl, SPM from rivers, CDOM,
! detritus (dead cells), fecal pellets ...
  real(rkind), dimension(km) :: OM1A,OM1Z,OM1R,OM1BC,d_sfc,totChl
  real(rkind), dimension(km) :: Chla_tot, CDOM_tot, OM1A_tot, OM1Z_tot, OM1R_tot, OM1BC_tot
  real(rkind), dimension(km) :: CDOM !After converting ppb to a490 (m-1)
  real(rkind), dimension(km) :: Chla_mass, CDOM_mass, OM1A_mass, OM1Z_mass, OM1R_mass, OM1BC_mass
  real(rkind) :: Z ! Angle of the sun
  real(rkind) :: calc_solar_zenith !function
  real(rkind) :: a490_mid, aSw_mid, aChl490_mid, aCDOM490_mid, bbChl490_mid, bb490_mid
  real(rkind) :: aOM1A490_mid, aOM1Z490_mid, aOM1R490_mid, aOM1BC490_mid
! Time variables  
  real(rkind), parameter :: OneD60     = 1.0/60.0  ! Convert 1/min to 1/sec
  real(rkind), parameter :: OneD3600     = 1.0/3600.0  ! Convert 1/min to 1/sec
  real(rkind)            :: rhr          ! Decimal hour of day
  integer         :: iYr,iMon,iDay,iHr,iMin,iSec !Time variables
  integer         :: jday     ! Holds Julian Day
  integer         :: k,isp

 !---------------------------------------------------------
 ! -- Convert units for light model 
 !    C_cf == conversion factor (mmol-C/m3 to g-C/m3) 
 ! Organic Matter from dead phytoplankton (mmol/m3) 
   OM1A(1:nz) = ff(1:nz,iOM1CA) * C_cf
 ! Organic Matter from fecal pellets      (mmol/m3)
   OM1Z(1:nz) = ff(1:nz,iOM1CZ) * C_cf
 ! Suspended Particulate Matter (SPM)    (mmol/m3) 
 ! CF_SPM is fraction of Organic Matter in SPM originating from the rivers, default=0.018
   OM1R(1:nz) = ff(1:nz,iOM1R) * C_cf / CF_SPM 
 ! Organic Matter from boundary conditions(mmol/m3) 
   OM1BC(1:nz)  = ff(1:nz,iOM1BC) * C_cf
!----------------------------------------------------------------
  !Calculate Yr/Mon/Day from Model time in seconds        
  call date_timestamp(iYrS,TC_8,iYr,iMon,iDay,iHr,iMin,iSec)
  !Hours in day
  rhr = real(iHr,8) + real(iMin,8)*OneD60 + real(iSec,8)*OneD3600
  !Day of the year
  jday = JDAY_IN_YEAR(iYr, iMon, iDay)
  !Solar zenith
  Z =  calc_solar_zenith(lat,lon,rhr,jday) !in rad

! First, convert CDOM(ppb) into CDOM, a490 (m-1)
! Once the CDOM (QSE ppb) is in the model domain, we advect and mix using the same 
! routine as for other dissolved constituents. However, to use the CDOM in the light 
! attenuation models, we need to calculate a490 (Penta et al. 2008). 
! 1) convert CDOM(QSE ppb) back to a312: a312 = (CDOM(QSE ppb)-0.538)/2.933 
! 2) convert a312 to a490: a490 = a312*exp(-0.016*(490-312)), where here S = 0.016 
! (mean value from D'Sa and DiMarco (2008)
   do k=1,nz
      CDOM(k) = (ff(k,iCDOM) - 0.538)/2.933 !ppb to a312
      CDOM(k) = CDOM(k) * exp(-0.016*(490.-312.))
      CDOM(k) = DMAX1(CDOM(k),0.0)
   enddo 

   do k = 1, nz 
    totChl(k) = 0.0
    do isp = 1, nospA
      totChl(k) =  totChl(k) + ff(iA(isp),k) * Qc(isp) * 12. * (1./CChla(isp))
    enddo ! isp = 1, nospA
   enddo ! k = 1, km

!Mass in each cell at layer k (area of volume part cancels out)
!The unit is (g C)
      Chla_mass(1:nz) = totChl(1:nz)*dz(1:nz)
      CDOM_mass(1:nz) = CDOM(1:nz)*dz(1:nz)
      OM1A_mass(1:nz) = OM1A(1:nz)*dz(1:nz)
      OM1Z_mass(1:nz) = OM1Z(1:nz)*dz(1:nz)
      OM1R_mass(1:nz) = OM1R(1:nz)*dz(1:nz)
      OM1BC_mass(1:nz) = OM1BC(1:nz)*dz(1:nz)

!Mass from surface to center of cell at layer k
!Is the sum of the mass of all previous k layers plus 
!half of the current k layer 
!Concentration is mass/dz 
      Chla_tot(1) = 0.5*Chla_mass(1)
      CDOM_tot(1) = 0.5*CDOM_mass(1)
      OM1A_tot(1) = 0.5*OM1A_mass(1)
      OM1Z_tot(1) = 0.5*OM1Z_mass(1)
      OM1R_tot(1) = 0.5*OM1R_mass(1)
      OM1BC_tot(1) = 0.5*OM1BC_mass(1)
      d_sfc(1) = 0.5*dz(1)

      do k=2,nz
        Chla_tot(k)  = 0.5*Chla_mass(k) + SUM(Chla_mass(1:k-1))
        CDOM_tot(k)  = 0.5*CDOM_mass(k) + SUM(CDOM_mass(1:k-1))
        OM1A_tot(k)  = 0.5*OM1A_mass(k) + SUM(OM1A_mass(1:k-1))
        OM1Z_tot(k)  = 0.5*OM1Z_mass(k) + SUM(OM1Z_mass(1:k-1))
        OM1R_tot(k)  = 0.5*OM1R_mass(k) + SUM(OM1R_mass(1:k-1))
        OM1BC_tot(k) = 0.5*OM1BC_mass(k)+ SUM(OM1BC_mass(1:k-1))
        d_sfc(k)     = 0.5*dz(k)+ SUM(dz(1:k-1))
      enddo

     if(debug.eq.4) then
        write(6,*) "solar zenith",Z
        write(6,*) "Chla",Chla_tot
        write(6,*) "CDOM",CDOM_tot
        write(6,*) "OM1A",OM1A_tot
        write(6,*) "OM1Z",OM1Z_tot
        write(6,*) "OM1R",OM1R_tot
        write(6,*) "OM1BC",OM1BC_tot
        write(6,*) "d_sfc",d_sfc
     endif

   do k=1,nz
      !Calculate absorption coefficients:
      aSw_mid       = aw490                          !Sea water absorption at mid cell
      aChl490_mid   = astar490 * Chla_tot(k) / d_sfc(k)   !Chla absorption at mid cell
      aCDOM490_mid  = CDOM_tot(k) / d_sfc(k)              !CDOM absorption at mid cell
      aOM1A490_mid  = astarOMA * OM1A_tot(k) / d_sfc(k)   !OM1A absorption at mid cell
      aOM1Z490_mid  = astarOMZ * OM1Z_tot(k) / d_sfc(k)   !OM1Z absorption at mid cell
      aOM1R490_mid  = astarOMR * OM1R_tot(k) / d_sfc(k)   !OM1R absorption at mid cell
      aOM1BC490_mid = astarOMBC * OM1BC_tot(k) / d_sfc(k) !OM1BC absorption at mid cell
      a490_mid = aSw_mid + aChl490_mid + aCDOM490_mid + aOM1A490_mid + aOM1Z490_mid + aOM1R490_mid + aOM1BC490_mid
      !Calculate backscattering coefficients:
      bbChl490_mid = 0.015 * (0.3*((Chla_tot(k) / d_sfc(k))**0.62)*(550./490.)) !Chla backscatter at mid cell
      bb490_mid    = bbChl490_mid !Only Chla backscatters for now
      if(debug.eq.4) write(6,*) "k,a490, bb490_mid",k,a490_mid, bb490_mid
      call IOP_PARattenuation(a490_mid, bb490_mid, PARsurf, Z, d_sfc(k), PARdepth(k)) 
   enddo

   end subroutine Call_IOP_PAR
