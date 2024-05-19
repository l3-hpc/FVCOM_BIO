!----------------------------------------------------------------------
      Subroutine getSolar( iYrS, TC_8, lon, lat, Rad) 
!----------------------------------------------------------------------
!     Written by  ::  D.S.Ko/NRL
!
!     Calculate visible solar radiation irradiance
!                        using just the solar zenith and ignore
!                        the effects of clouds and atmospheric
!                        absorption and scattering.
! ---------------------------------------------------------------------
     USE DATE_TIME
     use schism_glbl, only : rkind
!----------------------
! Interface variables:
!----------------------
      real(rkind)   , intent(in)  :: lon  ! longitude (deg E) at center of cell 
      real(rkind)   , intent(in)  :: lat  ! latitude (deg N) at center of cell 
      integer, intent(in)  :: iYrS
      integer, intent(in) :: TC_8 ! Current time in seconds since Model_dim::iYrS
      real(rkind)   , intent(out) :: Rad ! Solar Radiation
      integer  :: iYr      ! Year that Time_8 corresponds to
      integer  :: iMon     ! Month that Time_8 corresponds to
      integer  :: iDay     ! Day that Time_8 corresponds to
      integer  :: iHr      ! Hour that Time_8 corresponds to
      integer  :: iMin     ! Minute that Time_8 corresponds to
      integer  :: iSec     ! Second that Time_8 corresponds to
 
!-----------------------
! Local variables 
!-----------------------
      real(rkind), parameter :: cv        = 2.77e14 ! multiplicative factor used
                                                 ! to convert from watts/m2 
                                                 ! to photons/cm2/sec
                                                 ! Morel and Smith (1974)
      real(rkind), parameter :: OneD60   =  1./60.d0
      real(rkind), parameter :: OneD3600 =  1./3600.d0             
      integer  :: jday
      real(rkind)                :: rhr         ! decimal hr in the Julian Day 
      real(rkind)                :: Z   ! solar zenith angle 
      real(rkind) calc_solar_zenith
      real(rkind)                :: solconst      

         ! Note in next line that 1200 is the average clear-sky solar constant
         ! in watts/m^2
         solconst = 1200.00 * cv  ! in photons/cm2/sec
!-----------------------------------------------------------------------

         !Calculate Yr/Mon/Day from Model time in seconds	 
         call date_timestamp(iYrS,TC_8,iYr,iMon,iDay,iHr,iMin,iSec)

         !Hours in day
         rhr = real(iHr,8) + real(iMin,8)*OneD60 + real(iSec,8)*OneD3600

         !Day of the year
         jday = JDAY_IN_YEAR(iYr, iMon, iDay)
         Z =  calc_solar_zenith(lat,lon,rhr,jday) !in rad
         Rad = solconst * DMAX1( COS(Z), 0.d0)    ! COS(suna)<= 0 means night 
         !write(6,*) "Z,Rad",Z,Rad,iYr, iMon,iDay,lat,lon,rhr,jday

        RETURN
           
      END Subroutine getSolar
