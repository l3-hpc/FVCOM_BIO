!----------------------------------------------------------------------
      FUNCTION calc_solar_zenith(lat, lon, rhr, julianDay )    &      
      &        RESULT(sunangle) 
      
!-----------------------------------------------------------------------
!     Original MATLAB code Written by: Brad Penta/NRL
!
!     Translated into FORTRAN and 
!     Revised by                     : Barry Herchenroder/EMVL, April
!     2010
!                                                               June
!                                                               2011
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Calculate the solar zenith angle (radians) given the latitude 
!  (lat, degrees) longitude (lon, degrees) of a point in the Gulf of 
!  Mexico (GOM) as well as the Julian day of the year (julianDay),
!  the hour in that day (ihr-integer)
!------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! Definitions
  !
  ! sunangle     -> angle of the sun in radians, computed from latitude,
  !                 longitude, Julian day, and time--This is the zenith 
  !
  ! Compute the sun zenith angle in radians for a given Julian day
  ! (julianDay),
  ! GMT (ihr;in hours), latitude (lat), and longitude (lon).
  !
  ! Sun angle of pi/2 radians (90 deg) = sun at or below horizon i.e.
  ! dark
  ! Sun angle of 0 radians (0 deg)     = sun directly overhead.
!------------------------------------------------------------------------------

  IMPLICIT NONE 
   
!-------------------------------------------------
! Declare variables in the interface
!------------------------------------------------- 
  REAL*8   , INTENT(IN)    :: lat          ! Latitude (deg); lat > 0 means
                                         !                 North of
                                         !                 Equator
  REAL*8   , INTENT(IN)    :: lon          ! Longitude (deg E, 0 <= lon < 360); 
  REAL*8   , INTENT(IN)    :: rhr          ! decimal Hour of Julian Day  
  INTEGER, INTENT(IN)    :: julianDay    ! Julian Day GMT
  
!-------------------------------------------------
! Local variables
!------------------------------------------------- 
  REAL*8   :: costmp
  REAL*8   :: sdec
  REAL*8   :: rad
  REAL*8   :: rlat 
  REAL*8   :: rsdec   
  REAL*8   :: rthez 
  REAL*8   :: sunangle  ! returned as RESULT 
  REAL*8   :: tc  
  REAL*8   :: thez
  REAL*8   :: xha


#if defined(DEBUG)
write(6,*) "In calc_solar_zenith"
#endif



!---------------------------------------------------
!  Compute solar declination angle
  rad = 180.0/ (4.0 * atan(1.0))             ! rad = radian to degrees 
                                             !       conversion factor
  thez = 360.0*((FLOAT(julianDay-1)) )/365. ! thez = theta zero orbital 
                                                 !        position
                                                 !        (degs)
  rthez = thez/rad;  ! in radians

!  sdec = solar declination angle in degrees
   sdec = 0.396372                                                    &
   &    -22.91327 *cos(rthez)     + 4.02543*sin(rthez)                & 
   &    - 0.387205*cos(2.0*rthez) + 0.051967*sin(2.0*rthez)           &
   &    - 0.154527*cos(3.0*rthez) + 0.084798*sin(3.0*rthez)           
   rsdec = sdec/rad; ! in radians

#if defined(DEBUG)
!write(6,*) "In calc_solar_zenith,rsdec,rad=",rsdec,rad
#endif

!  Time correction (tc) for solar hour angle
   tc = 0.004297                                                      &
   &  + 0.107029*cos(rthez)     - 1.837877*sin(rthez)                 &
   &  - 0.837378*cos(2.0*rthez) - 2.342824*sin(2.0*rthez)

!  xha = solar hour angle in degrees
   xha = (rhr-12.0 )*15.0 + lon + tc    
     
   IF (xha > 180.0) THEN
      xha = xha - 360.0
   ENDIF 
   IF (xha > -180.0) THEN
      xha = xha + 360.0
   ENDIF 

!  Sun zenith angle sunangle. 
! Cosine(sunangle) == costmp
! rsdec            == sun declination angle
! xha              == local hour angle  

   rlat = lat/rad
   costmp = sin(rlat)*sin(rsdec)                                      & 
   &      + cos(rlat)*cos(rsdec) * cos(xha/rad)

! The next block if is for failsafe purposes
   IF (abs(costmp) > 1.0) THEN
       IF (costmp > 0.0) THEN
           costmp = 1.0
       ELSE
           costmp = -1.0
       ENDIF 
   ENDIF 

   sunangle = acos(costmp)  ! solar zenith angle in radians. Note that
                            ! 2.0*atan(1.0)= Pi/2 radians= 90 deg

   IF (sunangle > (2.0*atan(1.0))) THEN
       sunangle = 2.0*atan(1.0);
   ENDIF 

#if defined(DEBUG)
write(6,*) "End calc_solar_zenith"
#endif

 
   RETURN  
   END FUNCTION calc_solar_zenith   
