!----------------------------------------------------------------------
!*****************************************************************************************
      subroutine IOP_PARattenuation(a490, bb490, PARsurf, sun_zenith, d_sfc, PARdepth)

      use schism_glbl, only : rkind
      implicit none
! Penta et al., 2008 PAR penetration model: based on the IOP (inherent optical properties)
! absorption and backscatter at 490nm and the zenith angle of the sun; modified from Lee 
! et al., 2005 PAR penetration model: based on satellite absorption and backscater at 490nm

!    Lee, Z., K. Du, R. Arnone, S. Liew, and B. Penta, .Penetration of solar radiation
!         in the upper ocean: a numerical model for oceanic and coastal waters,. 
!         J. Geophys. Res. 110, C09019 doi:10.1029/2004JC002780 (2005).
!    Penta, B., Z. Lee, R.M. Kudela, S.L. Palacios, D.J. Gray, J.K. Jolliff, and 
!         I.G. Shulman, .An underwater light attenuation scheme for marine ecosystem 
!         models., Optics Express, 16, 16582-16591 (2008).

! Absorption (a490) is the total absorption at 490 nm calculated from all of the 
! constituents of the model that absorb light (Seawater, chlorophyll, CDOM, detritus, etc.
! Backscatter (bb490) is similar 

! Use the average value from the sea-surface to the depth of the calculation 
! (these calculations moved to a new subroutine(s) and the a490 and bb490 will be passed 
! into this subroutine  

! PAR (photosynthetically active radiation) just below the sea surface is used as a 
! starting value to be multiplied by the attenuation factor computed in this subroutine. 
! The starting value does not affect any of these calculations

! We generally use a factor of 0.4815 to represent the loss of PAR passing across the
! Air/Sea interface - this is included already in the data 

! Define coefficients for light attenuation model       
      real(rkind) a490, alpha0, alpha1, alpha2, bb490, chi0, chi1, chi2, d_sfc
      real(rkind) k1, k2, kpar, PARdepth, PARsurf, sun_zenith, zeta0, zeta1
      real(rkind) zeta2 
     
! sun_zenith => Solar Zenith Angle in radians for time and location
! d_sfc => depth of center of cell from elevated sea surface        

! Set values of coefficients for light attenuation model 
! (Lee et al., 2005; Penta et al., 2008)       
      parameter(alpha0=0.090, alpha1=1.465, alpha2=-0.667, chi0=-0.057, chi1=0.482, & 
     &          chi2=4.221, zeta0=0.183, zeta1=0.702, zeta2=-2.567)
     
! Calculate the attenuation (Equation 2 Penta et al., 2008 errata)
      k1 = ((chi0 + chi1*(a490**0.5) + chi2*bb490) & 
     &     * (1 + alpha0 * sin(sun_zenith)))
     
      k2 = ((zeta0 + zeta1*a490 + zeta2*bb490) &
     &     * (alpha1 + alpha2 * cos(sun_zenith) ))
     
      kpar = k1 + (k2 / (1+(d_sfc))**0.5)
     
      PARdepth = PARsurf * exp(-(kpar*(d_sfc)))

      return
      end subroutine IOP_PARattenuation

!---------------------------------------------------------------------- 
