module cgem

!CGEM STATE VARIABLES
use, intrinsic :: iso_fortran_env, only: stderr => error_unit
use schism_glbl, only : rkind

implicit none

save

!number of species
integer :: nospA
integer :: nospZ

!misc
logical     :: sinkwcgem,cgemcoords
integer     :: skipcgem,checkwindrad,debug
real(rkind) :: eps
!real(rkind) :: adjust_ws,adjust_fac
real(rkind) :: cgemlat,cgemlon
!Sinking
real(rkind), dimension(:), allocatable :: ws
real(rkind), dimension(:), allocatable :: fmin
real(rkind) :: x !for tiny

!mocsy/pH
real(rkind) :: pCO2

!Fixed Stoichiometry:
real(rkind) :: sx1R,sy1R,sx2R,sy2R
real(rkind) :: sx1BC,sy1BC,sx2BC,sy2BC

!Module Which_Flux
! =========================================================
! Define which fluxes shall be used
! =========================================================
integer, parameter :: iO2surf  = 1 !O2 surface flux
integer, parameter :: iDICsurf = 2 !DIC surface flux
integer, parameter :: iInRemin = 3 !Instant Remineralization in bottom layer
!integer, parameter :: iSOC     = 4 !Sediment Oxygen Consumption
!integer, parameter :: iMPB     = 5 !Microphytobethos
!integer, parameter :: iNutEx   = 6 !Sediment Nutrient Fluxes
!integer, parameter :: iCMAQ    = 7 !CMAQ surface deposition of NH4 and NO3
!integer, parameter :: iInRemin = 8 !Instant Remineralization in bottom layer
!integer, parameter :: iSDM     = 9 !Sediment Diagenesis Model

!Module CGEM_Flux
!---------------------------------------------------------      
!-A; Phytoplankton number density (cells/m3);
!---------------------------------------------------------  
      integer, dimension(:), allocatable :: iA(:)
!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
      integer, dimension(:), allocatable :: iQn(:)
!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
      integer, dimension(:), allocatable :: iQp(:)
!--------------------------------------------------------------------
!-Z: Zooplankton number density (individuals/m3);
!--------------------------------------------------------------------
      integer, dimension(:), allocatable :: iZ(:)
!-------------------------------
!-NO3; Nitrate (mmol-N/m3)
!-------------------------------
      integer :: iNO3
!--------------------------------      
!-NH4; Ammonium (mmol-N/m3)
!--------------------------------
      integer :: iNH4
!-------------------------------------------        
!-PO4: Phosphate (mmol-P/m3)
!--------------------------------------
      integer :: iPO4 
!---------------------------------------------------------
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) 
!---------------------------------------------------------
      integer :: iDIC 
!-------------------------------------------        
!-O2: Molecular Oxygen (mmol-O2/m3)
!------------------------------
      integer :: iO2 
!-------------------------------------------------------------
!-OM1CA: (mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton, Carbon content
!-------------------------------------------------------------
      integer :: iOM1CA
!-------------------------------------------------------------
!-OM1NA: (mmol-N/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton, Nitrogen content
!-------------------------------------------------------------
      integer :: iOM1NA
!-------------------------------------------------------------
!-OM1PA: (mmol-P/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton, Phosphorus content
!-------------------------------------------------------------
      integer :: iOM1PA
!-----------------------------------------------------------------
!-OM2CA: (mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton, Carbon content
!------------------------------------------------------------------
      integer :: iOM2CA
!-----------------------------------------------------------------
!-OM2NA: (mmol-N/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton, Nitrogen content
!------------------------------------------------------------------
      integer :: iOM2NA
!-----------------------------------------------------------------
!-OM2PA: (mmol-P/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton, Phosphorus content
!------------------------------------------------------------------
      integer :: iOM2PA
!-------------------------------------------------------------
!-OM1CZ:(mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets, Carbon content
!-------------------------------------------------------------
      integer :: iOM1CZ
!-------------------------------------------------------------
!-OM1NZ:(mmol-N/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets, Nitrogen content
!-------------------------------------------------------------
      integer :: iOM1NZ
!-------------------------------------------------------------
!-OM1PZ:(mmol-P/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets, Phosphorus content
!-------------------------------------------------------------
      integer :: iOM1PZ
!-------------------------------------------------        
!-OM2CZ:(mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets, Carbon content
!-----------------------------------------------
      integer :: iOM2CZ
!-------------------------------------------------        
!-OM2NZ:(mmol-N/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets, Nitrogen content
!-----------------------------------------------
      integer :: iOM2NZ
!-------------------------------------------------        
!-OM2PZ:(mmol-P/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets, Phosphorus content
!-----------------------------------------------
      integer :: iOM2PZ
!--------------------------------------------------------------------
!-OM1R: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter arising from river outflow
!--------------------------------------------------------------------
      integer :: iOM1R
!-------------------------------------------------      
!-OM2R: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter arising from river outflow
!--------------------------------------------------------------------
      integer :: iOM2R
!-------------------------------------------
!-CDOM: (ppb) 
!        -- Colored Dissolved Organic Matter
!-------------------------------------------
      integer :: iCDOM
!---------------------------------------------
!-Silica: (mmol-Si/m3) 
!        -- Silica
!-------------------------------------------
      integer :: iSi
!--------------------------------------------------------------------
!-OM1BC: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter in initial and boundary 
!            conditions 
!--------------------------------------------------------------------
      integer :: iOM1BC
!-------------------------------------------------
!-OM2BC: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter in initial and boundary
!            conditions
!--------------------------------------------------------------------
      integer :: iOM2BC
!-------------------------------------------
!-ALK:  (mmol-HCO3/m3)?
!        -- Alkalinity
!-------------------------------------------
      integer :: iALK
!-------------------------------------------
!-Tr: tracer
!-------------------------------------------
      integer :: iTr

!Total number of state variables
      integer :: nf

real(rkind) :: sinkTr
!----INPUT_VARS_CGEM
!--Switches in GEM---------
!integer Which_fluxes(8)
integer Which_fluxes(3)
integer Which_uptake
integer Which_quota
integer Which_photosynthesis
integer Which_growth
integer Which_temperature
integer Which_wind
integer Which_rad
!--Temperature
real(rkind), allocatable :: KTg1(:)
real(rkind), allocatable :: KTg2(:)
real(rkind), allocatable :: Tref(:)
real(rkind), allocatable :: Ea(:)
!--Optics-----------------------
real(rkind) Kw
real(rkind) Kcdom
real(rkind) Kspm
real(rkind) Kchla
!--in module LIGHT_VARS
real(rkind) astar490
real(rkind) aw490
real(rkind) astarOMA
real(rkind) astarOMZ
real(rkind) astarOMR
real(rkind) astarOMBC
real(rkind) PARfac
real(rkind) sinkCDOM
!---Phytoplankton 
real(rkind), allocatable :: ediblevector(:,:)
real(rkind), allocatable :: umax(:)
real(rkind), allocatable :: CChla(:)
real(rkind), allocatable :: alpha(:)
real(rkind), allocatable :: beta(:)
real(rkind), allocatable :: respg(:)
real(rkind), allocatable :: respb(:)
real(rkind), allocatable :: QminN(:)
real(rkind), allocatable :: QminP(:)
real(rkind), allocatable :: QmaxN(:)
real(rkind), allocatable :: QmaxP(:)
real(rkind), allocatable :: Kn(:)
real(rkind), allocatable :: Kp(:)
real(rkind), allocatable :: Ksi(:)
real(rkind), allocatable :: KQn(:)
real(rkind), allocatable :: KQp(:)
real(rkind), allocatable :: nfQs(:)
real(rkind), allocatable :: vmaxN(:)
real(rkind), allocatable :: vmaxP(:)
real(rkind), allocatable :: vmaxSi(:)
real(rkind), allocatable :: aN(:)
real(rkind), allocatable :: volcell(:)
real(rkind), allocatable :: Qc(:)
real(rkind), allocatable :: Athresh(:)
real(rkind), allocatable :: mA(:)
real(rkind), allocatable :: A_wt(:)
!---Zooplankton
real(rkind), allocatable :: Zeffic(:)
real(rkind), allocatable :: Zslop(:)
real(rkind), allocatable :: Zvolcell(:)
real(rkind), allocatable :: ZQc(:)
real(rkind), allocatable :: ZQn(:)
real(rkind), allocatable :: ZQp(:)
real(rkind), allocatable :: optNP(:)
real(rkind), allocatable :: ZKa(:)
real(rkind), allocatable :: Zrespg(:)
real(rkind), allocatable :: Zrespb(:)
real(rkind), allocatable :: Zumax(:)
real(rkind), allocatable :: Zm(:)

!---Organic Matter              
real(rkind) KG1
real(rkind) KG2
real(rkind) KG1R
real(rkind) KG2R
real(rkind) KG1BC
real(rkind) KG2BC
real(rkind) KNH4
real(rkind) nitmax
real(rkind) KO2
real(rkind) KstarO2
real(rkind) KNO3
real(rkind) KGcdom
real(rkind) CF_SPM
!----Other including Boundary Conditions------------
real(rkind) KH_coeff
integer Which_Outer_BC
real(rkind) wt_l, wt_o
real(rkind) wt_pl, wt_po
real(rkind) m_OM_init,m_OM_BC,m_OM_sh
real(rkind) KGbot

!Light curve parameters
real(rkind), allocatable :: alphad(:) ! Initial slope of photosynthesis-irradiance curve / Vmax
real(rkind), allocatable :: betad(:)  ! Photoinhibition constant / Vmax

!Initialize variables in cgem_init from cgem_read
real(rkind), dimension(:), allocatable, private :: sinkA
real(rkind), private :: sinkOM1A,sinkOM2A,sinkOM1Z,sinkOM2Z,sinkOM1R,sinkOM2R,sinkOM1BC,sinkOM2BC

contains

subroutine cgem_dim

  integer           :: istat,iunit
  character(len=1000) :: line
  !http://degenerateconic.com/namelist-error-checking.html
#if defined(DEBUG)
write(6,*) "Begin cgem_dim"
#endif

  namelist /nosp/ nospA,nospZ,skipcgem,checkwindrad,debug,sinkwcgem,cgemcoords,cgemlat,cgemlon

  open(action='read',file='cgem.nml',iostat=istat,newunit=iunit)

  !namelist /switches/
  read(nml=nosp,iostat=istat,unit=iunit)
  if (istat /= 0) then
   backspace(iunit)
   read(iunit,fmt='(A)') line
   write(6,'(A)') &
        'Invalid line in namelist nosp: '//trim(line)
   stop
  endif

  close(iunit)

  if(nospZ.ne.1.and.nospZ.ne.2) then
    write(6,*) "Z's Please, are 1 or 2 for now.  You chose:",nospZ
    write(6,*) "We're going to use nosp2=2 and keep going."
    nospZ=2
  endif

return
end subroutine cgem_dim


subroutine cgem_read

  integer           :: istat,iunit
  character(len=1000) :: line
  !http://degenerateconic.com/namelist-error-checking.html
  namelist /switches/ Which_fluxes,Which_temperature,Which_uptake,Which_quota,&
    Which_photosynthesis,Which_growth,Which_wind,Which_rad,sinkTr
  namelist /optics/ Kw,Kcdom,Kspm,Kchla,astar490,aw490,astarOMA,astarOMZ,astarOMR,astarOMBC,PARfac,sinkCDOM
  namelist /temperature/ Tref,KTg1,KTg2,Ea
  namelist /phytoplankton/ umax,CChla,alpha,beta,respg,respb,QminN,QminP,QmaxN,QmaxP,Kn,Kp,Ksi,KQn,&
     KQp,nfQs,vmaxN,vmaxP,vmaxSi,aN,volcell,Qc,Athresh,sinkA,mA,A_wt
  namelist /zooplankton/ ediblevector,Zeffic,Zslop,Zvolcell,ZQc,ZQn,ZQp,ZKa,Zrespg,Zrespb,Zumax,Zm
  namelist /OM/ KG1,KG2,KG1R,KG2R,KG1BC,KG2BC,KNH4,nitmax,KO2,KstarO2,KNO3,pCO2,&
     sx1R,sy1R,sx2R,sy2R,sx1BC,sy1BC,sx2BC,sy2BC,sinkOM1A,sinkOM2A,sinkOM1Z,sinkOM2Z,sinkOM1R,&
     sinkOM2R,sinkOM1BC,sinkOM2BC,KGcdom,CF_SPM,KGbot

  open(action='read',file='cgem.nml',iostat=istat,newunit=iunit)

  !namelist /switches/
  read(nml=switches,iostat=istat,unit=iunit)
  if (istat /= 0) then
   backspace(iunit)
   read(iunit,fmt='(A)') line
   write(6,'(A)') &
        'Invalid line in namelist: '//trim(line)
   stop
  endif

#if defined(DEBUG)
 write(6,nml=switches)
#endif

  !namelist /optics/
  read(nml=optics,iostat=istat,unit=iunit)
  if (istat /= 0) then
   backspace(iunit)
   read(iunit,fmt='(A)') line
   write(6,'(A)') &
        'Invalid line in namelist optics: '//trim(line)
   stop
  endif

#if defined(DEBUG)
 write(6,nml=optics)
#endif


  !namelist /temperature/
  read(nml=temperature,iostat=istat,unit=iunit)
  if (istat /= 0) then
   backspace(iunit)
   read(iunit,fmt='(A)') line
   write(6,'(A)') &
        'Invalid line in namelist temperature: '//trim(line)
   stop
  endif

#if defined(DEBUG)
 write(6,nml=temperature)
#endif


  !namelist /phytoplankton/
  read(nml=phytoplankton,iostat=istat,unit=iunit)
  if (istat /= 0) then
   backspace(iunit)
   read(iunit,fmt='(A)') line
   write(6,'(A)') &
        'Invalid line in namelist phytoplankton: '//trim(line)
   stop
  endif

#if defined(DEBUG)
 write(6,nml=phytoplankton)
#endif

  !namelist /zooplankton/
  read(nml=zooplankton,iostat=istat,unit=iunit)
  if (istat /= 0) then
   backspace(iunit)
   read(iunit,fmt='(A)') line
   write(6,'(A)') &
        'Invalid line in namelist zooplankton: '//trim(line)
   stop
  endif
  optNP = ZQn/ZQp    ! Optimal nutrient ratio for zooplankton

#if defined(DEBUG)
 write(6,nml=zooplankton)
#endif

  !namelist /OM/
  read(nml=OM,iostat=istat,unit=iunit)
  if (istat /= 0) then
   backspace(iunit)
   read(iunit,fmt='(A)') line
   write(6,'(A)') &
        'Invalid line in namelist OM: '//trim(line)
   stop
  endif

#if defined(DEBUG)
 write(6,nml=OM)
#endif

  close(iunit)

return
end subroutine cgem_read

subroutine cgem_allocate()

integer :: i,ierr
integer :: counter = 0

!---------------------------------------------------------      
!-A; Phytoplankton number density (cells/m3);
!---------------------------------------------------------  
       allocate (iA(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iA"
       do i=1,nospA
          counter = counter+1
          iA(i) = counter
       enddo
!----------------------------------------------------------------------
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!----------------------------------------------------------------------
       allocate (iQn(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iQn"
       do i=1,nospA
          counter = counter+1
          iQn(i) = counter
       enddo
!----------------------------------------------------------------------
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!----------------------------------------------------------------------
      allocate (iQp(nospA),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iQp"
       do i=1,nospA
          counter = counter+1
          iQp(i) = counter
       enddo
!--------------------------------------------------------------------
!-Z: Zooplankton number density (individuals/m3);
!--------------------------------------------------------------------
      allocate (iZ(nospZ),stat=ierr)
       if(ierr.ne.0) write(6,*) "error in allocating:iZ"
       do i=1,nospZ
          counter = counter+1
          iZ(i) = counter
       enddo
!-------------------------------
!-NO3; Nitrate (mmol-N/m3)
!-------------------------------
      iNO3 = counter+1
!--------------------------------      
!-NH4; Ammonium (mmol-N/m3)
!--------------------------------
      iNH4 = counter+2
!-------------------------------------------        
!-PO4: Phosphate (mmol-P/m3)
!--------------------------------------
      iPO4 = counter+3
!---------------------------------------------------------
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) 
!---------------------------------------------------------
      iDIC = counter+4
!-------------------------------------------        
!-O2: Molecular Oxygen (mmol-O2/m3)
!------------------------------
      iO2 = counter+5
!-------------------------------------------------------------
!-OM1_A: (mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton
!-------------------------------------------------------------
      iOM1CA = counter+6
      iOM1NA = counter+7
      iOM1PA = counter+8
!-----------------------------------------------------------------
!-OM2_A: (mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton 
!------------------------------------------------------------------
      iOM2CA = counter+9
      iOM2NA = counter+10
      iOM2PA = counter+11
!-------------------------------------------------------------
!-OM1_Z:(mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets.
!-------------------------------------------------------------
      iOM1CZ = counter+12
      iOM1NZ = counter+13
      iOM1PZ = counter+14
!-------------------------------------------------        
!-OM2_Z:(mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets.
!-----------------------------------------------
      iOM2CZ = counter+15
      iOM2NZ = counter+16
      iOM2PZ = counter+17
!--------------------------------------------------------------------
!-OM1_R: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter arising from river outflow
!--------------------------------------------------------------------
      iOM1R = counter+18
!-------------------------------------------------      
!-OM2_R: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter arising from river outflow
!--------------------------------------------------------------------
      iOM2R = counter+19
!-------------------------------------------
!-CDOM: (ppb) 
!        -- Colored Dissolved Organic Matter
!-------------------------------------------
      iCDOM = counter+20
!---------------------------------------------
!-Silica: (mmol-Si/m3) 
!        -- Silica
!-------------------------------------------
      iSi = counter+21
!--------------------------------------------------------------------
!-OM1_BC: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter in initial and boundary 
!            conditions 
!--------------------------------------------------------------------
      iOM1BC = counter+22
!-------------------------------------------------
!-OM2_BC: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter in initial and boundary
!            conditions
!--------------------------------------------------------------------
      iOM2BC = counter+23
!-------------------------------------------
!-ALK:  (mmol-HCO3/m3)?
!        -- Alkalinity
!-------------------------------------------
      iALK = counter+24
!-------------------------------------------
!-Tr
     iTr = counter + 25 

!How many state variables
      nf = iTr

!----allocate INPUT_VARS_CGEM

!---Phytoplankton 
allocate( ediblevector(nospZ,nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ediblevector"
allocate( umax(nospA),stat=ierr  )
if(ierr.ne.0) write(6,*) "error in allocating:umax"
allocate( CChla(nospA),stat=ierr  )
if(ierr.ne.0) write(6,*) "error in allocating:CChla"
allocate( alpha(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:alpha"
allocate( beta(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:beta"
allocate( respg(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:respg"
allocate( respb(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:respb"
allocate( QminN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QminN"
allocate( QminP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QminP"
allocate( QmaxN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QmaxN"
allocate( QmaxP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:QmaxP"
allocate( Kn(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Kn"
allocate( Kp(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Kp"
allocate( Ksi(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Ksi"
allocate( KQn(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:KQn"
allocate( KQp(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:KQp"
allocate( nfQs(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:nfQs"
allocate( vmaxN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxN"
allocate( vmaxP(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxP"
allocate( vmaxSi(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:vmaxSi"
allocate( aN(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:aN"
allocate( volcell(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:volcell"
allocate( Qc(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Qc"
allocate( Athresh(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Athresh"
allocate( mA(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:mA"
allocate( A_wt(nospA),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:A_wt"

!---Zooplankton
allocate( Zeffic(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zeffic"
allocate( Zslop(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zslop"
allocate( Zvolcell(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zvolcell"
allocate( ZQc(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQc"
allocate( ZQn(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQn"
allocate( ZQp(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZQp"
allocate( optNP(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:optNP"
allocate( ZKa(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:ZKa"
allocate( Zrespg(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zrespg"
allocate( Zrespb(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zrespb"
allocate( Zumax(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zumax"
allocate( Zm(nospZ),stat=ierr )
if(ierr.ne.0) write(6,*) "error in allocating:Zm"

!Light curve parameters
allocate( alphad(nospA),stat=ierr ) ! Initial slope of photosynthesis-irradiance curve / Vmax
if(ierr.ne.0) write(6,*) "error in allocating:alphad"
allocate( betad(nospA),stat=ierr )  ! Photoinhibition constant / Vmax
if(ierr.ne.0) write(6,*) "error in allocating:betad"


!Temperature parameters for growth rates
allocate(Tref(nospA+nospZ),stat=ierr)                   !Tref(nospA+nospZ): Optimum temperature for growth(C)
if(ierr.ne.0) write(6,*) "error in allocating:Tref"
allocate(KTg1(nospA+nospZ),stat=ierr)                   !KTg1(nospA+nospZ): Effect of T below Topt(C^2)
if(ierr.ne.0) write(6,*) "error in allocating:KTg1"
allocate(KTg2(nospA+nospZ),stat=ierr)                   !KTg2(nospA+nospZ): Effect of T above Topt(C^2)
if(ierr.ne.0) write(6,*) "error in allocating:KTg2"
allocate(Ea(nospA+nospZ),stat=ierr)                     !Ea(nospA+nospZ): Slope of Arrhenius plot(eV)
if(ierr.ne.0) write(6,*) "error in allocating:Ea"

!Sinking
allocate(ws(nf),stat=ierr)
if(ierr.ne.0) write(6,*) "error in allocating:ws"
allocate(fmin(nf),stat=ierr)
if(ierr.ne.0) write(6,*) "error in allocating:fmin"

!sinking
allocate(sinkA(nospA),stat=ierr)
if(ierr.ne.0) write(6,*) "error in allocating:sinkA"


return
end subroutine cgem_allocate

subroutine cgem_init

integer :: isp
real(rkind) tot

Athresh = Athresh*volcell   ! Threshold for grazing, um^3/m3
eps=0
do isp=1,nospA
   eps=0
   if(umax(isp).eq.0) eps=1.e-18
   alphad(isp) = alpha(isp)/(umax(isp)+eps) ! Initial slope of photosynthesis-irradiance curve / Vmax
   betad(isp)  = beta(isp)/(umax(isp)+eps)  ! Photoinhibition constant / Vmax
enddo

!Convert relative proportions of phytoplankton to percentage of total chlA
tot = SUM(A_wt)
if(tot.le.0) then
 write(6,*) "Error in A_wt, A_wt.le.0"
 stop
endif

do isp=1,nospA
   A_wt(isp) = A_wt(isp)/tot
enddo

fmin = tiny(x)
fmin(iA(:)) = 1.
fmin(iZ(:)) = 1.
fmin(iQn(:)) = QminN(:)
fmin(iQp(:)) = QminP(:)

ws = 0.
ws(iA(:))=sinkA(:)
ws(iQn(:)) = sinkA(:)
ws(iQp(:)) = sinkA(:)

ws(iCDOM) =  sinkCDOM

ws(iOM1CA) = sinkOM1A
ws(iOM1NA) = sinkOM1A
ws(iOM1PA) = sinkOM1A

ws(iOM2CA) = sinkOM2A
ws(iOM2NA) = sinkOM2A
ws(iOM2PA) = sinkOM2A

ws(iOM1CZ) = sinkOM1Z
ws(iOM1NZ) = sinkOM1Z
ws(iOM1PZ) = sinkOM1Z

ws(iOM2CZ) = sinkOM2Z
ws(iOM2NZ) = sinkOM2Z
ws(iOM2PZ) = sinkOM2Z

ws(iOM1R) = sinkOM1R
ws(iOM2R) = sinkOM2R

ws(iOM1BC) = sinkOM1BC
ws(iOM2BC) = sinkOM2BC

ws(iTr) = sinkTr

!Convert per m/d to m/s 
ws = ws / 86400. 

return
end subroutine cgem_init


end module cgem
