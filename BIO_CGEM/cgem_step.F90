!======================================================================     
  Subroutine cgem_step( ff,ff_new,dT, S, T, PAR, Wind, lat, dz, d_sfc, is_surface, is_bottom, Rad, inea )

!======================================================================
  use cgem
  use calcAgrow
  use cgem_utils
  use date_time
  use mvars
  use schism_glbl, only : rkind

  IMPLICIT NONE

!---------------------------------------------
! Interface variables
!---------------------------------------------------------------------
  logical, intent(in)  :: is_bottom ! Bottom fluxes 
  logical, intent(in)  :: is_surface ! Surface fluxes 
  integer, intent(in)  :: inea !element
  real(rkind), intent(in)     :: lat       ! For mocsy- latitude 
  real(rkind), intent(in)     :: d_sfc     ! For mocsy
  real(rkind), intent(in)     :: dz        ! For fluxes 
  real(rkind), intent(in)     :: PAR       ! PAR at cell center
  real(rkind), intent(in)     :: Wind 
  real(rkind), intent(in)     :: S,T,dT,Rad
  real(rkind),intent(in), dimension(nf) :: ff
  real(rkind),intent(out), dimension(nf) :: ff_new
!---------------------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------
  integer :: isp, isz ! Loop indicies, isp/isz is for phytoplankton/zooplankton species
!------------------------------------ 
! Phytoplankton parameters
! Phytoplankton uptake and growth
  real(rkind), dimension(nospA) :: A     ! Zooplankton number density (indv./m3)
  real(rkind), dimension(nospA) :: Agrow ! Phytoplankton growth (cells/m3/d)
  real(rkind), dimension(nospA) :: Aresp ! Total respiration from a phytoplankton group (cells/m3/d)
  real(rkind), dimension(nospA) :: uA    ! Specific growth rate (1/d)
  real(rkind), dimension(nospA) :: uN    ! Nitrogen Limited growth rate (1/d)
  real(rkind), dimension(nospA) :: uP    ! Phosphorus limited growth rate (1/d)
  real(rkind), dimension(nospA) :: uE    ! Light limited growth rate (1/d)
  real(rkind), dimension(nospA) :: uSi   ! Silica limited growth rate (1/d)
  real(rkind), dimension(nospA) :: Qn    ! Phytoplankton Nitrogen Quota (mmol-N/cell)
  real(rkind), dimension(nospA) :: Qp    ! Phytoplankton Phosphorus Quota (mmol-P/cell)
  real(rkind), dimension(nospA) :: f_Qn  ! Quota model for N
  real(rkind), dimension(nospA) :: f_Qp  ! Quota model for P
  real(rkind), dimension(nospA) :: vN    ! Phytoplankton uptake rate of Nitrogen (mmol-N/cell/d)
  real(rkind), dimension(nospA) :: vP    ! Phytoplankton uptake rate of Phosphorus (mmol-P/cell/d)
  real(rkind), dimension(nospA) :: vSi   ! Phytoplankton uptake rate of Silica (mmol-Si/cell/d)
  real(rkind)                   :: AupN  ! Total Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
  real(rkind)                   :: AupP  ! Total Phytoplankton uptake of Phosphorus (mmol-P/m3/d)
  real(rkind)                   :: AupSi ! Total Phytoplankton uptake of Silica (mmol-Si/m3/d)
! Monod equations for phytoplankton
  real(rkind), dimension(nospA) :: monodN  !Monod term in nitrogen uptake
  real(rkind), dimension(nospA) :: monodP  !Monod term in phosphorus uptake
  real(rkind), dimension(nospA) :: monodSi !Monod term in Si uptake
  real(rkind)                   :: Ntotal   ! Total N (mmol/m3)
! Phytoplankton nutrient loss
  real(rkind), dimension(nospA) :: Amort   ! Dead phytoplankton (cells/m3/day)
  real(rkind)  :: AexudN  ! Sum of Exudation of N from all phytoplankton groups (mmol-N/m3/d)
  real(rkind)  :: AexudP  ! Sum of Exudation of P from all phytoplankton groups (mmol-P/m3/d)
  real(rkind)  :: ArespC  ! Phytoplankton equivalent carbon loss from respiration (mmol-C/m3/d)
!------------------------------------------------------------------
! Zooplankton parameters
! Zooplankton uptake and growth
  real(rkind), dimension(nospZ)       :: Z          ! Zooplankton number density (indv./m3)
  real(rkind), dimension(nospZ)       :: Zgrow      ! Zooplankton growth (indv./m3/d)
  real(rkind), dimension(nospA,nospZ) :: Zgrazvol   ! Grazing rate in units of biovolume (um3/m3/d)
  real(rkind), dimension(nospA,nospZ) :: ZgrazA     ! Zooplankton grazing of phytoplankton (cells/m3/d)
  real(rkind), dimension(nospA)       :: ZgrazA_tot ! Total zooplankton grazing of phytoplankton (cells/m3/d)
  real(rkind), dimension(nospZ)       :: ZgrazN     ! Zooplankton grazing uptake of Nitrogen (mmol-N/m3/d)
  real(rkind), dimension(nospZ)       :: ZgrazP     ! Zooplankton grazing uptake of Phosphorus (mmol-P/m3/d)
  real(rkind), dimension(nospZ)       :: ZgrazC     ! Zooplankton grazing uptake of Carbon (mmol-C/m3/d)
  real(rkind), dimension(nospZ)       :: ZinN       ! Zooplankton ingestion of Nitrogen (mmol-N/m3/d)
  real(rkind), dimension(nospZ)       :: ZinP       ! Zooplankton ingestion of Phosphorus (mmol-P/m3/d)
  real(rkind), dimension(nospZ)       :: ZinC       ! Zooplankton ingestion of Carbon (mmol-C/m3/d)
!Monod equations for zooplankton ingestion of phytoplankton
  real(rkind), dimension(nospA,nospZ) :: monodZ     ! Monod term for zooplankton grazing
  real(rkind)                         :: Abiovol    ! Algae biovolume vector (um3/m3)
  real(rkind), dimension(nospA,nospZ) :: top_A      ! Monod numerator value for phytoplankton group
  real(rkind), dimension(nospA,nospZ) :: bottom_A   ! Monod Denominator value for phytoplankton group
  real(rkind), dimension(nospZ)       :: bottom     ! Sum of Monod Denominator value for all phytoplankton groups
!Zooplankton nutrient loss
  real(rkind), dimension(nospZ)       :: Zresp      ! Zooplankton respiration (individuals/m3/d)
  real(rkind)                         :: ZrespC     ! Carbon loss from zooplankton respiration (mmol-C/m3/day)
  real(rkind), dimension(nospZ)       :: ZunC       ! Unassimilated ingested Carbon (mmol-C/m3/d)
  real(rkind), dimension(nospZ)       :: ZunN       ! Unassimilated ingested Nitrogen (mmol-N/m3/d)
  real(rkind), dimension(nospZ)       :: ZunP       ! Unassimilated ingested Phosphorus (mmol-P/m3/d)
  real(rkind), dimension(nospZ)       :: ZunSi      ! Unassimilated ingested Silica (mmol-Si/m3/d)
  real(rkind), dimension(nospZ)       :: Zmort      ! Dead zooplankton (individuals/m3/d)
  real(rkind), dimension(nospZ)       :: ZmortC     ! Carbon released from dead zooplankton (mmol-C/m3/d)
  real(rkind), dimension(nospZ)       :: ZmortN     ! Nitrogen released from dead zooplankton (mmol-N/m3/d)
  real(rkind), dimension(nospZ)       :: ZmortP     ! Phosphorus released from dead zooplankton (mmol-P/m3/d)
  real(rkind), dimension(nospZ)       :: ZslopC     ! Carbon lost to sloppy feeding (mmol-C/m3/d)
  real(rkind), dimension(nospZ)       :: ZslopN     ! Nitrogen lost to sloppy feeding (mmol-N/m3/d)
  real(rkind), dimension(nospZ)       :: ZslopP     ! Phosphorus lost to sloppy feeding (mmol-P/m3/d)
  real(rkind), dimension(nospZ)       :: ZexN       ! Excretion from zooplankton (mmol-N/m3/d)
  real(rkind), dimension(nospZ)       :: ZexP       ! Excretion from zooplankton (mmol-P/m3/d)
  real(rkind), dimension(nospZ)       :: ZegC       ! Egestion from zooplankton (mmol-C/m3/d)
  real(rkind), dimension(nospZ)       :: ZegN       ! Egestion from zooplankton (mmol-N/m3/d)
  real(rkind), dimension(nospZ)       :: ZegP       ! Egestion from zooplankton (mmol-P/m3/d)
  real(rkind), dimension(nospZ)       :: ZegSi      ! Egestion from zooplankton (mmol-Si/m3/d)
  real(rkind)          :: OM1_Ratio, OM2_Ratio      ! Separates sloppy feeding into OM1 and OM2
!---------------------------------------------------------------------- 
! Time variables  
  real(rkind), parameter :: one_d_365  = 1.0/365.0  ! Convert 1/yr to 1/day
  real(rkind), parameter :: SDay       = 86400.0     ! Seconds in 1 day
!-----------------------------------------------------------------------
! Organic Matter Calculations
! Variables to calculate stoichiometry C:N:P ratios
  real(rkind)             :: OM1_CA, OM1_NA, OM1_PA ! OM from dead phytoplankton
  real(rkind)             :: OM2_CA, OM2_NA, OM2_PA   
  real(rkind)             :: OM1_CZ, OM1_NZ, OM1_PZ ! OM from zooplankton
  real(rkind)             :: OM2_CZ, OM2_NZ, OM2_PZ
  real(rkind)             :: sx1,sy1,sx2,sy2        ! x=C/P, y=N/P 
  real(rkind), parameter  :: sz = 1.0                ! z=P/P = 1
!---------------------------------------------------------------------------
! reaction and Nitrification subroutine variables
  real(rkind)    :: R_11                            ! Nitrification term
  real(rkind)    :: RNO3_A, RNO3_Z, RNO3_R, RNO3_BC ! Remineralization terms for NO3
  real(rkind)    :: RNH4_A, RNH4_Z, RNH4_R, RNH4_BC ! Remineralization terms for NH4
  real(rkind)    :: ROM1CA, ROM1CZ, ROM1_R, ROM1_BC ! Remineralization terms for POC
  real(rkind)    :: ROM2CA, ROM2CZ, ROM2_R, ROM2_BC ! Remineralization terms for DOC
  real(rkind)    :: ROM1NA, ROM1PA, ROM1NZ, ROM1PZ  !
  real(rkind)    :: ROM2NA, ROM2PA, ROM2NZ, ROM2PZ  !
  real(rkind)    :: RO2_A, RO2_Z, RO2_R, RO2_BC     ! Remineralization terms for O2
  real(rkind)    :: RPO4_A, RPO4_Z, RPO4_R, RPO4_BC ! Remineralization terms for PO4
  real(rkind)    :: RDIC_A, RDIC_Z, RDIC_R, RDIC_BC ! Remineralization terms for DIC
  real(rkind)    :: RSi_A, RSi_Z, RSi_R, RSi_BC     ! Remineralization terms for Si
  real(rkind)    :: RALK_A, RALK_Z, RALK_R, RALK_BC ! Remineralization terms for ALK
  real(rkind)    :: RN2_A, RN2_Z, RN2_R, RN2_BC     ! Remineralization terms for N2 
  real(rkind), dimension(10) :: RC   ! Array that returns remineralization terms for OM
  real(rkind)    :: RNO3,RNH4,RO2,RPO4,RDIC,RSi,RALK,RN2 ! Totals
!---------------------------------------------------------
  real(rkind)    :: OM1R,OM2R,OM1BC,OM2BC
  real(rkind)    :: CDOM,NO3,NH4,DIC,O2,PO4,Si,ALK
  real(rkind)    :: OM1CA,OM1NA,OM1PA,OM1CZ,OM1NZ,OM1PZ
  real(rkind)    :: OM2CA,OM2NA,OM2PA,OM2CZ,OM2NZ,OM2PZ
!-----------------------------------------------------------------------
! Other variables 
  real(rkind)            :: PrimProd                ! Primary production (photosynthesis)
  real(rkind), dimension(nospA+nospZ) :: Tadj       ! Temperature adjustment factor
  real(rkind) :: rKG1,rKG2
!------------------------------------------------------------------    
!timestep in days
  real(rkind) :: dTd
!------------------------------------------------------------------
!Output vars for alkalinity subroutine:
  real :: ph_calc(1), pco2_calc(1), fco2(1), co2(1), hco3(1), co3(1), omegaa(1), omegac(1), betad_calc(1) 
  real :: rhosw(1), p(1), tempis(1)
  real :: patm(1) = 1.
  real :: m_alk(1), m_dic(1), m_si(1), m_po4(1)
  real :: m_lat(1)
  real :: m_d_sfc(1),m_T(1),m_S(1)
  real :: pH

  !convert to timestep in days
  dTd = dt/SDay

  ! Renaming is for readability...
  A(:)  = ff(iA(:))
  ! After Advection and VMixing, return to Q's
  Qn(:) = ff(iQn(:)) / A(:) 
  Qp(:) = ff(iQp(:)) / A(:) 
  Z(:)  = ff(iZ(:))
  NO3   = ff(iNO3)
  NH4   = ff(iNH4)
  PO4   = ff(iPO4)
  DIC   = ff(iDIC)
  O2    = ff(iO2)
  OM1CA = ff(iOM1CA)
  OM1NA = ff(iOM1NA)
  OM1PA = ff(iOM1PA)
  OM2CA = ff(iOM2CA)
  OM2NA = ff(iOM2NA)
  OM2PA = ff(iOM2PA)
  OM1CZ = ff(iOM1CZ)
  OM1NZ = ff(iOM1NZ)
  OM1PZ = ff(iOM1PZ)
  OM2CZ = ff(iOM2CZ)
  OM2NZ = ff(iOM2NZ)
  OM2PZ = ff(iOM2PZ)
  OM1R  = ff(iOM1R)
  OM2R  = ff(iOM2R)
  CDOM  = ff(iCDOM)
  Si    = ff(iSi)
  OM1BC = ff(iOM1BC)
  OM2BC = ff(iOM2BC)
  ALK   = ff(iALK)
  Ntotal= NO3 + NH4

  !- calc_Agrow calculates: 
  !    Agrow: production (cells/m3/s)
  !    Aresp: sum of somatic and basal respirtion (cells/m3/s)
  call calc_Agrow(PAR,T,Qn,Qp,Si,A,Agrow,uA,Aresp,uN,uP,uE,uSi)
 
  !- ZgrazA_tot: total zooplankton grazing on Ai by all zooplankton groups (cells/m3/d)
  do isp = 1, nospA
     Abiovol         = A(isp)*volcell(isp) 
     top_A(isp,:)    = DMAX1((Abiovol-Athresh(isp))*ediblevector(:,isp),0.0)
     bottom_A(isp,:) = Abiovol * ediblevector(:,isp)
  enddo
 
  do isz = 1, nospZ
     bottom(isz) = SUM(bottom_A(:,isz))   ! sum over isp
  enddo
 
  do isp = 1, nospA
     monodZ(isp,:)  = top_A(isp,:)/(ZKa(:) + bottom(:))
  enddo 

  !--------------------------------------
  !-- Temperature adjustment factor 
  call func_T( T,Tadj )
  !-- Nutrient dependent growth function
  call func_Qs( Qn, Qp, f_Qn, f_Qp)

  !! Sum over phytoplankton groups for PrimProd, ArespC, AexudN, AexudP
  PrimProd = SUM(Agrow(:)*Qc(:)) ! Phytoplankton primary production (mmol-C/m3/d)
  ArespC   = SUM(Aresp(:)*Qc(:)) ! Phytoplankton respiration     (mmol-C/m3/d)
  AexudN   = SUM(Aresp(:)*Qn(:)) ! Total Phytoplankton exudation (mmol-N/m3/d)
  AexudP   = SUM(Aresp(:)*Qp(:)) ! Total Phytoplankton exudation (mmol-P/m3/d)
  Amort(:) = A(:) * mA(:)        ! Dead phytoplankton (cells/m3/d)

  !------------------------------------------------------------------------     
  ! Nutrient limited uptake:
  ! Find Rate Limiting Nutrient RLN for N, P, and Si:
  if(Rad.le.tiny(x)) then  !Nutrient uptake only takes place during the day
     vN = 0.
     vP = 0.
     vSi = 0.
  else!Day
      !Rate limiting nutrient is N
      !    Monod Equations
      ! Kx is half saturation coefficient for x
      do isp = 1, nospA
         monodN(isp)  = Ntotal/(Ntotal+Kn(isp))
         monodP(isp)  = PO4/(PO4+Kp(isp))
         monodSi(isp) = Si/(Si+Ksi(isp))
      enddo
     !--Rate limiting nutrient is N
      if(Ntotal.le.PO4.and.Ntotal.le.Si) then
         do isp = 1, nospA
            vN(isp) = Q10_T(T,vmaxN(isp))*monodN(isp)*f_Qn(isp)

            vP(isp) = Q10_T(T,vmaxP(isp))*monodP(isp)*f_Qp(isp) &
     &      *( Ntotal/(Ntotal+aN(isp)*Kn(isp)) )

            vSi(isp) = Q10_T(T,vmaxSi(isp))*monodSi(isp)        &
     &      *( Ntotal/(Ntotal+aN(isp)*Kn(isp)) )
         enddo

     !--Rate limiting nutrient is P
      elseif(PO4.le.Ntotal.and.PO4.le.Si) then
         do isp = 1, nospA
            vN(isp) = Q10_T(T,vmaxN(isp))*monodN(isp)*f_Qn(isp)&
     &      *( PO4/(PO4+aN(isp)*Kp(isp)) )

            vP(isp) = Q10_T(T,vmaxP(isp))*monodP(isp)*f_Qp(isp)

            vSi(isp) = Q10_T(T,vmaxSi(isp))*monodSi(isp)       &
     &      *( PO4/(PO4+aN(isp)*Kp(isp)) )
         enddo

     !--Rate limiting nutrient is Si 
      else
         do isp = 1, nospA
            vN(isp) = Q10_T(T,vmaxN(isp))*monodN(isp)*f_Qn(isp)&
     &      *( Si/(Si+aN(isp)*Ksi(isp)) )

            vP(isp) = Q10_T(T,vmaxP(isp))*monodP(isp)*f_Qp(isp)&
     &      *( Si/(Si+aN(isp)*Ksi(isp)) )

            vSi(isp) = Q10_T(T,vmaxSi(isp))*monodSi(isp)
         enddo
      endif !-End if rate limiting nutrient

  endif !End if day or night

  !! ----------------------------------------------------------------------
  !! Sum over phytoplankton groups for PrimProd, ArespC, AexudN, AexudP
  AupN  = SUM(A(:)*vN(:))     ! Phytoplankton uptake of Nitrogen (mmol-N/m3/d)
  AupP  = SUM(A(:)*vP(:))     ! Phytoplankton uptake of Phosphorus (mmol-P/m3/d)
  AupSi = SUM(A(:)*vSi(:))  ! Phytoplankton uptake of Silica (mmol-Si/m3/d)

  do isp=1,nospA
      ZgrazA(isp,:)   = (Z(:)*Zumax(:)*monodZ(isp,:))/volcell(isp)   ! Grazing of phytoplankton by plankton (cells/m3/d)
      ZgrazA_tot(isp) = SUM( ZgrazA(isp,:) )
  enddo

  do isz = 1,nospZ
      ZgrazC(isz) = SUM(ZgrazA(:,isz) * Qc(:)) ! Carbon uptake from grazing (mmol-C/m3/day)
      ZgrazN(isz) = SUM(ZgrazA(:,isz) * Qn(:)) ! Nitrogen uptake from grazing( mmol-N/m3/day)
      ZgrazP(isz) = SUM(ZgrazA(:,isz) * Qp(:)) ! Phosphorus uptake from grazing (mmol-P/m3/day)
  enddo

  !---------------------------------------------------------
  !-A; Phytoplankton number density (cells/m3);
  !---------------------------------------------------------
  ff_new(iA(:)) = DMAX1(A(:)        &
  & + ( Agrow(:) - Aresp(:) - ZgrazA_tot(:) - Amort(:) )*dTd,fmin(iA(:)))
  if(debug.eq.2.and.inea.eq.10) write(6,*) "PAR",PAR
  if(debug.eq.2.and.inea.eq.10) write(6,*) "A,Agrow,Aresp,ZgrazA_tot,Amort,dTd",A,Agrow,Aresp,ZgrazA_tot,Amort,dTd
  !----------------------------------------------------------------------
  !-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
  !----------------------------------------------------------------------
  !! Enforce minima for Droop, also enforce maxima if not equal Droop (which_quota=1)
  if(which_quota.eq.1) then
    do isp=1,nospA
       ff_new(iQn(isp)) = DMAX1(Qn(isp) + (vN(isp) - Qn(isp)*uA(isp))*dTd,QminN(isp))
    enddo
  !! , also enforce maxima if not equal Droop (which_quota=1)
  else
    do isp=1,nospA
        ff_new(iQn(isp)) = DMIN1(DMAX1(Qn(isp) + (vN(isp) - Qn(isp)*uA(isp))*dTd,QminN(isp)),QmaxN(isp))
    enddo
  endif
  if(debug.eq.2.and.inea.eq.10) write(6,*) "Qn,vN,Qp,uA",Qn,vN,Qp,vP,uA

  !----------------------------------------------------------------------
  !-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
  !----------------------------------------------------------------------
  !! Enforce minima for Droop, also enforce maxima if not equal Droop (which_quota=1)
  if(which_quota.eq.1) then
    do isp=1,nospA
        ff_new(iQp(isp)) = DMAX1(Qp(isp) + (vP(isp) - Qp(isp)*uA(isp))*dTd,QminP(isp))
    enddo
  !! , also enforce maxima if not equal Droop (which_quota=1)
  else
    do isp=1,nospA
        ff_new(iQp(isp)) = DMIN1(DMAX1(Qp(isp) + (vP(isp) - Qp(isp)*uA(isp))*dTd,QminP(isp)),QmaxP(isp))
    enddo
  endif
  !----------------------------------------------------------------------- 

  !----------------------------------------------------------------------
  ! ZgrazC, ZgrazN, and ZgrazP are total carbon, nitrogen, and phosphorus
  !   uptake of zooplankton from grazing all phytoplankton groups
  !---------------------------------------------------------------------
  !-------------------------------------------------------------------
  ! Now calculate the total ingested ZinC, ZinN, and ZinP of C, N, and P
  !-------------------------------------------------------------------
  ZslopC(:)  = Zslop(:)*ZgrazC(:)                      ! Sloppy feeding (mmol-C/m3/d)
  ZunC(:)    = (1.-Zeffic(:))*(ZgrazC(:)-ZslopC(:))    ! Unassimilated (mmol-C/m3/d)
  ZinC(:)    = ZgrazC(:) - ZslopC(:) - ZunC(:)         ! Ingested (mmol-C/m3/d)

  ZslopN(:)  = Zslop(:)*ZgrazN(:)                      ! Sloppy feeding (mmol-N/m3/d)
  ZunN(:)    = (1.-Zeffic(:))*(ZgrazN(:)-ZslopN(:))    ! Unassimilated (mmol-N/m3/d)
  ZinN(:)    = ZgrazN(:) - ZslopN(:) - ZunN(:)         ! Ingested (mmol-N/m3/d)

  ZslopP(:)  = Zslop(:)*ZgrazP(:)                      ! Sloppy feeding (mmol-P/m3/d)
  ZunP(:)    = (1.-Zeffic(:))*(ZgrazP(:)-ZslopP(:))    ! Unassimilated (mmol-P/m3/d)
  ZinP(:)    = ZgrazP(:) - ZslopP(:) - ZunP(:)         ! Ingested (mmol-P/m3/d)
  !-------------------------------------------------
  if(debug.eq.2.and.inea.eq.10) write(6,*) "ZslopC,ZunC,ZinC",ZslopC,ZunC,ZinC

  !------------------------------------         
  ! Liebigs Law for zooplankton group isz 
  !------------------------------------
  do isz=1,nospZ
  if (ZinN(isz) .gt. optNP(isz)*ZinP(isz)) then  
      Zgrow(isz)= ZinP(isz)/ZQp(isz)                   ! P-limited growth (indv./m3/d) 
      ZegN(isz) = ZinN(isz) - ZinP(isz)*optNP(isz)       ! P-limited N excretion (mmol-N/m3/d) 
                                                 ! determined by subtracting N-equivalent of ZinP
      ZegC(isz) = ZinC(isz) - ZinP(isz)/ZQp(isz)*ZQc(isz)  ! P-limited C excretion (mmol-C/m3/d)
      ZegP(isz) = 0.                        
  else
      Zgrow(isz)= ZinN(isz)/ZQn(isz)                   ! N-limited growth (indv./m3/d)
      ZegP(isz) = ZinP(isz) - ZinN(isz)/optNP(isz)       ! N-limited P excretion (mmol-P/m3/d)    
                                                 ! determined by subtracting P-equivalent of ZinN
      ZegC(isz) = ZinC(isz) - ZinN(isz)/ZQn(isz)*ZQc(isz)  ! N-limited C excretion (mmol-C/m3/d)
      ZegN(isz) = 0.
  endif
  enddo

  !------------------------------------------------

  !-----------------------------------------------------
  ! ZegC should not be negative 
  do isz=1,nospZ
      if(ZegC(isz).lt.0.0) then
          !ZegC(isz) = 0.0
          !ZegN(isz) = 0.
          !ZegP(isz) = 0.
          write(6,*) "ZegC.lt.0: ZegC,ZegN,ZegP",ZegC(isz),ZegN(isz),ZegP(isz)
          stop
      endif
  enddo

  ! Egestion and unassimilated for Si set equivalent to that of N
  ZegSi(:)  = ZegN(:)
  ZunSi(:)  = ZunN(:)
  Zresp(:)  = (Zgrow(:)*Zrespg(:) + Z(:)*Zrespb(:)) !Zooplankton respiration (indv./m3/d)
  ZexN(:)   = Zresp(:)*ZQn(:)               ! (mmol-N/m3/d)
  ZexP(:)   = Zresp(:)*ZQp(:)               ! (mmol-P/m3/d)
  Zmort(:)  = Zm(:) * Z(:) * Z(:)     ! (indv./m3/d)
  ZmortC(:) = Zmort(:)*ZQc(:)           ! (mmol-C/m3/d)
  ZmortN(:) = Zmort(:)*ZQn(:)           ! (mmol-N/m3/d)
  ZmortP(:) = Zmort(:)*ZQp(:)           ! (mmol-P/m3/d)

  ZrespC     = SUM(Zresp(:)*ZQc(:))

  !-------------------------------------------------------------------------

  !---------------------------------------------------------
  !-Z; Zooplankton number density (individuals/m3);
  !---------------------------------------------------------
  ff_new(iZ(:))  = DMAX1( Z(:)                         &
  &      + (Zgrow(:) - Zresp(:) - Zmort(:))*dTd, fmin(iZ(:)))
  if(debug.eq.2.and.inea.eq.10) write(6,*) "Z,Zgrow,Zresp,Zmort",Z,Zgrow,Zresp,Zmort

  !-----------------------------------------------------------
  ! Remineralization - reactions
  !---------------------------------------------------------------
  ! Instant Remineralization, if on bottom of shelf, redefine KG's
  rKG1 = KG1
  rKG2 = KG2
  if(is_bottom.and.which_fluxes(iInRemin).eq.1) then
           rKG1 = KGbot
           rKG2 = KGbot
  endif
  !------------------------------------------------------------
  ! Nitrification
  !--------------------------------------------------------------
  call Nitrification( O2, NH4, KO2, KNH4, nitmax, T, R_11 )

  !------------------------------------------------------------
  ! Carbon Chemistry
  !--------------------------------------------------------------
  !!! MOCSY alkalinity expressions:
  m_alk = ALK/1000.d0
  m_dic = DIC/1000.d0
  m_si  = Si/1000.d0
  m_po4 = PO4/1000.d0
  m_S = S
  m_T = T
  m_d_sfc = d_sfc
  m_lat = lat

  call vars(ph_calc, pco2_calc, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD_calc, rhoSW, p, tempis,&
  &    m_T, m_S, m_alk, m_dic, m_si, m_po4, patm, m_d_sfc, m_lat, 1, &
  &    'mol/m3', 'Tinsitu', 'm ', 'u74', 'l  ', 'pf ', 'Pzero  ')
  if(debug.eq.3) write(6,*) "ph, mocsy",ph_calc, pco2_calc, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD_calc, rhoSW, p, tempis,&
  &    m_T, m_S, m_alk, m_dic, m_si, m_po4, patm, m_d_sfc, m_lat
  pH = ph_calc(1)

  !------------------------------------------------------------
  ! Particulate and Dissolved dead phytoplankton, rate of remineralization
  !--------------------------------------------------------------
  sx1 = OM1CA/OM1PA
  sy1 = OM1NA/OM1PA
  sx2 = OM2CA/OM2PA
  sy2 = OM2NA/OM2PA
  call reaction( OM1CA, OM2CA, O2, NO3, rKG1, rKG2, KO2, KstarO2, KNO3,     &
  &  sx1, sy1, sz, sx2, sy2, sz, T, RC )
  RC     = one_d_365 * RC  !Change units from /year to /day
  ROM1CA = RC(1)           ! units are /m3/day
  ROM2CA = RC(2)
  ROM1NA = RC(1)*sy1/sx1
  ROM2NA = RC(2)*sy2/sx2
  ROM1PA = RC(1)/sx1
  ROM2PA = RC(2)/sx2
  RO2_A  = RC(3)
  RNO3_A = RC(4)
  RPO4_A = RC(5)
  RDIC_A = RC(6)
  RNH4_A = RC(7)
  RSi_A  = RC(8)
  RALK_A = RC(9)
  RN2_A  = RC(10)

  !------------------------------------------------------------
  ! Particulate and Dissolved fecal pellets, rate of remineralization
  !--------------------------------------------------------------
  sx1 = OM1CZ/OM1PZ
  sy1 = OM1NZ/OM1PZ
  sx2 = OM2CZ/OM2PZ
  sy2 = OM2NZ/OM2PZ
  call reaction( OM1CZ, OM2CZ, O2, NO3, rKG1, rKG2, KO2, KstarO2, KNO3,       &
  &  sx1, sy1, sz, sx2, sy2, sz, T, RC )
  RC     = one_d_365 * RC  !Change units from /year to /day
  ROM1CZ = RC(1)           ! units are /m3/day
  ROM2CZ = RC(2)
  ROM1NZ = RC(1)*sy1/sx1
  ROM2NZ = RC(2)*sy2/sx2
  ROM1PZ = RC(1)/sx1
  ROM2PZ = RC(2)/sx2
  RO2_Z  = RC(3)
  RNO3_Z = RC(4)
  RPO4_Z = RC(5)
  RDIC_Z = RC(6)
  RNH4_Z = RC(7)
  RSi_Z  = RC(8)
  RALK_Z = RC(9)
  RN2_Z  = RC(10)

  !------------------------------------------------------------
  ! Particulate and Dissolved riverine OM, rate of remineralization 
  !------------------------------------------------------------
  call reaction( OM1R, OM2R, O2, NO3, KG1R, KG2R, KO2, KstarO2, KNO3,               &
  &  sx1R, sy1R, sz, sx2R, sy2R, sz, T, RC )
  RC     = one_d_365 * RC  !Change units from /year to /day
  ROM1_R = RC(1)           ! units are /m3/day
  ROM2_R = RC(2)
  RO2_R  = RC(3)
  RNO3_R = RC(4)
  RPO4_R = RC(5)
  RDIC_R = RC(6)
  RNH4_R = RC(7)
  RSi_R  = RC(8)
  RALK_R = RC(9)
  RN2_R  = RC(10)

  !------------------------------------------------------------
  ! Particulate and Dissolved initial and boundary OM, rate of remineralization
  !------------------------------------------------------------
  call reaction( OM1BC, OM2BC, O2, NO3, KG1BC, KG2BC, KO2, KstarO2, KNO3,               &
  &  sx1BC, sy1BC, sz, sx2BC, sy2BC, sz, T, RC )
  RC      = one_d_365 * RC  !Change units from /year to /day
  ROM1_BC = RC(1)           ! units are /m3/day
  ROM2_BC = RC(2)
  RO2_BC  = RC(3)
  RNO3_BC = RC(4)
  RPO4_BC = RC(5)
  RDIC_BC = RC(6)
  RNH4_BC = RC(7)
  RSi_BC  = RC(8)
  RALK_BC = RC(9)
  RN2_BC  = RC(10)

  !--------------------------------------------------------------------
  ! Sum remineralization terms from dead phytoplankton, fecal pellets, and riverine particulate
  RO2   = RO2_A  + RO2_Z  + RO2_R  + RO2_BC  - 2.*R_11 ! (mmol-O2/m3/d)
  RNO3  = RNO3_A + RNO3_Z + RNO3_R + RNO3_BC + R_11    ! (mmol-NO3/m3/d)
  RNH4  = RNH4_A + RNH4_Z + RNH4_R + RNH4_BC - R_11    ! (mmol-NH4/m3/d)
  RPO4  = RPO4_A + RPO4_Z + RPO4_R + RPO4_BC           ! (mmol-PO4/m3/d)
  RDIC  = RDIC_A + RDIC_Z + RDIC_R + RDIC_BC           ! (mmol-DIC/m3/d)
  RSi   = RSi_A  + RSi_Z  + RSi_R  + RSi_BC            ! (mmol-Si/m3/d)
  RALK  = RALK_A + RALK_Z + RALK_R + RALK_BC - 2.*R_11 ! (mmol-HCO3/m3/d)
  RN2   = RN2_A  + RN2_Z  + RN2_R  + RN2_BC            ! (mmol-N2/m3/d)
  !--------------------------------------------------------------------
  if(debug.eq.2.and.inea.eq.10) write(6,*) "O,N,NH,P,D,S,N2,11",RO2,RNO3,RNH4,RPO4,RDIC,RSi,RALK,RN2,R_11
  !---------------------------------------------------------------------
  ! Stoichiometry - calculate C:N:P ratios for Remineralization equations
  !---------------------------------------------------------------------
  !-- Organic Matter from dead phytoplankton --------------------------
  OM1_CA = 0.0
  OM2_CA = 0.0
  OM1_NA = 0.0
  OM2_NA = 0.0
  OM1_PA = 0.0
  OM2_PA = 0.0

  do isp=1,nospA
   !If Nitrogen limited
   if ( uN(isp) .lt. uP(isp) ) then
    !Particulate
     OM1_CA = OM1_CA + Amort(isp)*(Qn(isp)-QminN(isp))/Qn(isp)*Qc(isp)
     OM1_NA = OM1_NA + Amort(isp)*(Qn(isp)-QminN(isp))
     OM1_PA = OM1_PA + Amort(isp)*(Qn(isp)-QminN(isp))/Qn(isp)*Qp(isp)
    !Dissolved
     OM2_CA = OM2_CA + Amort(isp)*QminN(isp)/Qn(isp)*Qc(isp)
     OM2_NA = OM2_NA + Amort(isp)*QminN(isp)
     OM2_PA = OM2_PA + Amort(isp)*QminN(isp)/Qn(isp)*Qp(isp)
   else !Phosphorus limited
    !Particulate
     OM1_CA = OM1_CA + Amort(isp)*(Qp(isp)-QminP(isp))/Qp(isp)*Qc(isp)
     OM1_NA = OM1_NA + Amort(isp)*(Qp(isp)-QminP(isp))/Qp(isp)*Qn(isp)
     OM1_PA = OM1_PA + Amort(isp)*(Qp(isp)-QminP(isp))
    !Dissolved
     OM2_CA = OM2_CA + Amort(isp)*QminP(isp)/Qp(isp)*Qc(isp)
     OM2_NA = OM2_NA + Amort(isp)*QminP(isp)/Qp(isp)*Qn(isp)
     OM2_PA = OM2_PA + Amort(isp)*QminP(isp)
   endif
  enddo

!!-- Organic Matter from fecal pellets ---------------------------------
  OM1_Ratio = SUM( (Qn(:)-QminN(:))/Qn(:)*A(:)) / SUM(A(:))
  OM2_Ratio = SUM(       (QminN(:)/Qn(:))*A(:)) / SUM(A(:))

  if(nospZ.eq.1) then 
    ! Particulate
     OM1_CZ  = .5*(ZegC(1) + ZunC(1) + SUM(ZmortC(:))) + OM1_Ratio*SUM(ZslopC(:)) !(mmol-C/m3/d)
     OM1_NZ  = .5*(ZegN(1) + ZunN(1) + SUM(ZmortN(:))) + OM1_Ratio*SUM(ZslopN(:)) !(mmol-N/m3/d)
     OM1_PZ  = .5*(ZegP(1) + ZunP(1) + SUM(ZmortP(:))) + OM1_Ratio*SUM(ZslopP(:)) !(mmol-P/m3/d)
    ! Dissolved
     OM2_CZ  = .5*(ZegC(1) + ZunC(1) + SUM(ZmortC(:))) + OM2_Ratio*SUM(ZslopC(:)) !(mmol-C/m3/d)
     OM2_NZ  = .5*(ZegN(1) + ZunN(1) + SUM(ZmortN(:))) + OM2_Ratio*SUM(ZslopN(:)) !(mmol-N/m3/d)
     OM2_PZ  = .5*(ZegP(1) + ZunP(1) + SUM(ZmortP(:))) + OM2_Ratio*SUM(ZslopP(:)) !(mmol-P/m3/d)
  else if(nospZ.eq.2) then
    ! Particulate
     OM1_CZ  = ZegC(1) + ZunC(1) + SUM(ZmortC(:)) + OM1_Ratio*SUM(ZslopC(:)) !(mmol-C/m3/d)
     OM1_NZ  = ZegN(1) + ZunN(1) + SUM(ZmortN(:)) + OM1_Ratio*SUM(ZslopN(:)) !(mmol-N/m3/d)
     OM1_PZ  = ZegP(1) + ZunP(1) + SUM(ZmortP(:)) + OM1_Ratio*SUM(ZslopP(:)) !(mmol-P/m3/d)
    ! Dissolved
     OM2_CZ  = ZegC(2) + ZunC(2) + OM2_Ratio*SUM(ZslopC(:))  !(mmol-C/m3/d)
     OM2_NZ  = ZegN(2) + ZunN(2) + OM2_Ratio*SUM(ZslopN(:))  !(mmol-N/m3/d)
     OM2_PZ  = ZegP(2) + ZunP(2) + OM2_Ratio*SUM(ZslopP(:))  !(mmol-P/m3/d)
  else 
    ! Particulate
     OM1_CZ  = ZegC(1) + ZunC(1) + SUM(ZmortC(:)) + OM1_Ratio*SUM(ZslopC(:))      !(mmol-C/m3/d)
     OM1_NZ  = ZegN(1) + ZunN(1) + SUM(ZmortN(:)) + OM1_Ratio*SUM(ZslopN(:))      !(mmol-N/m3/d)
     OM1_PZ  = ZegP(1) + ZunP(1) + SUM(ZmortP(:)) + OM1_Ratio*SUM(ZslopP(:))      !(mmol-P/m3/d)
    ! Dissolved
     OM2_CZ  = SUM(ZegC(2:nospZ)) + SUM(ZunC(2:nospZ)) + OM2_Ratio*SUM(ZslopC(:)) !(mmol-C/m3/d)
     OM2_NZ  = SUM(ZegN(2:nospZ)) + SUM(ZunN(2:nospZ)) + OM2_Ratio*SUM(ZslopN(:)) !(mmol-N/m3/d)
     OM2_PZ  = SUM(ZegP(2:nospZ)) + SUM(ZunP(2:nospZ)) + OM2_Ratio*SUM(ZslopP(:)) !(mmol-P/m3/d)
   endif

  !-------------------------------
  !-NO3; (mmol-N/m3)
  ff_new(iNO3) = NO3                            &
  &  + ( RNO3 - AupN*NO3/Ntotal)*dTd
  !--------------------------------
  !-NH4; Ammonium (mmol-N/m3)
  ff_new(iNH4) = NH4                            &
  & + ( RNH4 - AupN*NH4/(Ntotal) + AexudN + SUM(ZexN(:)) )*dTd
  !----------------------------
  !-Silica: (mmol-Si/m3)
  ff_new(iSi) =  Si                             &
  & + ( RSi - AupSi + SUM(ZegSi(:)) + SUM(ZunSi(:)) )*dTd
  !---------------------------------------------
  !-PO4: Phosphate (mmol-P/m3)
  ff_new(iPO4) = PO4                             &
  & + ( RPO4 - AupP + AexudP + SUM(ZexP(:)) )*dTd

  !---------------------------------------------------------
  !-DIC: Dissolved Inorganic Carbon (mmol-C/m3)
  ff_new(iDIC) = DIC                            &
  &  + ( RDIC - PrimProd + ArespC + ZrespC )*dTd
  !-----------------------------------------------------------------------      
  !-O2: Oxygen (mmol O2 m-3) 
  ff_new(iO2)  = O2                             &  
  &  + ( RO2  + PrimProd - ArespC - ZrespC)*dTd
  !-----------------------------------------
  !-OM1_A: (mmol-C/m3-- Dead Phytoplankton Particulate)
  ff_new(iOM1CA) = OM1CA + (ROM1CA + OM1_CA)*dTd
  ff_new(iOM1NA) = OM1NA + (ROM1NA + OM1_NA)*dTd
  ff_new(iOM1PA) = OM1PA + (ROM1PA + OM1_PA)*dTd
  if(debug.eq.2.and.is_bottom) write(6,*) "Bottom OM1A"
  if(debug.eq.2) write(6,*) "f,fnew,rOM1A,OM1C",OM1CA,ff_new(iOM1CA),ROM1CA,OM1_CA
  !-----------------------------------------------
  !-OM2_A: (mmol-C/m3-- Dead Phytoplankton Dissolved)
  ff_new(iOM2CA) = OM2CA + (ROM2CA + OM2_CA)*dTd
  ff_new(iOM2NA) = OM2NA + (ROM2NA + OM2_NA)*dTd
  ff_new(iOM2PA) = OM2PA + (ROM2PA + OM2_PA)*dTd
  if(debug.eq.2.and.is_bottom) write(6,*) "Bottom OM2A"
  if(debug.eq.2) write(6,*) "f,fnew,rOM2A,OM2C",OM2CA,ff_new(iOM2CA),ROM2CA,OM2_CA
  !-----------------------------------------------
  !-OM1_Z:(mmol-C/m3--G particulate)
  ff_new(iOM1CZ) = OM1CZ + (ROM1CZ + OM1_CZ)*dTd
  ff_new(iOM1NZ) = OM1NZ + (ROM1NZ + OM1_NZ)*dTd
  ff_new(iOM1PZ) = OM1PZ + (ROM1PZ + OM1_PZ)*dTd
  if(debug.eq.2.and.is_bottom) write(6,*) "Bottom OM1Z"
  if(debug.eq.2) write(6,*) "f,fnew,rOM1Z,OM1Z",OM1CZ,ff_new(iOM1CZ),ROM1CZ,OM1_CZ
  !-----------------------------------------------
  !-OM2_Z:(mmol-C/m3--G dissolved)
  ff_new(iOM2CZ) = OM2CZ + (ROM2CZ + OM2_CZ)*dTd
  ff_new(iOM2NZ) = OM2NZ + (ROM2NZ + OM2_NZ)*dTd
  ff_new(iOM2PZ) = OM2PZ + (ROM2PZ + OM2_PZ)*dTd
  if(debug.eq.2.and.is_bottom) write(6,*) "Bottom OM2Z"
  if(debug.eq.2) write(6,*) "f,fnew,rOM2Z,OM2Z",OM2CZ,ff_new(iOM2CZ),ROM2CZ,OM2_CZ
  !---------------------------------------------------------------------
  !-OM1_R: (mmol-C/m3--SPM particulate)
  ff_new(iOM1R) = OM1R + ROM1_R*dTd
  !---------------------------------------------------------------------
  !-OM2_R: (mmol-C/m3--SPM dissolved)
  ff_new(iOM2R) = OM2R + ROM2_R*dTd
  !---------------------------------------------------------------------
  !-OM1_BC: (mmol-C/m3--initial and boundary condition OM particulate)
  ff_new(iOM1BC) = OM1BC + ROM1_BC*dTd
  !---------------------------------------------------------------------
  !-OM2_BC: (mmol-C/m3--initial and boundary condition OM dissolved)
  ff_new(iOM2BC) = OM2BC + ROM2_BC*dTd
  !---------------------------------------------------------------------
  !-CDOM: (ppb) 
  ff_new(iCDOM) =  CDOM - CDOM*KGcdom*dTd
  !!---------------------------------------------------------------------
  !!-ALK: (mmol-HCO3/m3)
  ff_new(iALK) =  ALK +                  &
  & (RALK + AupN*NO3/(Ntotal)            &
  &          - AupN*NH4/(Ntotal)         &
  &          + AupP + 4.8*AupP)*dTd
  !Tracer
  ff_new(iTr) = ff(iTr)

  ! Before transport, Combine A/Q's
  ff_new(iQn(:)) = ff_new(iQn(:)) * ff_new(iA(:))
  ff_new(iQp(:)) = ff_new(iQp(:)) * ff_new(iA(:))

  if(is_surface.and.debug.eq.3) write(6,*) "before",ff_new(iO2),dz
  if(is_surface.and.debug.eq.3) write(6,*) "before",ff_new(iDIC)
  if(is_surface) call surface_flux(ff_new,dT,dz,T,S,Wind,pH)
  if(is_surface.and.debug.eq.3) write(6,*) "after",ff_new(iO2)
  if(is_surface.and.debug.eq.3) write(6,*) "after",ff_new(iDIC)
  return
  end subroutine cgem_step
!---------------------------------------------------------------------- 
