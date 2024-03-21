Subroutine Read_InputFile_GD(filename,myid,numprocs)

USE Model_dim
USE INPUT_VARS
USE INPUT_VARS_GD
USE STATES
USE EUT
USE FLAGS
USE LIGHT_VARS
USE DATE_TIME
!FV use mpi_interface

IMPLICIT NONE

integer icent_jcent_units
character(120), intent(in) :: filename
integer, intent(in) :: myid, numprocs
!FV integer mpierr
real icent_in, jcent_in

ws = 0.

if(myid.eq.0) then
!--Code Identifier--------------
open(unit=999,file=filename,form='formatted',status='old')
read(999,*) code_ID
read(999,*)
!--Simulation specifics------
read(999,*)
read(999,*) iYrS,iMonS,iDayS,iHrS,iMinS,iSecS
read(999,*) iYrE,iMonE,iDayE,iHrE,iMinE,iSecE
read(999,*) dT, dT_out
read(999,*) icent_jcent_units
read(999,*) icent_in, jcent_in
icent = INT(icent_in)
jcent = INT(jcent_in)

read(999,*)
!--Switches in GEM---------
read(999,*)
read(999,*) Which_Fluxes
read(999,*) !Comment Line
read(999,*) Read_Solar,Read_Wind,Read_T,Read_Sal
read(999,*) InitializeHow
read(999,*) 
!--Optics-----------------------
read(999,*)
read(999,*) Which_irradiance
read(999,*) astar490
read(999,*) aw490
read(999,*) astarOMA
read(999,*) astarOMZ
read(999,*) astarOMR
read(999,*) astarOMBC
read(999,*) PARfac
read(999,*)
!----River Params-----------------------------------
read(999,*) 
read(999,*) rcNO3
read(999,*) rcNH4
read(999,*) rcPO4
read(999,*) rcSi
read(999,*)
!----Other including Boundary Conditions-------------
read(999,*)
read(999,*) Which_Vmix, Which_Adv   
read(999,*) KH_coeff  
read(999,*) Which_Outer_BC 
read(999,*) m_OM_init,m_OM_bc,m_OM_sh 
read(999,*) DoDroop
!---------------------------------------------------

!--GOMDOM params------------------------------------------------
read(999,*)
read(999,*)             !C
read(999,*)             !C  FIREAD file
read(999,*)             !C
read(999,*) 
read(999,*) ALPHA_DIA     ! ALPHA_DIA: Initial slope of diatom light sat curve 
read(999,*) ALPHA_GRE     ! ALPHA_GRE: Initial slope of greens light sat curve
read(999,*) 
read(999,*) ANCP          ! ANCP: N:C ratio
read(999,*) APCP          ! APCP: P:C ratio
read(999,*) ASCD          ! ASCD: Si:C ratio
read(999,*) AVFRAC        ! AVFRAC: Available DOP fraction
read(999,*) AVFRACDON     ! AVFRACDON: Available DON fraction
read(999,*) 
read(999,*) BMRD          ! BMRD: diatom mortality
read(999,*) BMRG          ! BMRG: greens mortality
read(999,*) 
read(999,*) CCHLD         ! CCHLD: C:chl ratio for diatoms
read(999,*) CCHLG         ! CCHLG: C:chl ratio for greens
read(999,*) 
read(999,*) CGZ           ! CGZ: zoo grazing rate
read(999,*)
read(999,*) DENIT_CN_RATIO ! DENIT_CN_RATIO 
read(999,*)
read(999,*) FCDD          ! FCDD: Frac mort = DOC (d)
read(999,*) FCDG          ! FCDG: Frac mort = DOC (g)
read(999,*) FCDP          ! FCDP: Frac pred = DOC
read(999,*) FCDZ          ! FCDZ: Frac z-mort = DOC
read(999,*) FCLD          ! FCLD: Frac mort = LOC (d)
read(999,*) FCLG          ! FCLG: Frac mort = LOC (g)
read(999,*) FCLP          ! FCLP: Frac LPDC (pred)
read(999,*) FCLZ          ! FCLZ: Frac LPDC z-mort
read(999,*) FCRD          ! FCRD: Frac mort = ROC (d)
read(999,*) FCRG          ! FCRG: Frac mort = ROC (g)
read(999,*) FCRP          ! FCRP: Frac RPDC (pred)
read(999,*) FCRZ          ! FCRZ: Frac RPDC z-mort
read(999,*) FNDD          ! FNDD: Frac mort = DON(d)
read(999,*) FNDG          ! FNDG: Frac mort = DON (g)
read(999,*) FNDP          ! FNDP: Frac pred= DON
read(999,*) FNDZ          ! FNDZ: Frac z-mort = DON
read(999,*) FNID          ! FNID: Frac mort=DIN (d)
read(999,*) FNIG          ! FNIG: Frac mort = DIN(g)
read(999,*) FNIP          ! FNIP: Frac pred = DIN
read(999,*) FNIZ          ! FNIZ: Frac DIN z-mort
read(999,*) FNLD          ! FNLD: Frac mort=LON(d)
read(999,*) FNLG          ! FNLG: Frac mort= LON(g)
read(999,*) FNLP          ! FNLP: Frac pred=LON
read(999,*) FNLZ          ! FNLZ: Frac z-mort=LON
read(999,*) FNRD          ! FNRD: Frac mort=RON(d)
read(999,*) FNRG          ! FNRG: Frac mort=RON(g)
read(999,*) FNRP          ! FNRP: Frac pred=RON
read(999,*) FNRZ          ! FNRZ: Frac z-mort= RON
read(999,*) FPDD          ! FPDD: Frac mort=DOP(d)
read(999,*) FPDG          ! FPDG: Frac mort = DOP(g)
read(999,*) FPDP          ! FPDP: Frac pred=DOP
read(999,*) FPDZ          ! FPDZ: Frac z-mort DOP
read(999,*) FPID          ! FPID: Frac mort=DIP(d)
read(999,*) FPIG          ! FPIG: Frac mort=DIP(g)
read(999,*) FPIP          ! FPIP: Frac pred=DIP
read(999,*) FPIZ          ! FPIZ: Frac z-mort=DIP
read(999,*) FPLD          ! FPLD: Frac mort=LOP(d)
read(999,*) FPLG          ! FPLG: Frac mort= LOP(g)
read(999,*) FPLP          ! FPLP: Fract pred=LOP
read(999,*) FPLZ          ! FPLZ: Frac z-mort LOP
read(999,*) FPRD          ! FPRD: Frac mort=ROP(d)
read(999,*) FPRG          ! FPRG: Frac mort=ROP(g)
read(999,*) FPRP          ! FPRP: Frac pred=ROP
read(999,*) FPRZ          ! FPRZ: Frac z-mort ROP
read(999,*) FSAP          ! FSAP: Frac pred = DSi
read(999,*) 
read(999,*) GREFF         ! GREFF: zoo grazing coeff
read(999,*) 
read(999,*) ILMUL         ! ILMUL: scaling factor for surface short wave radiation
read(999,*) 
read(999,*) KDC           ! KDC: min DOC mineraliz
read(999,*) KDCALG        ! KDCALG: DOC alg dependance
read(999,*) KDN           ! KDN: min DON mineraliz
read(999,*) KDNALG        ! KDNALG: DON alg depedence
read(999,*) KDP           ! KDP: min DOP mineral
read(999,*) KDPALG        ! KDPALG: DOP alg dependence
read(999,*) 
read(999,*) KDWD          ! KDWD: Specifies light attenuation eqn, default=0
read(999,*) KE            ! KE: Background light attenuation [m^-1]
read(999,*) KECHL         ! KECHL: Light attenuation factor for chlorophyll-a [m^2-kg^-1]
read(999,*) 
read(999,*) KHDONT_SED    ! KHDONT_SED
read(999,*) KHN           ! KHN: Organic N Decay Half Saturation Constant
read(999,*) KHND          ! KHND: mean N half sat (dia)
read(999,*) KHNG          ! KHNG: mean N half sat (gre)
read(999,*) KHNNT         ! KHNNT: half sat of NH3 for nitr
read(999,*) KHODOC_SED    ! KHODOC_SED
read(999,*) KHP           ! KHP: Organic P Decay Half Saturation Constant
read(999,*) KHPD          ! KHPD: mean P half sat (dia)
read(999,*) KHPG          ! KHPG: mean P half sat (gre)
read(999,*) KHSD          ! KHSD: mean Si half sat (dia)
read(999,*) KLC           ! KLC: min hyd for LOC
read(999,*) KLCALG        ! KLCALG: LOC alg dependence
read(999,*) KLN           ! KLN: min hyd for LON
read(999,*) KLNALG        ! KLNALG: LON alg dependence
read(999,*) KLP           ! KLP: min hyd for LOP
read(999,*) KLPALG        ! KLPALG: LOP alg dependence
read(999,*) KRC           ! KRC: min hyd for ROC
read(999,*) KRCALG        ! KRCALG: POC alg dependence
read(999,*) KRN           ! KRN: min hyd for RON
read(999,*) KRNALG        ! KRNALG: RON alg dependence
read(999,*) KRP           ! KRP: min hyd for ROP
read(999,*) KRPALG        ! KRPALG: ROP alg dependence
read(999,*) KSUA          ! KSUA: BSi diss rate
read(999,*) KSZ           ! KSZ: zoo half sat phyto
read(999,*) KTBD          ! KTBD: dia mort temp coeff 
read(999,*) KTBG          ! KTBG: gre mort temp coeff 
read(999,*) KTGD1         ! KTGD1: temp coeff < (dia)
read(999,*) KTGD2         ! KTGD2: temp coeff > (dia)
read(999,*) KTGG1         ! KTGG1: temp coeff < (gre)
read(999,*) KTGG2         ! KTGG2: temp coeff > (gre)
read(999,*) KTHDR         ! KTHDR: hyd temp dep
read(999,*) KTMNL         ! KTMNL: min temp dep
read(999,*) KTNT1         ! KTNT1: temp coeff < nitr
read(999,*) KTNT2         ! KTNT2: temp coeff > nitr
read(999,*) KTSUA         ! KTSUA: silic diss temp coeff
read(999,*) NTM           ! NTM: nitrif rate 
read(999,*) PBMAX_DIA     ! PBMAX_DIA: photo rate at opt illum
read(999,*) PBMAX_GRE     ! PBMAX_GRE: photo rate at opt illum
read(999,*) PMD           ! PMD: dia prod
read(999,*) PMG           ! PMG: gre prod
read(999,*) 
read(999,*) RSODNTR       ! RSODNTR
read(999,*)
read(999,*) SILIM         ! SILIM: Minimum nutrient limits growth
read(999,*) 
read(999,*) TMD           ! TMD: opt temp dia
read(999,*) TMG           ! TMG: opt temp gre
read(999,*) TMNT          ! TMNT: opt temp nitr
read(999,*) TRD           ! TRD: opt metab temp (dia)
read(999,*) TRG           ! TRG: opt metab temp (gre)
read(999,*) TRHDR         ! TRHDR: opt hydr temp
read(999,*) TRMNL         ! TRMNL: opt mineral temp
read(999,*) TRSUA         ! TRSUA: opt Si diss temp
read(999,*) TZREF         ! TZREF: opt pred temp
read(999,*) TREF          ! TREF: SOD reference temp
read(999,*) ZDTH          ! ZDTH: zoo mort
read(999,*) ZTHET         ! ZTHET: temp coeff pred
read(999,*) 
read(999,*) !  Dissolved Oxygen related parameters
read(999,*) 
read(999,*) KCOD            ! KCOD: DOC oxidation rate
read(999,*) KDENITR         ! KDENITR: Max denitrification rate coefficient
read(999,*) KHDENITR        ! KHDENITR: Half-sat [NO3] for dentrification
read(999,*) KHDONT          ! KHDONT: Half-sat [D.O] required for nitrification
read(999,*) KHOCOD          ! KHOCOD: Half-sat [D.O] required for exertion of COD
read(999,*) KHODOC          ! KHODOC: Half-sat [D.O] required for oxic respiration
read(999,*) KRDO            ! KRDO: D.O. reaeration coefficient (m/s)
read(999,*) RCDO            ! RCDO: D.O. to carbon ratio in resp 
read(999,*) RNTO            ! RNTO: ratio of O2 consumed to nitrif
read(999,*) 
read(999,*) !Droop Kinetics
read(999,*)
read(999,*) FINTNID
read(999,*) FINTNDD
read(999,*) FINTNLD
read(999,*) FINTNRD
read(999,*)
read(999,*) FINTNIG
read(999,*) FINTNDG
read(999,*) FINTNLG
read(999,*) FINTNRG
read(999,*)
read(999,*) FINTLUXNIP
read(999,*) FINTSTRNIP
read(999,*) FINTLUXNDP
read(999,*) FINTSTRNDP
read(999,*) FINTLUXNLP
read(999,*) FINTSTRNLP
read(999,*) FINTLUXNRP
read(999,*) FINTSTRNRP
read(999,*)
read(999,*) FINTPID
read(999,*) FINTPDD
read(999,*) FINTPLD
read(999,*) FINTPRD
read(999,*)
read(999,*) FINTPIG
read(999,*) FINTPDG
read(999,*) FINTPLG
read(999,*) FINTPRG
read(999,*)
read(999,*) FINTLUXPIP
read(999,*) FINTSTRPIP
read(999,*) FINTLUXPDP
read(999,*) FINTSTRPDP
read(999,*) FINTLUXPLP
read(999,*) FINTSTRPLP
read(999,*) FINTLUXPRP
read(999,*) FINTSTRPRP
read(999,*)
read(999,*) KHINTND
read(999,*) KHINTNG
read(999,*) KHINTPD
read(999,*) KHINTPG
read(999,*)
read(999,*) QMINND
read(999,*) QMINNG
read(999,*) QMINPD
read(999,*) QMINPG
read(999,*) QMAXND
read(999,*) QMAXNG
read(999,*) QMAXPD
read(999,*) QMAXPG
read(999,*) UPNMAXD
read(999,*) UPNMAXG
read(999,*) UPPMAXD
read(999,*) UPPMAXG
read(999,*)
read(999,*) ! Settling rates
read(999,*) 
read(999,*) ws(JDIA)  ! VDIA
read(999,*) ws(JDIAN)  ! VDIAN
read(999,*) ws(JDIAP)  ! VDIAP
read(999,*) ws(JGRE)  ! VGRE	
read(999,*) ws(JGREN)  ! VGREN
read(999,*) ws(JGREP)  ! VGREP
read(999,*) ws(JLOC)  ! VLOC
read(999,*) ws(JROC)  ! VROC
read(999,*) ws(JLON)  ! VLON
read(999,*) ws(JRON)  ! VRON
read(999,*) ws(JLOP)  ! VLOP
read(999,*) ws(JROP)  ! VROP
read(999,*) ws(JSU)   ! VSU
read(999,*) ws(JTR)   ! VTR
endif
!FV if(numprocs.gt.1) then
!FV call MPI_BCAST(code_ID,50,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
!FV call MPI_BCAST(iDayE,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV call MPI_BCAST(iDayS,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV
!FV ...more MPI_BCAST... 
!FV
!FV endif
!FV
    ! Compute starting time of run in seconds since Model_dim::iYr0:
      START_SECONDS = &
      TOTAL_SECONDS( iYr0, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS )

      END_SECONDS = &
      TOTAL_SECONDS( iYr0, iYrE, iMonE, iDayE, iHrE, iMinE, iSecE )

#ifdef DEBUG
      write(6,*) "iYr0",iYr0, iYrS, iMonS, iDayS, iHrS, iMinS, iSecS
      write(6,*) "START_SECONDS",START_SECONDS,START_SECONDS/3600./24./365.
#endif

!FV not defined, causes integer devide by zero 
       nstep = ( END_SECONDS - START_SECONDS ) / real(dT) !number of timesteps in a run
       iout = dT_out/real(dT) !output time-interval in timesteps

! --- sinking speed: converted from m/s downward positive to m/s negative
      ws = -ws
      if(Which_Adv.eq.2.or.Which_adv.eq.3) ws = 0.
      write(6,*) "Which_Adv=",Which_Adv
      write(6,*) "Which_VMix=",Which_VMix

      return

END Subroutine Read_InputFile_GD
