MODULE MOD_1D
!******************************************************************
!***        DEFINE all PARAMETERS needed for every bio          ***
!***   PARAMETERS ARE LISTED ALPHABETICALLY FOR EACH CATEGORY   ***
!******************************************************************
#    if !defined (DOUBLE_PRECISION)
     INTEGER, PARAMETER    :: SPP = SELECTED_REAL_KIND(6,30)
#    else
     INTEGER, PARAMETER    :: SPP = SELECTED_REAL_KIND(12,300)
#    endif
     INTEGER, PARAMETER    :: DPP = SELECTED_REAL_KIND(12,300)
     INTEGER               :: KBV              ! VERTICAL LAYERS
     INTEGER               :: KBVM1            ! VERTICAL LAYERS MINUS 1
     INTEGER               :: NTT              ! TOTAL NUMBER OF STATE VARIALBES
     !Number of iterations
     INTEGER               :: ITS
     !Node number, for debugging
     INTEGER               :: INEA


     CHARACTER(LEN=50)               :: STRUCTURE     ! MODEL STRUCTURE
     CHARACTER(LEN=50), ALLOCATABLE  :: BIO_NAME(:,:) ! VARIABLE NAME FOR NETCDF OUTPUT

!     REAL(SPP), PARAMETER  :: ZERO             = 0.0_SPP
     REAL(SPP), ALLOCATABLE :: BIO_VAR(:,:)     ! ENSEMBLE VARIABLES FOR MIXING
     REAL(SPP), ALLOCATABLE :: DELTA_D(:)       ! LAYER THICKNESS
     REAL(SPP), ALLOCATABLE :: DELTA_Z(:)       ! LAYER BETWEEN SIGMA CENTRE
     REAL(SPP), ALLOCATABLE :: DEPTH_Z(:)       ! DEPTH AT EACH LAYER CENTRE
     REAL(SPP), ALLOCATABLE :: KM_BIO(:)        ! VERTICAL MIXING COEFFICIENT
     REAL(SPP), ALLOCATABLE :: RADIATION(:)     ! SURFACE RADIATION
     REAL(SPP), ALLOCATABLE :: T_BIO(:)         ! TEMPERATURE
     REAL(SPP), ALLOCATABLE :: S_BIO(:)         ! Salinity
     REAL(SPP)              :: T_STEP           ! TIME STEP
END MODULE MOD_1D
