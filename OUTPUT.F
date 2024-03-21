       Module OUTPUT 

        SAVE

        CHARACTER(LEN=6),DIMENSION(:),allocatable :: VARIABLE_NAMES
        LOGICAL,DIMENSION(:), allocatable :: WRITE_VARIABLE
        CHARACTER(LEN=100),DIMENSION(:), allocatable :: VARIABLE_DESCRIPTIONS
        CHARACTER(LEN=100),DIMENSION(:), allocatable :: VARIABLE_UNITS
        CHARACTER(LEN=19),DIMENSION(:), allocatable :: EXTRA_VARIABLE_NAMES
        LOGICAL,DIMENSION(:), allocatable :: WRITE_EXTRA_VARIABLE
        CHARACTER(LEN=100),DIMENSION(:),allocatable :: EXTRA_VARIABLE_DESCRIPTIONS
        CHARACTER(LEN=100),DIMENSION(:),allocatable :: EXTRA_VARIABLE_UNITS
        INTEGER,DIMENSION(:),allocatable :: EXTRA_VAR ! NetCDF IDs for extra vars.
        INTEGER,DIMENSION(:),allocatable :: F_VAR ! NetCDF IDs for each variable.
        INTEGER,save :: EXTRA_VARIABLES
        INTEGER,save :: STATE_VARIABLES

       END Module OUTPUT 
