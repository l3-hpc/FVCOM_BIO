Module Model_dim

!FV use serial
!FV use mpi_interface
!FV use mpi

IMPLICIT NONE


! =========================================================
! Dimensions for Lisa's Fish Tank - USER modified part
! =========================================================
      INTEGER, SAVE :: IM, JM, KM
      INTEGER, SAVE :: NSL                   !Extra layer for cell faces
      INTEGER,ALLOCATABLE, SAVE :: NZA(:,:) 
      INTEGER, SAVE :: nospA  !Read in by main.F90 in data/Model_dim.txt
      INTEGER, SAVE :: nospZ 
      INTEGER, SAVE :: nf 
      INTEGER, SAVE :: Which_gridio
      CHARACTER(200), SAVE :: DATADIR
      INTEGER, SAVE :: iYr0 ! Reference year all timestamps are relative to.
      INTEGER, SAVE :: nRiv
      INTEGER, SAVE :: nBC
!For parallel runs
      INTEGER, SAVE :: myim !Number of cells in i direction on the processor
      INTEGER, SAVE :: myi_start !location of first array index WRT full grid
      INTEGER, SAVE :: myi_end !location of last array index WRT full grid
      INTEGER, SAVE :: myimp2  !Number of i in f array (myim + 2 boundary cells)

CONTAINS

Subroutine Set_Model_dim(myid,numprocs)


IMPLICIT NONE


      integer, intent(in) :: myid,numprocs
      character(200) filename
      integer mpierr

if(myid.eq.0) then

! --- Read in location of data:
      open(unit=19,file="./data/MyFiles.inp", form="formatted", status="old")
      read(19,'(a)') DATADIR
      close(19)
!"
      write(filename,'(A, A)') trim(DATADIR),'/Model_dim.txt'


! --- Read in Model_dim parameters
      open(unit=19,file=filename, form="formatted", status="old")
      read(19,*) !Header
      read(19,*) im
      read(19,*) jm
      read(19,*) km 
      read(19,*) nospA
      read(19,*) nospZ
      read(19,*) Which_gridio
      read(19,*) iYr0
      read(19,*) nRiv
      read(19,*) nBC
      close(19)
      
      if(Which_gridio.eq.2 .or. Which_gridio.eq.3) then
         nsl = km + 1 
      else
         nsl = km
      endif

endif

!FV if(numprocs.gt.1) then
!FV       call MPI_BCAST(im,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
!FV       call MPI_BCAST(jm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(km,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(nospA,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(nospZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(Which_gridio,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(iYr0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(DATADIR,200,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(nsl,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(nRiv,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV       call MPI_BCAST(nBC,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
!FV endif

!Divide grid among processors
!FV      call Decomp1D(im, numprocs, myid, myi_start, myi_end)

!Define Loop Bounds
      myim = myi_end - myi_start + 1
      myimp2 = myim + 2

      call Model_dim_allocate()

#ifdef DEBUG
write(6,*) "---Set_Model_dim----"
write(6,*) "Model_dim.txt=",filename
write(6,*) "im,jm,km,nsl",im,jm,km,nsl
write(6,*) "im,myim",im,myim
write(6,*) "myi_start,end",myi_start,myi_end
write(6,*)
#endif

Return

End Subroutine


Subroutine Model_dim_allocate()

!FV USE Fill_Value

IMPLICIT NONE

integer ierr

ALLOCATE(nza(im,jm),stat=ierr)
!FV if(ierr.ne.0) call error("nza",ierr)

!Initialize for netCDF
!FV nza=fill(0)
nza=9999.

#ifdef map_code
write(6,*) "---Model_dim_allocate---"
!FV write(6,*) "allocating nza, fill with",fill(0)
write(6,*)
#endif

Return

END Subroutine Model_dim_allocate


END Module Model_dim
