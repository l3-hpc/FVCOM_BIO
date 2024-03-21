#!/bin/tcsh
#BSUB -n 128 
#BSUB -W 2:00
#BSUB -R span[ptile=32]
##BSUB -x
#BSUB -q standard_ib 
#BSUB -J fvcom_mi 
#BSUB -oo out
#BSUB -eo err
module load cmaq-libs/intel2018.4-ncdf4
mpirun ./fvcom --casename=mi > mi.out
