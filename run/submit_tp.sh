#!/bin/bash
#SBATCH --job-name="letp"
#SBATCH --output="qletp-out.%j"
#SBATCH --error="qletp-err.%j"
##SBATCH --partition=shared
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
##SBATCH --ntasks-per-node=1
#SBATCH --account=ncs126
#SBATCH -t 12:00:00
module load intel/19.1.3.304/6pv46so
module load intel-mpi/2019.10.317/ezrfjne
module load netcdf-c/4.8.1/ej3w5cy
module load netcdf-fortran/4.5.3/vqnicf7
#1node
#/usr/bin/time -v mpirun -n 128 ./fvcom_tp --casename=letp > letp.out
#2nodes
/usr/bin/time -v mpirun -n 256 ./fvcom_tp --casename=letp > letp.out
#mpirun -n 1 ./fvcom_tp --casename=letp > letp.out
echo "seff"
seff $SLURM_JOB_ID 
echo "sacct"
sacct --jobs=$SLURM_JOB_ID --format=nnodes,ncpu,ntasks,maxrss,cputime,avecpu,maxdiskread,maxdiskwrite,maxpages
