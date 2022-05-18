#!/bin/bash

#SBATCH --partition=nufgroup
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=72:00:00
#SBATCH --oversubscribe

#SBATCH --job-name=ib_parallel
#SBATCH --output=status.out

##module load petsc/3.11.1-i
module load mvapich2/2.3-intel-18.0
module load fftw-3.3.8/intel-18.0/mvapich2-2.3

# Change to the directory from which the batch job was submitted
cd $SLURM_SUBMIT_DIR

mpirun -np 4 bin/ib > ib.out


