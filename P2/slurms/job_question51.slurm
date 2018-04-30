#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem=50000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="poisson_test"
#SBATCH --output=output/question51.out
#SBATCH --mail-user=sondenge@ntnu.no
#SBATCH --mail-type=ALL

module load GCC OpenMPI
#CC=mpicc FC=mpif90
cd $PBS_O_WORKDIR
n=$((2**12))
for i in `seq 1 3 36`
do
    mpiexec -np 12 git/TMA4280v2018/P2/poisson_paralell/poisson $n 1 5 5 #simply request rhs = two points , -1 and 1. 
done

echo 'jobs done'
