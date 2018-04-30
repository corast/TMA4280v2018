#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:15:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem=12000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="poisson_test"
#SBATCH --output=output/question1_delay.out
#SBATCH --mail-user=sondenge@ntnu.no
#SBATCH --mail-type=ALL

module load GCC OpenMPI
#CC=mpicc FC=mpif90
cd $PBS_O_WORKDIR
#mpiexec git/TMA4280v2018/P2/poisson_paralell/poisson 1

echo "Time it takes to divide the work between the different process on process 0"
for i in `seq 1 18`
do 
    mpiexec -np $i git/TMA4280v2018/P2/poisson_paralell/poisson $[2**12] 1 0
done

echo 'jobs done'
