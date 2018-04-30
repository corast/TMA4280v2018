#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem=50000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="poisson_test"
#SBATCH --output=output/question4_p10.csv
#SBATCH --mail-user=sondenge@ntnu.no
#SBATCH --mail-type=ALL

module load GCC OpenMPI
#CC=mpicc FC=mpif90
cd $PBS_O_WORKDIR
#create an file, which we can read like an csv file.
n=$((2**10))
echo "np;time"
for i in `seq 1 64`
do
    mpiexec -np $i git/TMA4280v2018/P2/poisson_paralell/poisson $n 1 41
done
