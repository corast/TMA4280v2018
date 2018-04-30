#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem=50000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="poisson_test"
#SBATCH --output=output/question3_new_pt.out
#SBATCH --mail-user=sondenge@ntnu.no
#SBATCH --mail-type=ALL

module load GCC OpenMPI
#CC=mpicc FC=mpif90
cd $PBS_O_WORKDIR
#mpiexec git/TMA4280v2018/P2/poisson_paralell/poisson 1
arrayp=(1 2 3 4 6 9 12 18 36)
arrayt=(36 18 12 9 6 4 3 2 1)
len=${#arrayp[@]}
n=$((2**12))
#do the different combinations of p*nt(per process) = 36

echo "Testing run time on pt = 36"
for ((i=0; i<${len}; i++));
do
    mpiexec -np ${arrayp[$i]} git/TMA4280v2018/P2/poisson_paralell/poisson $n ${arrayt[$i]} 0 
done

echo 'jobs done'
