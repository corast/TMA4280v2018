#!/bin/bash
for i in `seq 2 10`
do
    #echo $[2**$i] #do the power of two series
    mpiexec -np 4 poisson $[2**$i] 1
done
echo "jobs done"