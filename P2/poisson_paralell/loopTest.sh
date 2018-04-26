#!/bin/bash
for i in `seq 0 8`
do
    #echo $[2**$i] #do the power of two series
    mpirun -np 3 poisson $[2**$i] 1
done