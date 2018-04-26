#!/bin/bash
for i in `seq 0 10`
do
    #echo $[2**$i] #do the power of two series
    mpirun -np 4 poisson $[2**$i] 2
done