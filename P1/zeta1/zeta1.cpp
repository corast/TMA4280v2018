#include "zeta1.h"
#include <math.h> //for pow function.
#include <cstdio>

double zeta(int start, int n, int myid){ 
    double S = 0.0;
    printf("Process %d calculate interval: [%d , %d]\n",myid, start, start+n-1);
    for(int i = start; i < n+start; i++){//each process will do an equal amout of itterations.
        S += 1.0/pow(i,2);
    }
    return S;
}