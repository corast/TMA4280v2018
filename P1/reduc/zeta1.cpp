#include "zeta1.h"
#include <math.h> //for pow function.
#include <cstdio>

double zeta(int start, int n, int myid){ 
    double S = 0.0;
    int end_interval = n == 0 ? 0 : start+n-1; //correctly label the end_interval in printout.
    printf("Process %d calculate interval: [%d , %d]\n",myid, start, end_interval);
    for(int i = start; i < n+start; i++){//each process will do an equal amout of itterations.
        S += 1.0/pow(i,2);
    }
    return S;
}