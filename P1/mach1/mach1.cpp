#include "mach1.h"
#include <math.h> //for pow function.
#include <cstdio>

double arctan(int start, int n, double x, int myid);

double mach(int start, int n, int myid){ 
    double S = 0.0;
    printf("Process %d calculate interval: [%d , %d]\n",myid, start, start+n-1);
    for(int i = start; i < n+start; i++){//each process will do an equal amout of itterations.
        S += 1.0/pow(i,2);
    }
    return S;
}

double arctan(int start, int n, double x, int myid){
    printf("Process %d calculate interval: [%d , %d]\n",myid, start, start+n-1);
    double S = 0.0;
    for (int i = start; i < n+start; i++){
        double V = pow((-1),i-1);
        double V2 = (2*i)-1;
        double V3 = pow(x,V2);
        S += V*(V3/V2);
    }
    return S;
}
