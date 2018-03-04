#include "mach1.h"
#include <math.h> //for pow function.
#include <cstdio>

double arctan(int start, int n, double x, int myid){
    int end_interval = n == 0 ? 0 : start+n-1; //correctly label the end_interval in printout.
    printf("x = %f Process %d calculate interval: [%d , %d]\n",x, myid, start, end_interval);
    double S = 0.0;
    for (int i = start; i < n+start; i++){
        double V = pow((-1),i-1);
        double V2 = (2*i)-1;
        double V3 = pow(x,V2);
        S += V*(V3/V2);
    }
    return S;
}
