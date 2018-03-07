#include "mach1.h"
#include <math.h> //for pow function.
#include <cstdio>

double arctan(int start, int n, double x, int myid){
    double S = 0.0;
    for (int i = start; i < n+start; i++){
        double V = pow((-1),i-1);
        double V2 = (2*i)-1;
        double V3 = pow(x,V2);
        S += V*(V3/V2);
    }
    //printf("process %i S= %d",);
    return S;
}
