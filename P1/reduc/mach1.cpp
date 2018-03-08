#include "mach1.h"
#include <math.h> //for pow function.
#include <cstdio>

double arctan(int start, int n, double x){
    double S = 0.0;
    for (int i = start; i < n+start; i++){
        double V = pow((-1),i-1); //1 or -1, depending on i.
        double V2 = (2*i)-1;
        double V3 = pow(x,V2);
        S += V*(V3/V2);
    }
    return S;
}
