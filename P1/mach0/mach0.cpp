#include "mach0.h"
#include <math.h>

double arctan(int n, double x);

double mach(int n){

    double arctan_1 = arctan(n, (double)1/5);
    double arctan_2 = arctan(n, (double)1/239);

    double pi = 4*(4*arctan_1 - arctan_2);

    return pi;
}


double arctan(int n, double x){
    double S = 0.0;
    for(int i = 1; i <= n; i++){
        double V = pow((-1),i-1);
        double V2 = (2*i)-1;
        double V3 = pow(x,V2);
        S += V*(V3/V2);
    }
    return S;
}

