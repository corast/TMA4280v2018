#include "mach0.h"
#include <math.h>

double mach(int n, int x){

    double S = 0.0;

    for(int i = 1; i <= n; i++){
        double V = pow((-1),i-1);
        double V2 = (2*i)-1;
        double V3 = pow(x,V2);
        S += V*(V3/V2);
    }

    return S;
}

