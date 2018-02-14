#include "zeta0.h"
#include <math.h>

double zeta0(int n)
{
    double S = 0.0;

    for(int i = 1; i <= n; i++){
        S += 1/(i*i);
    }

    //S_n = pi^2/6; 
    return sqrt(S*6);
}