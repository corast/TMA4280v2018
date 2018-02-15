#include "zeta0.h"
#include <math.h> //for sqrt function

double zeta0(int n)
{
    double S = 0.0;

    for(int i = 1; i <= n; i++){
        S += (double)1/(i*i);
    }
    //Now that we have calculated S, we need to multiply by 6, and find the square root.
    //S_n = pi^2/6; 
    return sqrt(S*6);
}