#include "zeta1.h"
#include <math.h> //for pow function.
#include <cstdio>

double zeta(int start, int n){ 
    double S = 0.0;
    for(int i = start; i < n+start; i++){//each process will do an equal amout of itterations.
        S += 1.0/pow(i,2);
    }
    return S;
}