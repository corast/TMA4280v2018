#include "zeta2.h"
#include <math.h> //for pow function.
#include <cstdio>
#include <omp.h>

double zeta(int n, int threads){ 
    double S = 0.0;
    #pragma omp paralell for num_threads(threads) schedule(static) reduction(+:S)
    for(int i = 0; i < n; i++){//each process will do an equal amout of itterations.
        //omp_get_thread_num();
        S += 1.0/pow(i,2);
    }

    printf("S %f ",S);
    return S;
}