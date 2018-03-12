#include "zeta2.h"
#include <math.h> //for pow function.
#include <cstdio>
#include <omp.h>

double zeta(int n, int threads){ 
    double S = 0.0;
    int tid;
    printf("Start running in paralell\n");
    #pragma omp paralell for num_threads(threads) reduction(+:S) schedule(static)
    for(int i = 1; i <= n; i++){//each process will do an equal amout of itterations.
        tid = omp_get_thread_num();
        S +=1.0 /pow(i,2);
        printf(" thread %d working on %d \n",tid, i);
    }
    return S;
}