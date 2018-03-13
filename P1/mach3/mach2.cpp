#include "mach2.h"
#include <math.h> //for pow function.
#include <cstdio>
#include <omp.h>

double arctan(int start, int n, double x, int threads){
    //int tid;
    double S = 0.0;

    #pragma omp parallel for num_threads(threads) reduction(+:S) schedule(static)
    for (int i = start; i < n+start; i++){
        double V = pow((-1),i-1); //1 or -1, depending on i.
        double V2 = (2*i)-1;
        double V3 = pow(x,V2);
        S += V*(V3/V2);
        //tid = omp_get_thread_num();
        //int num = omp_get_num_threads();
        //printf(" thread %d of %d working on %d \n",tid+1, num, i);
        if(i == 1){
            //printf(" num of threads = %d\n",omp_get_num_threads());
        }
    }
    return S;
}
