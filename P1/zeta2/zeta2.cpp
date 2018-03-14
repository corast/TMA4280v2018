#include "zeta2.h"
#include <math.h> //for pow function.
#include <cstdio>
#include <omp.h>

double zeta(int n, int threads){ 
    double S = 0.0;
<<<<<<< HEAD
    //int tid;
    //omp_set_dynamic(0);
    //omp_set_num_threads(4);
    //printf(" max threads %d\n",omp_get_max_threads());
    #pragma omp parallel for num_threads(threads) reduction(+:S) schedule(static)
    for(int i = 1; i <= n; i++){//each process will do an equal amout of itterations.
        //tid = omp_get_thread_num();
        if(i == 1){
            printf(" num of threads = %\n",omp_get_num_threads());
        }
        
        S +=1.0 /pow(i,2);
        //printf(" thread %d of %d working on %d \n",tid+1, num, i);
=======
    int tid;
    printf("RUnning");
    //omp_set_dynamic(0);
    //omp_set_num_threads(8);
    printf(" max threads %d\n",omp_get_max_threads());
    #pragma omp parallel for num_threads(4) schedule(static)
    for(int i = 1; i <= n; i++){//each process will do an equal amout of itterations.
        tid = omp_get_thread_num();
        int num = omp_get_num_threads();
        printf(" thread %d of %d working on %d \n",tid+1, num, i);
        S +=1.0 /pow(i,2);
>>>>>>> bf84f160e56e857a05d4f480ef2598c1b78c06f2
    }

    return S;
}