#include "zeta1.h"
#include <math.h> //for pow function.

/*  Distributed calculation of vi, where each v value, is to be stored on an array in memory
    n is the number of elements this process should calculate*/

void zeta(int s_val, int n, double *ptr_array){ //We require an pointer to an array, and 
    //compute the partial sum of zeta given n.
    for(int i = s_val; i < n+s_val; i++){
        *(ptr_array) = 1/pow((double)i, 2); //Take the value from the array positon and calculate v_i. 
    }
}