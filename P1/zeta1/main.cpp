#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "zeta1.h"

double mpi_zeta(int n);

int main(int argc, char* argv[])
{
       //check if we got one argument
    if(argc != 2){
        std::printf("Please type one number n as argument to this program\n");
        return 1;
    }
    int n = std::stoi(argv[1]);
    //This is process 0, we need to make sure it divides the work that is to be done.
    //we need to allocate  an array with n elements of double precision.

    
    return 0;
}

double mpi_zeta(int n){
    //return calculated pi value
    /*
        We know that n is the number of elements that we need to calculate
    */
    double *pnt_array;
    pnt_array = (double*) calloc(n, sizeof(double));

    free(pnt_array);

    /*
    if(pid == 0)
        ChildProcess(); //zeta(...)
    else
        ParentProcess() //calcualte Pi and return.
    */
    double sum = 0.0;
    for(int i = 0; i < n; i++){
        sum += *(pnt_array+(i)); //itterate tru all values.
    }
    return sqrt(6*sum);
}