#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <mpi.h>

#include "zeta1.h"

double mpi_zeta(int n, double* pnt_array);

int main(int argc, char* argv[])
{   
    int numprocs, n;
    int myid;
    MPI_Status status;

    MPI_Init(&argc, &argv); //init MPI

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//Get the number of processors.
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); //Get my rank(id)

    //Only want to run this on master once.
    if(myid == 0){
        //check if we got one argument
        if(argc != 2){
            std::printf("Please type one number n as argument to this program and\n");
            return 1;
        }
        n = std::stoi(argv[1]);
    }

    //Check if the number of processes are not a power of 2.(1,2,..,2^n)
    if( (numprocs & (numprocs - 1)) != 0 || numprocs == 1) { 
        std::printf("Number of processes are not a number of 2\n");
        MPI_Finalize();
        return 0;
    }
    
    /*
        Main continues from here as process 0, whist the rest start from the top in main.
    */
    if(myid == 0){
        //We are master.
        printf("We have %d processes\n", numprocs);
        double *pnt_array;
        pnt_array = (double*) calloc(n, sizeof(double)); //Alloate an array with enough space for n double elements.
    }

    //double pi = mpi_zeta(n);

    MPI_Finalize();
    return 0;
}

double mpi_zeta(int n, double* pnt_array){
    //return calculated pi value
    /*
        We know that n is the number of elements that we need to calculate
    */
    


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
    free(pnt_array); //Free the memory, since we are done using it. We don't want memory leaks.
    return sqrt(6*sum);
}