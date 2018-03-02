#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <mpi.h>

#include "zeta1.h"

double mpi_zeta();

//Global variables
int numprocs, n;
int myid;
double time_start;

int main(int argc, char* argv[])
{   
    MPI_Status status;

    if(argc != 2){
        std::printf("Please type one number n as argument to this program and\n");
        return 1;
    }
    n = std::stoi(argv[1]);

    MPI_Init(&argc, &argv); //init MPI

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//Get the number of processors.
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); //Get my rank(id)

    //Check if the number of processes are not a power of 2.(1,2,..,2^n)
    if( (numprocs & (numprocs - 1)) != 0) { 
        printf("Number of processes are not a number of 2\n");
        MPI_Finalize();
        return 0;
    }

    if(myid == 0){ //We are master.
        time_start =  MPI_Wtime(); //Initialize a time, to measure the duration of the processing time.
    }else{ //We are a slave
        ; //nothing to do
    }

    mpi_zeta();

    MPI_Finalize();
    return 0;
}

double mpi_zeta(){
    int divide = (double)n/numprocs; //compute the minimum amount of work for each process. 
    double sum = 0.0;
    int rest_tasks = n%numprocs;
    if(rest_tasks != 0){ //If we can't cleanly divide the number of tasks between the processes
        //We need to check how many extra tasks there are.
        //We simply need to give rest_task processes an extra amout of work to complete
        int tasks = n; //We divide this out to the tasks.

        if(rest_tasks > myid){ //we want the tasks with id higher than the work to do the extra task.
            int start = myid*divide + myid + 1;
            int n = divide + 1;
            sum = zeta(start, n, myid); //each process calculate it's sum.
        }else
        {   //We need to shift the work intervals by the rest_tasks.
            int start = myid*divide + 1 + rest_tasks;
            int n = divide;
            sum = zeta(start, n, myid); //each process calculate it's sum.
        }
    }else{//it divides cleanly between processes.
        
        sum = zeta(myid*divide + 1, divide, myid); //each process calculate it's sum.
    }

    double sum_all = 0.0;

    MPI_Reduce(&sum, &sum_all, 1 , MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); //we sum the sum variable from each process and store in pi.

    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        double pi = sqrt(sum_all*6); 
        double error = fabs(pi - 4.0 * atan(1.0));
        printf("pi = %e, error=%e, duration=%e ms\n",pi, error, duration*1000);
    }
}