#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <tuple>
#include <mpi.h>

#include "zeta2.h"

void mpi_zeta();
std::tuple<int, int> calculate_work();

//Global variables
int numprocs, n, threads;
int myid = -1; //Default value.
std::clock_t start;
double time_start;

int main(int argc, char* argv[])
{   
    //We need to pass more arguments.
    if(argc != 3){
        std::printf("Please type number of iteration, number of process and number of threads \n");
        return 1;
    }
    n = std::stoi(argv[1]);//Number of itterations.
    threads = std::stoi(argv[2]);//Number of threads
    
    MPI_Init(&argc, &argv); //init MPI

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//Get the number of processors.
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); //Get my rank(id)

    //Check if the number of processes are not a power of 2.(1,2,..,2^n)
    if( (numprocs & (numprocs - 1)) != 0 && numprocs) { 
        printf("Number of processes are not a number of 2\n");
        MPI_Finalize();
        return 0;
    }   

    if(myid == 0){
        time_start =  MPI_Wtime(); //Initialize a time, to measure the duration of the processing time.
    }     
    

    mpi_zeta();
    MPI_Finalize();
    return 0;
}

/* Return how many tasks a given process should do, an attept at load balancing */
std::tuple<int, int> calculate_work(){
    if(myid == -1){
        return std::make_tuple(0,n);
    }
      if(numprocs > n){//special case we need to handle, if there are more processes than tasks
        //e.g n=4 and np=8. the first 4 processes should get one task each, the rest 0.
        if(n > myid){ 
            return std::make_tuple(myid+1,1); //one task for first n processes.
        }else{
            return std::make_tuple(0,0);//no work to be done
        }
    }
    int division = n/numprocs;
    int remainder = n%numprocs;
    int start = myid*division + 1; //start position from the n tasks. Shift as needed to not overlap work area.
    int m = division; //minimum amount of work for each process.
    if(remainder != 0){
        if(remainder > myid){//We just give the first remainder processes one extra task.
            m += 1; //add one extra task.
            start += myid; //We need to shift start position by myid.
            return std::make_tuple(start, m);
        }else{
            start += remainder; //We need to shift start position by remainder.
            return std::make_tuple(start, m);
        }
    }else{//we can divide tasks cleanly
        return std::make_tuple(start, m);
    }
}

void mpi_zeta(){
    auto work = calculate_work();
    int start = std::get<0>(work);
    int n = std::get<1>(work);
    //int end_interval = n == 0 ? 0 : start+n-1; //correctly label the end_interval in printout.

    double sum = zeta(start, n, threads);
    //printf("Process %d calculate interval: [%d , %d]\n",myid, start, end_interval);
    double sum_all = 0.0;
    
    MPI_Allreduce(&sum, &sum_all, 1 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); //we sum the sum variable from each process and store in pi.
    
    printf("Process %d has sum %f\n",myid, sum);

    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        double pi = sqrt(sum_all*6); 
        double error = fabs(pi - 4.0 * atan(1.0));
        printf("pi = %.15g, error=%.15f, duration=%f ms\n",pi, error, duration*1000);
    }
}