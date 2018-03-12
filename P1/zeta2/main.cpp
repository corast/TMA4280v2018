#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include "zeta2.h"

void zeta();
std::tuple<int, int> calculate_work();

void recursive_doubling_mach(double *,double *);
void recursive_doubling_zeta(double *,double *);

//Global variables
int numprocs, n, method, type, threads;
int myid;
std::clock_t start;

int main(int argc, char* argv[])
{   
    //We need to pass more arguments.
    if(argc != 3){
        std::printf("Please type number of iteration, method of computation(0 = zeta, 1 = mach) and type of summation(0 = MPI_Allreduce, 1 = recursive-doubling) \n");
        return 1;
    }
    n = std::stoi(argv[1]);//Number of itterations.
    threads = std::stoi(argv[2]);//Number of threads

    /*method = std::stoi(argv[2]); //What method of computation we use, zeta or mach.
    type = std::stoi(argv[3]); //What type of summation we use, allreduce or recursive-doubling sum.
    
    MPI_Init(&argc, &argv); //init MPI

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//Get the number of processors.
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); //Get my rank(id)

    //Check if the number of processes are not a power of 2.(1,2,..,2^n)
    if( (numprocs & (numprocs - 1)) != 0 && numprocs) { 
        printf("Number of processes are not a number of 2\n");
        MPI_Finalize();
        return 0;
    } */

    
    start = std::clock();
    //MPI_Finalize();
    zeta();
    return 0;
}

/* Return how many tasks a given process should do, an attept at load balancing */
std::tuple<int, int> calculate_work(){
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

void zeta(){
    double sum = zeta(n, threads);
    double duration  = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
    double pi = sqrt(sum*6);
    double error = fabs(pi - 4.0 * atan(1.0));
    printf("pi = %e, error=%e, duration=%e ms\n",pi, error, duration*1000);
}