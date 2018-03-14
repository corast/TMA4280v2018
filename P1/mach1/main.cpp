#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <tuple> //For returning multple variables from an function, could use an array too.
#include <mpi.h>

#include "mach1.h"
 
/* Function decleration */
double mpi_mach();
std::tuple<int, int> calculate_work();

/* Global variables */
int numprocs, n;
int myid;
double time_start;
double *pnt_array;

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
    if( (numprocs & (numprocs - 1)) != 0 && numprocs) { 
        printf("Number of processes are not a number of 2\n");
        MPI_Finalize();
        return 0;
    }
    
    if(myid == 0){ //We are master.
        time_start =  MPI_Wtime(); //Initialize a time, to measure the duration of the processing time.
    }else{ //We are a slave
        ; //nothing to do
    }

    mpi_mach();

    MPI_Finalize();
    return 0;
}
/* Return how many tasks a given process should do */
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
}//std::get<0>(tuple_name) for element 0.

double mpi_mach(){
    //MPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuff, recvcount, recvtype, 0, MPI_COOM_WORLD);

    double arctans[2]; //each process should store their corresponding calculation in this array.

    auto work = calculate_work(); //calculate work-load for this process.
    int start = std::get<0>(work);
    int m = std::get<1>(work);
    int end_interval = m == 0 ? 0 : start+m-1; //correctly label the end_interval in printout for processes with no work.
    //printf("Process %d calculate interval: [%d , %d]\n",myid, start, end_interval);

    //Do the work
    arctans[0] = arctan(start,m,(double)1/5);
    arctans[1] = arctan(start,m,(double)1/239);
    //To more cleanly represent the results, remove this output, to see the result on each process, uncomment.
    //printf("Process %d interval [%d,%d] calculated arctan_1/5 %.16f, arctan_1/239 %.16f\n", myid, start, end_interval, arctans[0], arctans[1]);
    double arctans_all[2]; //Hold the final sum
    //MPI_Reduce on the adresses of arctans and arctans_all
    MPI_Reduce(arctans,arctans_all, 2 ,MPI_DOUBLE, MPI_SUM, 0 , MPI_COMM_WORLD);

    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        double pi = 4*(4*arctans_all[0] - arctans_all[1]);
        double error = std::abs(pi - (4.0 * atan(1.0)));
        printf("cpi =%.16g pi = %.16g, error=%.16f, duration=%f ms\n", 4.0 * atan(1.0),pi, error, duration*1000);
        //std::cout <<"pi = "<< pi <<", error= "<< error <<", duration= " << duration*1000 <<" ms" <<std::endl;
    }
}