#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <mpi.h>

#include "mach1.h"

double mpi_mach();

//Global variables
int numprocs, n;
int myid;
double time_start;
double *pnt_array;

int *global_i;
int li[2];

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
        pnt_array = (double*) calloc(n,sizeof(double));
        //initialize this array with the elements that we want to share between the processes.
       /* for(int i=0; i<n; i++){
            //++pnt_array;
            printf("setting %f to %d\n",*pnt_array, i+1);
            *(pnt_array+i) = (double)i+1; //set the value of this index as 
        }
        for(int i = 0; i<n; i++){
            printf("%f ",*(pnt_array+i));
        }printf("\n");
        */
    }else{ //We are a slave
        ; //nothing to do
    }

    mpi_mach();

    MPI_Finalize();
    return 0;
}

double mpi_mach(){

    MPI_Scatter(global_i, 2, MPI_INT, li, 2 ,MPI_INT, 0, MPI_COMM_WORLD);
    /*
    double arctans[2]; //an array that holds arctan1 and arctan2.
    
    int divide = (double)n/numprocs; //compute the minimum amount of work for each process.
    int rest_tasks = n%numprocs;
    if(rest_tasks != 0){ //If we can't cleanly divide the number of tasks between the processes
        //We need to check how many extra tasks there are.
        //We simply need to give rest_task processes an extra amout of work to complete
        int tasks = n; //We divide this out to the tasks.

        if(rest_tasks > myid){ //we want the tasks with id higher than the work to do the extra task.
            int start = myid*divide + myid + 1;
            int n = divide + 1;
            arctans[0] = arctan(start, n,1/8 , myid); //each process calculate it's sum.
            arctans[1] = arctan(start, n,1/256, myid);
        }else
        {   //We need to shift the work intervals by the rest_tasks.
            int start = myid*divide + 1 + rest_tasks;
            int n = divide;
            arctans[0] = arctan(start, n,1/8, myid); //each process calculate it's sum.
            arctans[1] = arctan(start, n,1/256, myid);
        }
    }else{//it divides cleanly between processes.
        
        arctans[0] = arctan(myid*divide + 1, divide,1/8, myid); //each process calculate it's sum.
        arctans[1] = arctan(myid*divide + 1, divide,1/256, myid); //each process calculate it's sum.
    }
    */

 

    //MPI_Reduce(&sum, &sum_all, 2 , MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); //we sum the sum variable from each process and store in pi.

    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        //double arctan_1 = arctans[0];
        //double arctan_2 = arctans[1];
        
        //double pi = 4*(4*arctan_1 - arctan_2);
        double pi = 1;
        double error = fabs(pi - 4.0 * atan(1.0));
        free(pnt_array);
        printf("pi = %e, error=%e, duration=%e ms\n",pi, error, duration*1000);
    }
}