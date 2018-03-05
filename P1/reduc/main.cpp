#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <mpi.h>
#include "mach1.h"
#include "zeta1.h"

void mpi_zeta();
void mpi_mach();
std::tuple<int, int> calculate_work();
void recursive_doubling(void *recvbuf);

//Global variables
int numprocs, n, method, type;
int myid;
double time_start;

int main(int argc, char* argv[])
{   
    MPI_Status status;
    //We need to pass more arguments.
    if(argc != 4){
        std::printf("Please type number of iteration, method of computation(0 = zeta, 1 = mach) and type of summation(0 = MPI_Allreduce, 1 = recursive-doubling) \n");
        return 1;
    }
    n = std::stoi(argv[1]);//Number of itterations.
    method = std::stoi(argv[2]); //What method of computation we use, zeta or mach.
    type = std::stoi(argv[3]); //What type of summation we use, allreduce or recursive-doubling sum.
    
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
        std::string m[2] = {"zeta","mach"};
        std::string t[2] = {"Allreduce","Recursive-doubling"};
        std::cout << "Using method "<< m[method] << " with summation type " << t[type] << std::endl;
        time_start =  MPI_Wtime(); //Initialize a time, to measure the duration of the processing time.
    }

    switch(method){
        case 0:{
            mpi_zeta();
        }break;

        case 1:{
            mpi_mach();
        }break;

        default: {
            
        }
    }

    MPI_Finalize();
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

void mpi_zeta(){
    auto work = calculate_work();

    double sum = zeta(std::get<0>(work),std::get<1>(work), myid);
    double sum_all = 0.0;
    
    if(!type){//If we want all_reduce
        MPI_Allreduce(&sum, &sum_all, 1 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); //we sum the sum variable from each process and store in pi.
    }else{
        recursive_doubling(&sum_all); //Where we want the final sum to be stored.
    }

    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        double pi = sqrt(sum_all*6); 
        double error = fabs(pi - 4.0 * atan(1.0));
        printf("pi = %e, error=%e, duration=%e ms\n",pi, error, duration*1000);
    }
}

void mpi_mach(){
    double arctans[2]; //each process should store their corresponding calculation in this array.

    auto work = calculate_work(); //calculate work-load for this process.

    int start = std::get<0>(work);
    int m = std::get<1>(work);
    int end_interval = m == 0 ? 0 : start+m-1; //correctly label the end_interval in printout for processes with no work.
    printf("Process %d calculate interval: [%d , %d]\n",myid, start, end_interval);

    //Do the work
    arctans[0] = arctan(start,m,(double)1/5, myid);
    arctans[1] = arctan(start,m,(double)1/239, myid);

    //Do the work
    arctans[0] = arctan(std::get<0>(work),std::get<1>(work),(double)1/5, myid);
    arctans[1] = arctan(std::get<0>(work),std::get<1>(work),(double)1/239, myid);

    double arctans_all[2]; //Hold the final sum
    //MPI_Reduce on the adresses of arctans and arctans_all
    if(!type){
        MPI_Allreduce(arctans,arctans_all, 2 ,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }else{
        recursive_doubling(arctans_all); //Where we want the final sum to be stored.
    }
    

    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        double pi = 4*(4*arctans_all[0] - arctans_all[1]);
        double error = std::abs(pi - (4.0 * atan(1.0)));
        printf("pi=%.15g, error=%e, duration=%f ms\n",pi, error, duration*1000);
        //std::cout <<"pi = "<< pi <<", error= "<< error <<", duration= " << duration*1000 <<" ms" <<std::endl;
    }
}

void recursive_doubling(void *recvbuf){
    //How does this work?
}