#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <bitset>
#include <tuple>
#include <mpi.h>

#include "mach2.h"

void mpi_mach();
std::tuple<int, int> calculate_work();

//Global variables
int numprocs, n, threads;
int myid = -1; 
double time_start;

int main(int argc, char* argv[])
{   
    //MPI_Status status;
    //We need to pass more arguments.
    if(argc != 3){
        std::printf("Please type number of iteration, method of computation(0 = zeta, 1 = mach) and type of summation(0 = MPI_Allreduce, 1 = recursive-doubling) \n");
        return 1;
    }
    n = std::stoi(argv[1]);//Number of itterations.
    threads = std::stoi(argv[2]);//Number of threads.

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
    }

    mpi_mach();
    MPI_Finalize();
    return 0;
}

/* Return how many tasks a given process should do, an attept at load balancing */
std::tuple<int, int> calculate_work(){
    if(myid == -1){
        return std::make_tuple(1,n);
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

void mpi_mach(){
    double arctans[2]; //each process should store their corresponding calculation in this array.

    auto work = calculate_work(); //calculate work-load for this process.

    int start = std::get<0>(work);
    int m = std::get<1>(work);
    //int end_interval = m == 0 ? 0 : start+m-1; //correctly label the end_interval in printout for processes with no work.
    //printf("Working interval [%d, %d] \n", start, end_interval);
    //Do the work
    arctans[0] = arctan(start, m, (double)1/5, threads);
    arctans[1] = arctan(start, m, (double)1/239, threads);
    //printf("Process %d calculate interval: [%d , %d]\n",myid, start, end_interval);

    double arctans_all[2]; //Hold the final sum
    //MPI_Reduce on the adresses of arctans and arctans_all
    MPI_Allreduce(arctans,arctans_all, 2 ,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    if(myid == 0 || myid == -1){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        //double duration  = ( std::clock() - start_time ) / (double)CLOCKS_PER_SEC;
        double pi = 4*(4*arctans[0] - arctans[1]);
        double error = std::abs(pi - (4.0 * atan(1.0)));
        printf("Threads=%d pi=%.15f, error=%.15f, duration=%f ms\n",threads ,pi, error, duration*1000);
        //std::cout <<"pi = "<< pi <<", error= "<< error <<", duration= " << duration*1000 <<" ms" <<std::endl;
    }
}

bool forwarSendDirection(int i){
    /* Return whether or not a given process should MPI_sendrecv forward the given distance or not   
    */
    //We need to create an mask corresponding to itteration, to check a given bit.
    int mask = 1;
    if(i > 0){
        mask = mask << i;//we need to shift our mask by i bits. mask*2^î 
    }
    int sendDir = mask ^ myid; //xor the mask with myid.
    //if the bit position corresponding to the mask is set, we need to send the sum a given distance.
    std::bitset<8> bit (sendDir); //8 bits should be enough.
    //if bit nr i has value 1 we should send forward.
    if(bit.test(i)){
        return true;     
    }//else we need to send backwards.
    return false;
}

void recursive_doubling_mach(double *arctans, double *arctans_all){
    //Recursive doubling: ech process pass on the sum and calculate a partial sum from this. First all processes send their sum once, 
    /*
        IF we have 8 processes. Everyone process hold a partial sum and pass on the value. 
        P0<->P1 P2<->P3 P4<->P5 P6<->P7
        P0<->P2 P1<->P3 P4<->P6 P5<->P7
        P0<->P4 P1<->P5 P2<->P6 P3<->P7

        #Ather this exchante og sums, every process should process should be able to compute pi.
        MPI_Sendrecv(sendbuffer, sendcount, Datatype, destination, sendtag, receivbuffer, recievcount, Datatype, source , recvtag, comm );
    */
    //MPI_Status status;
   
    double recbuff[2] = {0,0}; //Local reciev buffer
    double sum[2] = {arctans[0],arctans[1]}; //Local sum
    for(int i = 0; i < log(numprocs); i++){
        int distance = (int)pow(2,i); //This is the distance between who is sending and who is receiving.
        if(myid == 0){
            //printf("Distance %d at i = %d\n",distance,i);
        } 
        int sendId = myid;

        if(forwarSendDirection(i)){
            sendId += distance; //send distance forward 
        }else{
            sendId -= distance; //send diastance backward
        }

        //MPI_Send(&sum, 2, MPI_DOUBLE, sendId, 0, MPI_COMM_WORLD);
            //double recievBuffer = 0;
        //MPI_Recv(&recbuff, 2, MPI_DOUBLE, sendId, 0, MPI_COMM_WORLD, &status);

        for(int j = 0; j<2; j++){
            //update the values.
            sum[j] += recbuff[j];
        }
        
    }
    //Write to answer to the buffer, this is simpy so that we don't need to change anything else
    for(int i = 0; i<2; i++){
        //Update with final valyes
        arctans_all[i] = sum[i];
    }
}