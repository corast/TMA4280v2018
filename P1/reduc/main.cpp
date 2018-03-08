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
void recursive_doubling(void *,void *);

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
    int start = std::get<0>(work);
    int n = std::get<1>(work);
    int end_interval = n == 0 ? 0 : start+n-1; //correctly label the end_interval in printout.

    double sum = zeta(start, n);
    //printf("Process %d calculate interval: [%d , %d]\n",myid, start, end_interval);
    double sum_all = 0.0;
    
    if(!type){//If we want all_reduce
        MPI_Allreduce(&sum, &sum_all, 1 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); //we sum the sum variable from each process and store in pi.
    }else{
        recursive_doubling(&sum,&sum_all); //Where we want the final sum to be stored.
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


    //Do the work
    arctans[0] = arctan(start, m, (double)1/5);
    arctans[1] = arctan(start, m, (double)1/239);
    printf("Process %d calculate interval: [%d , %d]\n",myid, start, end_interval);

    double arctans_all[2]; //Hold the final sum
    //MPI_Reduce on the adresses of arctans and arctans_all
    if(!type){
        MPI_Allreduce(arctans,arctans_all, 2 ,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }else{
        recursive_doubling(arctans,arctans_all); //Where we want the final sum to be stored.
    }
    

    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        double pi = 4*(4*arctans_all[0] - arctans_all[1]);
        double error = std::abs(pi - (4.0 * atan(1.0)));
        printf("pi=%.15g, error=%e, duration=%f ms\n",pi, error, duration*1000);
        //std::cout <<"pi = "<< pi <<", error= "<< error <<", duration= " << duration*1000 <<" ms" <<std::endl;
    }
}

bool turnToSend(int distance, int x){
    //to generate the sequece of processes to send this turn
    // 2x*(distance to reciever process)
    //We need to calculate the values, and check if its equal to our process, if its higher we return false.
    /*
       x = 0      1, 3, 5, 7, 9, 11, ... (2*i + 1)
       x = 1      2, 6, 10, 14, 18, ... (2*2*i + 2)
       x = 2      4, 12, 20, ... (2*2*2*i + 4)
    */ 
    int proc = 0;
    for(int i = 0; i<numprocs; i++){ 
        proc = (pow(2,x+1)*i)+distance; //Sequence is (2^x)*i+distance
        if(myid == 4){
            printf("itteration=%d i=%d proc=%d distance=%d :  %d\n",x,i,proc,distance, 2*i+distance);
        }
        if(proc == myid ){//We are one of the senders this turn
            return true;
        }else if(proc > myid){//we are not one of the senders.
            return false;
        }
    }  
}

void recursive_doubling(void *psum,void *recvbuf){
    MPI_Status status;

    //Recursive doubling: ech process pass on the sum and calculate a partial sum from this. First all processes send their sum once, 
    /*
        IF we have 8 processes. Everyone process hold a partial sum and pass on the value. 
        P0<-P1 P2<-P3 P4<-P5 P6<-P7
        P0<----P2     P4<----P7
        P0<-----------P4
        We can now calculate pi on P0 with the final sum.
    */
    for(int i = 0; i < log(numprocs); i++){
        //Every odd numbered process need to send to its -1 process

        //Very uneffective test
        /*int ps = numprocs/pow(2,i+1);
        int array[ps]; */

        int distance = (int)pow(2,i); //This is the distance between who is sending and who is receiving.
        
        if(turnToSend(distance,i)){
            printf("i %d, process %d to send to %d, distance %d \n", i, myid, myid-distance, distance);
            if(method){//if we come from zeta
                MPI_Send(&psum,1,MPI_DOUBLE, myid-distance,0,MPI_COMM_WORLD);
            }else{
                //MPI_Sendrecv(&psum,1, )
            }
            
        }

        /*
        //now we need to figure out if we are a sender or a reciever.
        if(i+distance == myid){
            printf("i %d, sender:%d\n",i , (int)pow(2,i));
        }
        if( (myid-distance)%2 == 0 && i == 0){ //This apply only for the first itteration.
            //First itteration this work.
        }else{
            //distance is a number of 2 from here on out.
            if(myid%distance == 0){
                //then we need to send.
            }

        }*/
    }
    //We need to write the sum
    //How does this work?
}