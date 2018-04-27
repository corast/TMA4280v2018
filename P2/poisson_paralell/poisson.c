/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

int myid = -1;
int numprocs, n, threads;
int time_start;
double time_divide_work;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y, int mode);

real u(real x, real y);
real f(real x, real y);

void transpose_paralell(real **bt, real **b, size_t m);
void calculate_sendcounts(size_t m);
void calculate_workingarea(size_t m);

void fill_rec_recdiv(size_t m, int *recvc, int *recvdis, int *sendc, int *sdispls);

void divide_work(size_t m);

void create_mpi_datatype(size_t m);
void free_mpi_datatype();
/*####Validation functions etc*/
double findGlobalUmax(real **b, size_t m);
double calcualteGlobalError(real **b, size_t m, real *grid);

void fillSolutionMatrix(real **U, real *grid ,size_t m);
//tuple calculate_work(size_t rows);
// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);


//debug functions.
void testTranspose(real** start, real**end, size_t m);

void printMatrix(real** matrix, size_t n, size_t m, char* c);
void printVector(real* vector, int size, char* c);
void printArray(int *array, int size, char* c);
real **createTransposeMatrix(size_t m);
//end debug functions.

int *sendcounts;
int *recvcounts;
int *sdisplacements;
int *recvdisplacements;

size_t start; //from what row we are responsible for handling the data.
size_t end; //to what row we stop.

int mode = 2;

//used as sendtype in MPI_Alltoallv, gives us the ability to send whole matrix.
MPI_Datatype type_matrix;
MPI_Datatype column;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv); //init MPI

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//Get the number of processors.
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); //Get my rank(id)
    
    if (argc < 3) {
        if(myid == 0){
            printf("Usage:\n");
            printf("  poisson n t\n\n");
            printf("Arguments:\n");
            printf(" n: the problem size (must be a power of 2)\n");
            printf(" t: the number of threads \n");
        }
        MPI_Finalize();
        return 0;
    }

    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */
    
    int n = atoi(argv[1]); //Problem size 

    //Assume the number of threads are an integer, not some char
    int threads = atoi(argv[2]); //Get number of threads per process.
    int m = n - 1; //Set degrees of freedom.
    real h = 1.0 / n; //Set mesh size.

 
    if( (n & (n - 1)) != 0 && n) { //Check that the problem size is a power of 2.
        if(myid == 0){//remove clutter, don't want np outputs.
            printf("the problem size n must be a power of 2\n");
        }
    
        MPI_Finalize();
        return 0;
    }

    if(myid == 0){
        //printf("Running with %d processes and %d threads on each process\n",numprocs,threads);
        time_start =  MPI_Wtime(); //Initialize a time, to measure the duration of the processing time.
    } 

    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads.  Note: default shedule is static(each is given a fixed amount)
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }


    /*
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    real *diag = mk_1D_array(m, false);
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }
    //printVector(diag, m, "d1 line 103");

    /*
     * Allocate the matrices b and bt which will be used for storing value of
     * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
     */
    real **b = mk_2D_array(m, m, false);
    real **bt = mk_2D_array(m, m, false);

    /*
     * This vector will holds coefficients of the Discrete Sine Transform (DST)
     * but also of the Fast Fourier Transform used in the FORTRAN code.
     * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
     * - Fourier coefficients are complex so storage is used for the real part
     *   and the imaginary part.
     * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while 
     *   DST coefficients are defined for j [[ 0, n-1 ]].
     * As explained in the Lecture notes coefficients for positive j are stored
     * first.
     * The array is allocated once and passed as arguments to avoid doings 
     * reallocations at each function call.
     */
    int nn = 4 * n;
    /*
        We must give each tread their own z vector to work with. 
    */
    real *z[threads];
    for(int i = 0; i < threads; i++){ 
        z[i] = mk_1D_array(nn, false);
    }
    //real *z = mk_1D_array(nn, false);
    //printVector(z, nn,"z1 line 129");

    /*
     * Initialize the right hand side data for a given rhs function.
     * Note that the right hand-side is set at nodes corresponding to degrees
     * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
     * 
     */
    #pragma omp parallel for num_threads(threads) collapse(2)//Parellalize the code with threads, 
    //Paralellalize to use threads on nested loop as well.
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            b[i][j] = h * h * rhs(grid[i+1], grid[j+1], mode);
        }
    }

    create_mpi_datatype(m);//create datatype 

    time_divide_work =  MPI_Wtime(); //Initialize a time, to measure the duration of the processing time.
    divide_work(m); //divide the work, filling the nececery arrays
    
    
    start = recvdisplacements[myid]; //from what row we are responsible for in the solution matrix b, calculating error etc.
    end = recvdisplacements[myid]+recvcounts[myid]; //to what row we need stop.

    time_divide_work = MPI_Wtime() - time_divide_work; //end time. 
    /*
     * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
     * Instead of using two matrix-matrix products the Discrete Sine Transform
     * (DST) is used.
     * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
     * The array zz is used as storage for DST coefficients and internally for 
     * FFT coefficients in fst_ and fstinv_.
     * In functions fst_ and fst_inv_ coefficients are written back to the input 
     * array (first argument) so that the initial values are overwritten.
     */
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        fst_(b[i], &n, z[omp_get_thread_num()], &nn);
    }

    transpose_paralell(b,bt,m); //transpose b matrix and put into bt.

    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        fstinv_(bt[i], &n, z[omp_get_thread_num()], &nn);
    }
    //done in O(n² log n)

    /*
     * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
     */
    #pragma omp parallel for num_threads(threads) collapse(2)//Parellalize the code with threads.
    //collapse to run the nested loop with threads aswell. 
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[i] + diag[j]);
        }
    }
    //done in O(n²)

    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        fst_(bt[i], &n, z[omp_get_thread_num()], &nn);
    }
    transpose_paralell(bt, b, m); //transpose bt and put into b.
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        fstinv_(b[i], &n, z[omp_get_thread_num()], &nn);
    }
 
    /*
     * Compute maximal value of solution for convergence analysis 
     */
    double global_umax = findGlobalUmax(b,m);

    //calculate error
    double global_error = calcualteGlobalError(b,m,grid);

    #if 0 //set to one for testing transposing one matrix, which is easy to vertify by eye. 
    real** holderMatrix = mk_2D_array(m,m,true);
    real** transposeMatrix = createTransposeMatrix(m);

    testTranspose(transposeMatrix, holderMatrix, m);
    #endif

    #if 0 //adjust to 1 to print b and solution matrix.
    printMatrix(b,m,m,"U");
    /*
     * Solution matrix, to compare the values.
     */
    real **solU = mk_2D_array(m,m, false); 
    fillSolutionMatrix(solU,grid , m);
    printMatrix(solU,m,m, "solutions U");
    #endif
    
    #if 1 //question3-2
    if(myid == 0){
        printf("np =%3d,myid = %d, m =%6d, work_divide_duration %8.3f ms \n",numprocs, myid , m, time_divide_work*1000);
    }
    #endif

    #if 0 //question 4
    //we want to createa a csv file as output, which we can plot in python.
    //we want to plot time against number of processes, treads is set to 1.
    if(myid == 0){
        double duration  = MPI_Wtime() - time_start;
        printf("thr_p:%3d, np =%3d, n =%6d, duration = %8.2f ms, u_max = %8f, error_max = %15.15f \n", threads, numprocs, n, duration*1000, global_umax, global_error);
    }
    #endif

    #if 0
    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        printf("thr_p:%3d, np =%3d, n =%6d, duration = %8.2f ms, u_max = %8f, error_max = %15.15f \n", threads, numprocs, n, duration*1000, global_umax, global_error);
    }
    #endif
    
    //free some memory.
    free_mpi_datatype();

    MPI_Finalize();
    return 0;
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs(real x, real y, int mode) {
    switch(mode){
        case 0:{
            return 2 * (y - y*y + x - x*x);
        }break;
        case 1:{
            return 1;
        }break;
        case 2:{ //to test our solution with an exact amount.
            return 5*PI*PI*sin(PI*x)*sin(2*PI*y);
        }
        default:
            return 2 * (y - y*y + x - x*x);
    }
}

real solution(real x, real y, int mode){
    switch(mode){
        case 0:{
            return 1;
        }
        case 1:{
            return 1;
        }break;
        case 2:{ //exact solution to the problem.
            return sin(PI*x)*sin(2*PI*y);
        }
    }
}


void transpose_paralell(real **b, real **bt, size_t m){
    //m is the amount of data points this process should send.
    //Note, sending doubles, resulting in number of elements to send in sendcount ect, but receiving in matrix columns(rows)
    MPI_Alltoallv(b[0],sendcounts, sdisplacements, MPI_DOUBLE, bt[0], recvcounts, recvdisplacements, type_matrix, MPI_COMM_WORLD);
}

/*##############################
    Dividing matrix functions
################################*/
void divide_work(size_t m){
    //we know that every process has to send an equal amount of data(elemnts to every other process)

    //first we calculate how many rows each process should be responsible for.
    /*
        If m = nprocs, then nrows = 1.
        If m < nprocs, then first m procs get to send one.(rest shouldnt do anything other than recieve)
        If m > nprocs, then we divide up the rows to m/nprocs each, and give the first m%nprocs one extra row to send.
    */

    sendcounts = (int*) malloc(numprocs*sizeof(int));
    recvcounts = (int*) calloc(numprocs,sizeof(int));
    sdisplacements = (int*) malloc(numprocs*sizeof(int));
    recvdisplacements = (int*) calloc(numprocs,sizeof(int));

    fill_rec_recdiv(m, recvcounts, recvdisplacements, sendcounts, sdisplacements);

    #if 0 //write the different buffers needed for MPI_alltoallv, for debugging.
    if(myid == 0){
        printArray(recvcounts,numprocs, "revcounts");
        printArray(recvdisplacements,numprocs, "revcountsdisplacements");  
    }
    printArray(sendcounts, numprocs, "sendcounts");
    printArray(sdisplacements, numprocs, "sdisplacements");
    #endif
}

void fill_rec_recdiv(size_t m, int *recvc, int *recvdis, int *sendc, int *sdispls){
    //basicly use the same function as last project, to load balance the work, but we need to know everyones work.
    int division = m/numprocs;
    int remainder = m%numprocs;
    for(int p = 0; p<numprocs; p++){//fill the array for each process.
        if(numprocs >= m){ //we need to divide to the m first prosesse
            if(p < m){ //give the first m processes a task
                recvc[p] = 1;
                recvdis[p] = p;
            }
        }else if(remainder != 0){
            if(remainder > p){//first p processes are given one extra row to work with.
                recvc[p] = division + 1; //we given extra rows
                recvdis[p] = p*division + (p*1) ; //we need to offset with this extra row.
            }else{//the processes that do not recieve extra work, should offset by one.
                recvc[p] = division;
                recvdis[p] = p*division + remainder; //we know that the previous processes recieved one extra task, 
            }
        }else{ //we know that each prosess has an equal amount of rows.
            recvc[p] = division;
            recvdis[p] = p*division;
        }
        //Else this process have no rows to send.(handled in calloc)
        //else should have all zeros as recvc[p] and recvdis[p]
    }
    //cant update send and sdispls in the same loop, unless we are process 0. Cant be paralalized.
    for(int p = 0; p<numprocs; p++){
        sendc[p] = recvc[myid]*m; //number of rows we are reponsible of * elements per row.
        sdispls[p] = recvdis[myid]*m; //number of rows we need to diplace * elements per row.
    }
}

void create_mpi_datatype(size_t m){//creat the custom datatypes for storing matrix as vectors, columwise.
    //Pack the data into custom datatypes for MPI, one column at a time. 
    //We want every process to be responsible for one row at a time. 
    // count = m, block_length = 1 (block per element), stride = m (how far to same element column in next row)
    /*  count = 3, block_length = 1, stride = 3. 
        a b c
        a b c
        a b c
    */
    MPI_Type_vector(m, 1, m, MPI_DOUBLE , &column);
    MPI_Type_commit(&column);//commit the datatype
    //lb = 0, extend = sizeof(double)
    MPI_Type_create_resized(column, 0, sizeof(double),&type_matrix); //duplicates the matrix datatype and changes the upper bound, lower bound and extent.
    MPI_Type_commit(&type_matrix); //commit the datatype
}


void free_mpi_datatype(){ // Free the created types after use from memory from each process.
    MPI_Type_free(&column);
    MPI_Type_free(&type_matrix);
}


/*##########################
     MK array functions
############################*/
/*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    real **ret = (real **)malloc(n1 * sizeof(real *));

    // 2
    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }
    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

/*##########################
        Print functions
############################*/

//Printout function, to debug and figure out what the code actualy does
void printMatrix(real** matrix, size_t n, size_t m, char* c){
    if(myid != 0){
        return;
    }
    printf("%s ->\n",c);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            printf("%f ", matrix[i][j]);
            //u_max = u_max > b[i][j] ? u_max : b[i][j];
        }
    printf("\n");
    }
}

//Printout function, to debug and figure out what the code actualy does
void printVector(real* vector, int size, char* c){ 
    if(myid != 0){
        return;
    }
    printf("%s ->\n",c);
    for (size_t i = 0; i < size; i++) {
        printf("%f ", vector[i]);
    }
    printf("\n");
}

void printArray(int *array, int size, char* c){
    printf("%s ->\n",c);
    for (size_t i = 0; i < size; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");
}
/*###########################
    Debugging functions.
#############################*/

void testTranspose(real** start, real** end, size_t m){
    /* Test the transpose function, by transposing twice, if the resulting matrix is the same as the origin, we know it works. */
    printMatrix(start, m, m,"before transpose");
    transpose_paralell(start, end, m);
    printMatrix(end, m,m, "after one transpose");
    transpose_paralell(end, start, m);
    printMatrix(start, m,m, "after two transposes");
    //we want to check whether or not start is equal to start.
}

/* create a matrix where every column is the same and increasing for every column */
real **createTransposeMatrix(size_t m){
    //We create an unit matrix and transpose it.
    //Simple test to check if our transpose is working correctly.

    real **transposeMatrix = mk_2D_array(m, m, true);//instanceiate.

    for(size_t i = 0; i < m; i++){
        for(size_t j = 0; j < m; j++){
            transposeMatrix[i][j] = i;
        }
    }

    return transposeMatrix;
}

/*####################################
        Erroc checking functions.
######################################*/

void fillSolutionMatrix(real **U, real *grid ,size_t m){
    real x,y;
    for(int i = 0; i<m; i ++){
        for(int j = 0; j<m; j++){
            x = grid[i+1];
            y = grid[j+1];
            U[i][j] = solution(x, y, mode);

        }
    }
}


double findGlobalUmax(real **b, size_t m){
    //each process find their u_max value and send to process 0. log2(p) time with one thread.
    double u_max = 0.0;
    double global_umax = 0.0;
    #pragma omp parallel for num_threads(threads) collapse(2)//Parellalize the code with threads. 
    for (size_t i = start; i < end; i++) { //each process only calculate their u_max from the assigned rows from the matrix.
        for (size_t j = 0; j < m; j++) {
            u_max = u_max > b[i][j] ? u_max : b[i][j];
        }
    }
    //need to MPI_reduce, with the max from every process of u_max.
    MPI_Reduce(&u_max, &global_umax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return global_umax;

}

double calcualteGlobalError(real **b, size_t m, real *grid){
/*  As x and y is a value between 0 and 1, we have to use the grid instead of index of array.
    Using U from assignment paper as reference to grid points.
        U       Grid (x)      Grid(y)   Grid = [0 1/4 2/4 3/4 1] Column_major.
    [a b c] [1/4 1/4 1/4] [1/4 2/4 3/4]
    [d e f] [2/4 2/4 2/4] [1/4 2/4 3/4]
    [g h i] [3/4 3/4 3/4] [1/4 2/4 3/4]
*/  
    double error = 0;
    double global_error = 0.0;
    real x;
    real y;
    double local_error = 0;
    //dont bother paralizing this one with threads.
    //loop trou every row and every column.
    for(size_t i = start; i < end; i++){ //each process only calculate their max error from the assigned rows from the matrix.
        for(size_t j = 0; j < m; j++){
            x = grid[i+1];//see example grids
            y = grid[j+1];
            local_error = fabs(solution(x, y, mode) - b[i][j]);
            error = error > local_error ?  error :  local_error; //update the error with the highest number.
        }
    }

    //need to MPI_reduce, with the max from every process of u_max.
    MPI_Reduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return global_error;
}