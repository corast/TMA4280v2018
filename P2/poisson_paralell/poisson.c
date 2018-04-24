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

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);

void transpose_paralell(real **bt, real **b, size_t m);
void calculate_sendcounts(size_t m);
void calculate_workingarea(size_t m);

void fill_rec_recdiv(size_t m, int *recvc, int *recvdis, int *sendc, int *sdispls);

void divide_work(size_t m);

void create_mpi_datatype(size_t m);
void free_mpi_datatype();
//tuple calculate_work(size_t rows);
// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);


//debug functions.
void printMatrix(real** matrix, int size, char* c);
void printVector(real* vector, int size, char* c);
void printArray(int *array, int size, char* c);
void printDebug(real** matrix, int size);
//end debug functions.

bool dbug = true;



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
        if(myid == 0){
            printf("the problem size n must be a power of 2\n");
        }
    
        MPI_Finalize();
        return 0;
    }

    if(myid == 0){
        printf("Running with %d processes and %d threads on each process\n",numprocs,threads);
        time_start =  MPI_Wtime(); //Initialize a time, to measure the duration of the processing time.
    } 


    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }
    //printVector(grid, n, "grid");

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
    real *z = mk_1D_array(nn, false);
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
            b[i][j] = h * h * rhs(grid[i+1], grid[j+1]);
        }
    }

    create_mpi_datatype(m);//create datatype 

    divide_work(m); //divide the work, filling the nececery arrays
    

    //calculate_sendcounts(m);
    //calculate_workingarea(m);

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
        fst_(b[i], &n, z, &nn);
    }
    //printMatrix(b,m,"b2 line 159");
    transpose_paralell(bt,b,m);
    //transpose(bt, b, m);
    //printMatrix(bt,m,"bt1 line 162");
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        fstinv_(bt[i], &n, z, &nn);
    }
    //printMatrix(bt,m,"bt2 line 166");
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
    //printMatrix(bt,m,"bt3 line 177");
    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        fst_(bt[i], &n, z, &nn);
    }
    //printMatrix(bt,m,"bt4 line 186");
    transpose_paralell(b, bt, m);
    #pragma omp parallel for num_threads(threads)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        fstinv_(b[i], &n, z, &nn);
    }
 
    /*
     * Compute maximal value of solution for convergence analysis in L_\infty
     * norm.
     */
    double u_max = 0.0;
    #pragma omp parallel for num_threads(threads) collapse(2)//Parellalize the code with threads. 
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            u_max = u_max > b[i][j] ? u_max : b[i][j];
        }
    }

    if(myid == 0){//process zero should do the final calculation.
        double duration  = MPI_Wtime() - time_start;
        printf("duration = %f ms\n", duration*1000);
        printf("u_max = %e\n", u_max);
    }

    //free some memory.
    free_mpi_datatype();

    MPI_Finalize();
    return 0;
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs(real x, real y) {
    return 2 * (y - y*y + x - x*x);
}

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

void transpose(real **bt, real **b, size_t m)
{
    //printMatrix(b,m,"transpose b");

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
        }
    }
    //printMatrix(bt,m,"transpose bt");
}

int *sendcounts;
int *recvcounts;
int *sdisplacements;
int *recvdisplacements;

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
    if(myid == 0){
        //printArray(recvcounts_2,numprocs, "revcounts_2");
        //printArray(recvdisplacements_2,numprocs, "revcountsdisplacements_2");  
    }
    //printArray(sendcounts_2, numprocs, "sendcounts_2");
    //printArray(sdisplacements_2, numprocs, "sdisplacements_2");
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
        }
        //Else this process have no rows to send.(handled in calloc)
        //else should have all zeros as recvc[p] and recvdis[p]
    }
    //cant update send and sdispls in the same loop, unless we are process 0.
    for(int p = 0; p<numprocs; p++){
        sendc[p] = recvc[myid]*m; //number of rows we are reponsible of * elements per row.
        //TODO: need to hande when the number of elements to send ist the same as previous rows.
        sdispls[p] = recvdis[myid]*m;
    }
}

MPI_Datatype type_matrix;
MPI_Datatype column;

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


void transpose_paralell(real **bt, real **b, size_t m){
    //m is the amount of data points this process should send.
    //printMatrix(b,m,"b");
    //Note, sending doubles, resulting in number of elements to send in sendcount ect, but receiving in matrix columns(rows)
    MPI_Alltoallv(b[0],sendcounts, sdisplacements, MPI_DOUBLE, bt[0], recvcounts, recvdisplacements, type_matrix, MPI_COMM_WORLD);

    //printMatrix(bt,m,"transpose b");
}

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


//Printout function, to debug and figure out what the code actualy does
void printMatrix(real** matrix, int size, char* c){
    if(myid != 0){
        return;
    }
    printf("%s ->\n",c);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
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

void printDebug(real** matrix, int size){
     if(myid != 0){
        return;
    }
    printf("\n");
    for(int i = 0; i<size; i++){
        printf("%p ",matrix[i]);
        //printf("%f ", **(matrix++));
        printf("\n");
    }
    printf("\n%p",*matrix);
    printf("\n%p",matrix++);
    /*
    printf(" 1 %p ",matrix[i]);
        printf(" 2 %p ",matrix[0+1]);
        printf(" 3 %p ", *matrix);
        printf(" 4 %p ", *(matrix++));
        printf("\n");
    */
}