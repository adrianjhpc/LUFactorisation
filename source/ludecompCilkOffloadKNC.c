/*
 * Mateusz Iwo Dubaniowski, EPCC, 2015
 * Adrian Jackson, EPCC, 2016
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define TINY 1e-20 
#define _dim int n

//Declaration of matrix type
typedef double ** __restrict__  __attribute__((align_value (64))) mat;
 
/*
 * mat_zero function
 * This function populates a matrix with zeroes.
 */
__attribute__ ((target(mic))) //offloadable
void mat_zero(mat x, int n) {
	int i, j;
__assume_aligned(x, 64);
#pragma omp parallel for private(i,j) 
	for (i=0; i < n; i++)
#pragma vector
#pragma ivdep
			for (j=0; j < n; j++)
				x[i][j] = 0;
}
 
/*
 * mat_new function
 * This function allocates a new matrix.
 */
__attribute__ ((target(mic))) //offloadable
mat mat_new(_dim)
{
	int i;
	mat x = (double**)_mm_malloc(sizeof(double*) * n,64);		//Allocate array of pointers
	x[0]  = (double *)_mm_malloc(sizeof(double) * n * n,64);	//Allocate data array to the first pointer of the pointer array
        double *start;
        start = x[0];
	//Assign the array of pointers
#pragma omp parallel for private(i)
	for (i=1; i < n; i++)
		x[i] = start + n * i;
	mat_zero(x, n);		//Populate with zeroes
 
	return x;
}
 
/*
 * mat_copy function
 * This function copies matrix and returns a copy.
 */
mat mat_copy(void *s, _dim)
{
	int i, j;
	mat x = mat_new(n);
	double *ss=(double *) s;
	double *stemp=NULL;
	memcpy(x[0], s, n*n*sizeof(double));
	return x;
}
 
/*
 * mat_del function
 * Delete an array and free memory.
 */
__attribute__ ((target(mic))) //offloadable
void mat_del(mat x) { 
	_mm_free(x[0]);
	_mm_free(x);
}

 /*
 * mat_show function
 * Print a matrix to standard output.
 */
#define _QUOT(x) #x
#define QUOTE(x) _QUOT(x)
#define _show(a) printf(QUOTE(a)" =");mat_show(a, 0, n)
void mat_show(mat x, char *fmt, _dim)
{
	int i, j;

	if (!fmt) fmt = "%3.4g";
	for (i=0; i < n; i++){
		printf(i ? "      " : " [ ");
		for (j=0; j < n; j++){
			printf(fmt, x[i][j]);
			printf(j < n - 1 ? "  " : i == n - 1 ? " ]\n" : "\n");
		}
	}
}

 /*
 * mat_pivot function
 * Matrix pivoting using simple matrix transformations.
 */
__attribute__ ((target(mic))) //offloadable
void mat_pivot(mat a, int min, _dim){
__assume_aligned(a, 64);
  int i, j;
  double *temp;
  temp=_mm_malloc(n*sizeof(double),64);
  
  i = min;	
  int max_j = i;
  for (j=i; j < n; j++)
    if (fabs(a[j][i]) > fabs(a[max_j][i])) max_j = j;
  
  // Performing pivoting
  if (max_j != i){
    // Intel Cilk array notation introduced
__assume_aligned(temp, 64); 
__assume_aligned(a, 64);
__assume(i%16==0);
__assume(n%16==0);
    for(j=0;j<n;j++){
    temp[j]=a[i][j];
    }
__assume_aligned(a, 64);
__assume(i%16==0);
__assume(max_j%16==0);
__assume(n%16==0);
    for(j=0;j<n;j++){
    a[i][j]=a[max_j][j];
    }
__assume_aligned(a, 64);
__assume_aligned(temp, 64);
__assume(n%16==0);
    for(j=0;j<n;j++){
    a[max_j][j]=temp[j];
    }
  }
  
  _mm_free(temp);
}

 /*
 * lup_od_omp function
 * The key LU decomposition function.
 */
__attribute__ ((target(mic))) //offloadable
void lup_od_omp(int n, mat a, int min){
__assume_aligned(a, 64);
  int i,k, j, vec;
  k = min;
  //Ensuring numerical stability
  if(fabs(a[k][k])<=TINY){

  }else{
    
    //Intel Cilk array notation introduced
    //a[k+1:n-(k+1)][k] /= a[k][k];
    
   vec =  ((k+1)%16 == 0 && n-(k+1) >= 16 && (n-(k+1))%16 == 0);

    //OpenMP parallel for pragma introduced with static schedule for regular workload
#pragma omp parallel for default(none) shared(a,n,k,vec) private(i) schedule(static,4)
    for(i = k + 1; i < n; i++) {
      a[i][k] /= a[k][k];
      if(a[i][k]!=0){
	const double aik = a[i][k];
        __assume_aligned(a, 64);
       if(vec){
       __assume((k+1)%16==0);
       __assume((n-(k+1))%16==0);
       #pragma vector aligned
  	a[i][k+1:n-(k+1)] -= aik * a[k][k+1:n-(k+1)];	//Intel Cilk array notation
       }else{
        a[i][k+1:n-(k+1)] -= aik * a[k][k+1:n-(k+1)];   //Intel Cilk array notation
       }

      }
    }
  }
}
 
int main( int argc, char *argv[] ){
  int n = 3;
  int i, j;
  double start_time, end_time;
  double	*A_ptr, *read_ptr;
  mat A=NULL;
  
  n=read_arr(argv[1], &read_ptr);
  
  printf("Input ok: %d\n",n);
  
  //Initializing variables
  A_ptr=read_ptr;	
  A = (double**)_mm_malloc(sizeof(double*) * n,64);
  A[0] =A_ptr;
  for (i=0; i < n; i++)
    A[i] = A[0] + n * i;
  A_ptr=A[0];
  
  //Starting the timer
  start_time=omp_get_wtime();
  
  //Offloading to Intel MIC (Intel Xeon Phi) - this line was modified when executing on the host
#pragma offload target(mic) in(n) inout(A_ptr : length(n*n) alloc_if(1)) nocopy(end_time) nocopy(start_time)
  {
    int i, num_threads, j;
    double start_time, end_time;
#pragma omp parallel
    {
      num_threads=omp_get_num_threads();
    }	
    //Starting timer for "on Phi" timings
    start_time=omp_get_wtime();
    //Initializing local variables
    mat x = (double**)_mm_malloc(sizeof(double*) * n,64);

    x[0]=A_ptr;
    double *start;
    start = x[0]; 
    // calculation of the pivot
#pragma omp parallel for
    for (i=1; i < n; i++)
      x[i] = start + n * i;
    for(j=0; j < n - 1; j++){
      __assume_aligned(x, 64);
      __assume(j%16==0);
      __assume(n%16==0);
      mat_pivot(x, j, n);
      
      // invoking the actual LU decomposition
      __assume_aligned(x, 64);
      __assume(j%16==0);
      __assume(n%16==0);
      lup_od_omp(n, x, j);
    }	
    //Copying matrix back to A_ptr
    _mm_free(x);
    end_time=omp_get_wtime();
    printf("Running time on %d threads on MIC: %lf\n", num_threads, end_time-start_time);
  }
  end_time=omp_get_wtime();
  printf("Running time with offload to MIC: %lf\n", end_time-start_time);
  //Uncomment to show the matrix
  //_show(A); 
  
  //Cleaning up
  mat_del(A); 
  
  printf("\n");
  
  return 0;
}
