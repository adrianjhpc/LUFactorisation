/*
 * Mateusz Iwo Dubaniowski, EPCC, 2015
 * Adrian Jackson, EPCC, 2016
 */
 
 /*
  * ludecompSimd.c
  * This file contains LU decomposition code using OpenMP4.0 simd pragma
  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define TINY 1e-20 
#define _dim int n
 
//Declaration of matrix type
typedef double **mat;
 
/*
 * mat_zero function
 * This function populates a matrix with zeroes.
 */
void mat_zero(mat x, int n) {
	int i, j;
	for (i=0; i < n; i++)
			for (j=0; j < n; j++)
				x[i][j] = 0;
}
 
/*
 * mat_new function
 * This function allocates a new matrix.
 */
mat mat_new(_dim){
  int i;
  mat x = (double**)malloc(sizeof(double*) * n);		//Allocate array of pointers
  x[0]  = (double *)malloc(sizeof(double) * n * n);	//Allocate data array to the first pointer of the pointer array
  
  //Assign the array of pointers
  for (i=0; i < n; i++)
    x[i] = x[0] + n * i;
  mat_zero(x, n);			//Populate with zeroes
  
  return x;
}
 
/*
 * mat_copy function
 * This function copies matrix and returns a copy.
 */
mat mat_copy(void *s, _dim){
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
void mat_del(mat x) { 
  free(x[0]);
  free(x);
}

 /*
 * mat_show function
 * Print a matrix to standard output.
 */
#define _QUOT(x) #x
#define QUOTE(x) _QUOT(x)
#define _show(a) printf(QUOTE(a)" =");mat_show(a, 0, n)
void mat_show(mat x, char *fmt, _dim){
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
void mat_pivot(mat a, int min, _dim){
  int i, j, k;
  double *temp;//, *temp2;
  temp=malloc(n*sizeof(double));	
  
  i = min;
  int max_j = i;
  for (j=i; j < n; j++)
    if (fabs(a[j][i]) > fabs(a[max_j][i])) max_j = j;
  
  //Performing pivoting, using simd pragma
  if (max_j != i){
#pragma omp simd
    for (k=0; k < n; k++) {
      temp[k]=a[i][k];
    }
#pragma omp simd
    for (k=0; k < n; k++) {
      a[i][k]=a[max_j][k];
    }
#pragma omp simd
    for (k=0; k < n; k++) {
      a[max_j][k]=temp[k];
    }
  }
  
  free(temp);
}

 /*
 * lup_od_omp function
 * The key LU decomposition function.
 */
void lup_od_omp(int n, mat a, int min){
  int i,k,j;

  k = min;
  //Ensuring numerical stability
  if(fabs(a[k][k])<=TINY){
  }else{
    
    //OpenMP parallel for pragma introduced with static schedule for regular workload
#pragma omp parallel for shared(a,n,k) private(i) schedule(static)
    for(i = k + 1; i < n; i++) {
      a[i][k] /= a[k][k];
      if(a[i][k]!=0){
	const double aik = a[i][k];
	//OpenMP 4.0's simd pragma used to aid vectorisation
#pragma omp simd
	for(j = k + 1; j < n; j++) {
	  a[i][j] -= aik * a[k][j];
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
  
  printf("Input ok\n");
  
  //Initializing variables
  A_ptr=read_ptr;	
  A = (double**)malloc(sizeof(double*) * n);
  A[0]=A_ptr;
  for (i=0; i < n; i++)
    A[i] = A[0] + n * i;
  A_ptr=A[0];
  
  //Starting the timer
  start_time=omp_get_wtime();
  int num_threads;	
#pragma omp parallel
  {
    num_threads=omp_get_num_threads();
  }
  //Starting timer for "on Phi" timings
  start_time=omp_get_wtime();
  //Initializing local variables
  mat x = (double**)malloc(sizeof(double*) * n);
  x[0]=A_ptr;
  
  //calculation of the pivot
#pragma omp parallel for
  for (i=0; i < n; i++)
    x[i] = x[0] + n * i;
  
  for(j=0; j < n - 1; j++){
    mat_pivot(x, j, n);
    
    // invoking the actual LU decomposition
    lup_od_omp(n, x, j);
  }
  end_time=omp_get_wtime();
  printf("Running time on %d threads on host: %lf\n", num_threads, end_time-start_time);
  A[0]=A_ptr;
  for (i=0; i < n; i++)
    A[i] = A[0] + n * i;
  
  //Uncomment to show the matrix
  //_show(A); 
  
  //Cleaning up
  mat_del(A);
  
  printf("\n");
  
  return 0;
}
