/*
 * Mateusz Iwo Dubaniowski, EPCC, 2015
 * Adrian Jackson, EPCC, 2016
 */
  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mkl.h>

#define TINY 1e-20 
#define _dim int n
 
//Declaration of matrix type
typedef double **mat;
   
/*
 * mat_del function
 * Delete an array and free memory.
 */
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

int main( int argc, char *argv[] )
{
  int n = 3;
  int i, j, k;
  double start_time, end_time;
  double	*A_ptr, *read_ptr;
  int *ipiv;
  int info;
  mat A=NULL;
  double *temparr;

  n=read_arr(argv[1], &read_ptr);
  
  printf("Input ok\n");
  
  //Initializing variables
  A_ptr=read_ptr;	
  A = (double**)_mm_malloc(sizeof(double*) * n,64);
  A[0] =A_ptr;
  for (i=0; i < n; i++)
    A[i] = A[0] + n * i;
  A_ptr=A[0];
  
  
  ipiv = (int *)malloc(sizeof(int) * n);
  temparr = (double *)malloc(sizeof(double*) * n * n);
  
  start_time=omp_get_wtime();	
  k = 0;
  for(i=0; i < n; i++){
    for(j=0; j < n; j++){
      temparr[k] = A[i][j];
      k++;
    }
  }
  
  int num_threads;
#pragma omp parallel
  {
    num_threads=omp_get_num_threads();
  }	
  
  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,n,n,temparr,n,ipiv);
  
  k = 0;
  for(i=0; i < n; i++){
    for(j=0; j < n; j++){
      A[i][j] = temparr[k];
      k++;
    }
  }
  end_time=omp_get_wtime();

  printf("Running time on %d threads on host and phi: %lf\n", num_threads, end_time-start_time);
  
  //Uncomment to show the matrix
  //_show(A); 
  
  //Cleaning up
  mat_del(A); 
  
  free(ipiv);
  free(temparr);
  
  printf("\n");

  return 0;
}
