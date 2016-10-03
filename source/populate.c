/*
 * Mateusz Iwo Dubaniowski, EPCC, 2015
 * Adrian Jackson, EPCC, 2016
 */
 
 /*
  * populate.c
  * This file defined a program used to create a randomised dense matrices for benchmarking purposes.
  * The program expects a single integer as input, on the command line, that defines the size of the 
  * dense matrix (we are generating square matrices so the input integer defines the size of one side of
  * the matrix).
  * The file prints the matrix values to screen, it is designed to be run on the command line and the 
  * output piped to a file to store the matrix data, i.e.:
  * ./populate 100 > dense100x100.mtx
  * The matrix is outputted in matrix market format, i.e. row, column, data value
  */
  
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + ((double)rand() / div);
}

int main(int argc, char *argv[])
{
   int i, j, n;
   time_t t;

   if(argc != 2){
      printf("This program is expecting a single command line option, the size of the matrix to be created.  This should be a single number (we are generating square matrices) that is the size of one side of the matrix.\n");
      printf("Aborting now because matrix size was not provided\n");
      return(1);
   }else{
      sscanf (argv[1],"%d",&n);
      if(n < 1){
         printf("Invalid size of matrix specified, expecting a number greater than zero but got %d\n",n);
         return(1);
      }
   } 
   
   /* Intializes random number generator */
   srand((unsigned) time(&t));

   /* Print 5 random numbers from 0 to 50 */
   for( i = 1 ; i <= n ; i++ ) 
   for( j = 1 ; j <= n ; j++ )
   	{
	      printf("%d %d %.5lf\n", i, j, randfrom(0, 100));
	}
   
   return(0);
}
