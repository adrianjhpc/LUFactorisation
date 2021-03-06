/*
 * Mateusz Iwo Dubaniowski, EPCC, 2016
 * Adrian Jackson, EPCC, 2016
 */
 
 /*
  * readmatrix.c
  * This file contains methods to read a matrix stored in matrix exchange format
  */

#include <stdio.h>
#include <stdlib.h>

#define ROW_LENGTH 1024
/*
 * process_header function
 * This function reads the header of a matrix exchange (.mtx) format file.
 */
static int process_header(FILE *f, int *m, int *n, int *mxn)
{
	char line[ROW_LENGTH];
	int read;
	line[0] = 0;
	*m = *n = *mxn = 0;

	do 
	{
		if (fgets(line, ROW_LENGTH-1, f) == NULL) return 1;
	} while (line[0] == '%'); // Skipping lines until we pass the header comments
	// read the header
#ifdef INT64__
#ifdef _WIN32
	if (sscanf(line, "%I64u %I64u %I64u", m, n, mxn) == 3)
#else
	if (sscanf(line, "%llu %llu %llu", m, n, mxn) == 3)
#endif
#else
	if (sscanf(line, "%u %u %u", m, n, mxn) == 3)
#endif
	{
		return 0;
	}
	else
	{
		do
		{ 
#ifdef INT64__
#ifdef _WIN32
			read = fscanf(f, "%I64u %I64u %I64u", m, n, mxn);
#else
			read = fscanf(f, "%llu %llu %llu", m, n, mxn);
#endif
#else
			read = fscanf(f, "%u %u %u", m, n, mxn);
#endif
			if (read == EOF) return 1;
		} while (read != 3);
	}

	return 0;
}

/*
 * read_arr function
 * This function reads the body of a matrix exchange (.mtx) format file.
 */
int read_arr(char* file, double **arr){
	int m, n, mxn, k;
	unsigned int i, j;
	double temp;

	//Open the file
	FILE *fp;
	printf("Start...\n");
	fp = fopen(file, "r");
	
	//Read the header
	process_header(fp, &m, &n, &mxn);
	*arr=calloc(m*n, sizeof(double));	//allocate the memory based on the header

	// Reading code line by line to populate the matrix
	for (k=0; k<mxn; ++k)
	{
		fscanf(fp, "%u %u %lf", &i, &j, &temp);
		
		(*arr)[(i-1)*n+(j-1)]=temp;
		
	}

	fclose(fp);
	printf("File read completed.\n");
	return n;
}

// Commented code below was used to test the functions in this file and their behaviour
/*
int main(){
	int n, i, j;
	double *arr;
	n=read_arr("Hamrle1.mtx", &arr);
	printf("n: %d, arr[0]: %.1lf\n", n, arr[0]);
	for(i=0; i<n; i++){
                for(j=0; j<n; j++)
                        printf("%.1lf ", arr[i*n+j]);
                printf("\n");
        }
	printf("Done!\n");
	return 0;
}
*/
