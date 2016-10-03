CC=icc
CFLAGS=-qopenmp -O3 
MKLFLAGS= -L${MKLROOT}/lib/intel64/ -Wl,--start-group -lmkl_intel_lp64 \
-lmkl_core -lmkl_intel_thread -Wl,--end-group

all: cilk ivdep simd mkl

cilk: readmatrix.c ludecompCilk.c 
	$(CC) $(CFLAGS) -o $@ readmatrix.c ludecompCilk.c

ivdep: readmatrix.c ludecompIvdep.c
	$(CC) $(CFLAGS) -o $@ readmatrix.c ludecompIvdep.c

simd: readmatrix.c ludecompSimd.c
	$(CC) $(CFLAGS) -o $@ readmatrix.c ludecompSimd.c

mkl: readmatrix.c ludecompMkl.c
	$(CC) $(CFLAGS) -o $@ readmatrix.c ludecompMkl.c $(MKLFLAGS)
	
populate: populate.c
	$(CC) -o populate populate.c
	
clean:
	rm -fr mkl simd cilk ivdep populate *.o
