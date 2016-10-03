CC=icc
CFLAGS=-qopenmp -O3 
MKLFLAGS= -L${MKLROOT}/lib/intel64/ -Wl,--start-group -lmkl_intel_lp64 \
-lmkl_core -lmkl_intel_thread -Wl,--end-group

ALL = cilk ivdep simd mkl
all: $(ALL)
.PHONY: all
KNCOFFLOAD =  knccilk kncivdep kncsimd kncmkl
kncoffload:  $(KNCOFFLOAD)
.PHONY: kncoffload

cilk: readmatrix.c ludecompCilk.c 
	$(CC) $(CFLAGS) -o $@ $^

ivdep: readmatrix.c ludecompIvdep.c
	$(CC) $(CFLAGS) -o $@ $^

simd: readmatrix.c ludecompSimd.c
	$(CC) $(CFLAGS) -o $@ $^

mkl: readmatrix.c ludecompMkl.c
	$(CC) $(CFLAGS) -o $@ $^ ${MKLFLAGS}
	
knccilk: readmatrix.c ludecompCilkOffloadKNC.c
	$(CC) $(CFLAGS) -o $@ $^

kncivdep: readmatrix.c ludecompIvdepOffloadKNC.c
	$(CC) $(CFLAGS) -o $@ $^ 

kncsimd: readmatrix.c ludecompSimdOffloadKNC.c
	$(CC) $(CFLAGS) -o $@ $^

kncmkl: readmatrix.c ludecompMklOffloadKNC.c
	$(CC) $(CFLAGS) -o $@ $^ $(MKLFLAGS)

populate: populate.c
	$(CC) -o $@ $<
	
clean:
	rm -fr $(ALL) $(KNCOFFLOAD) populate *.o
