# LUFactorisation
A LU Factorisation implementation used to benchmark various directives and languages for performance programming (i.e. effective use of vectorisation hardware).  The respository contains code for both standard CPUs and Intel's Xeon Phi KNC, which has an offload mode.

The code is designed to read in a matrix, parallelise the work of performing an LU factorisation on that matrix using OpenMP, then optimise the work of the indivdual threads using a variety of different options.

## Authors
The code in this repository was written by Adrian Jackson, 2016, based on some work done by Mateusz Iwo Dubaniowski as part of his dissertation on the MSc in HPC at EPCC in 2015.
  

## Building the code
To compile the code please use Makefile:

```
make <option>
```

The following options are available:
* `all` - compiles the targets below
  *    `cilk` - compiles an executable with Cilk array notation
  *       `simd` - compiles an executable with simd OpenMP pragmas
  *      `ivdep` - compiles ane executable with ivdep pragmas
  *       `mkl`  - compiles a version of the benchmark that uses MKL to do the LU factorisation rather than an of our computational kernels.  This requires either MKL to be installed on the system you are compiling on, or another implementation of LAPACK.
* `kncoffload` - compiles the same implementations as above, but with Intel LEO offloading functionality to use the Xeon Phi KNC as an offload device.  As this requires the KNC software enviroment (MPSS) to be installed this is not done by default
  *    `knccilk`
  *       `kncsimd` 
  *       `kncivdep` 
  *       `kncmkl`
*       `populate` - creates an executable that will generate a dense matrix for benchmarking
*       `clean` - cleans the directory of executables and `.o` files

make without an option builds `cilk`, `simd`, `ivdep`, and `mkl`.

## Running the code
Compiled code is run the following manner:
```
./codename <matrix>
```
Where codename is replaced by `cilk`, `simd`, `ivdep`, or `mkl` respectively.
`<matrix>`, the input matrix, is expected to be in [matrix market format](http://math.nist.gov/MatrixMarket/formats.html), that is a file that has each entry in the matrix on a separate line in the format row, column, data value.

The code will run with sparse or dense matrices.  For our benchmarks we used a range of sparse matrices from the [University of Florida Sparse Matrix Collection] (http://www.cise.ufl.edu/research/sparse/matrices/index.html) and we generated our own dense matrices.

## Generating dense matrices
The `populate` executable will generate square dense matrices.  You can use it as follows:
```
./populate X > matrixfilename.mtx
```

Where `X` is an integer which specifies the size of one side of the matrix (i.e. for a 10x10 matrix, `X` would be replaced by 10 above)

## License
LU Factorisation benchmark
Copyright (C) 2015 Mateusz Iwo Dubaniowski, 2016 Adrian Jackson

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details (LICENSE.txt).

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

