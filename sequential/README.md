# Sequential AtA
In this project a sequential version of the AtA algorithm is provided, together with a sequential, highly optimized implementation of the Strassen's algorithm.

## Prerequisite
In order to build the sequential implementation of the AtA algorithm, you need a working installation of *Intel MKL* (which you can find [here](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html)).  
As an alternative, you can build the project using any *CBLAS* backend. For example, you can use [OpenBlas](https://www.openblas.net/), or the distribution provided from the *Linux* package manager, which you can easily install with the command
```
`sudo apt-get install libblas-dev
```

## Building the project
Create a directory for building the project and move in
```
mkdir build
cd build
```

Configure using CMake
```
cmake <options> ..
```
where available options are:
 - `USE_DOUBLES`=**ON**|OFF, which determines if the project must be compiled using double or single floating point precision.
 - `WHICH_BLAS`=**MKL**|BLAS, which determined if the BLAS backend for the project is the *Intel MKL* implementation or the system default *BLAS* implementation

Compile
```
make
```

## Running the examples
The building process produces, in the building directory, two executables
 - `TestStrassen`
 - `TestAtA`

### Testing Strassen's algorithm
The program `TestStrassen` generates two random matrices, `A` and `B`, respectively with dimensions `k` by `m` and `k` by `n`, and computes a third matrix `C`, with dimension `m` by `n`, result of the operation `A^T * B`.

The command line execution for `TestStrassen` is the following
```
./TestStrassen [-M m] [-N n] [-K k] [-c] [-t]
```
where the options are:
 - `-M` to select the number of columns of `A` and rows of `C`. Default is `1000`;
 - `-N` to select the number of columns of `B` and `C`. Default is `1000`;
 - `-K` to select the number of rows of `A` and `B`. Default is `1000`;
 - `-c` to determine if the result of the operation must be tested against the result from the *BLAS* routine `?gemm`;
 - `-t` to determine if the execution time of the algorithm must be printed on screen. In combination with `-c`, also prints the execution time for *BLAS*.


### Testing AtA algorithm
The program `TestAtA` generates a random matrix, `A`, with dimensions `k` by `n`, and computes another matrix `C`, with dimension `n` by `n`, result of the operation `A^T * A`.

The command line execution for `TestAtA` is the following
```
./TestAtA [-N n] [-K k] [-c] [-t]
```
where the options are:
 - `-N` to select the number of columns of `A` and the dimensions of `C`. Default is `1000`;
 - `-K` to select the number of rows of `A`. Default is `1000`;
 - `-c` to determine if the result of the operation must be tested against the result from the *BLAS* routine `?syrk`;
 - `-t` to determine if the execution time of the algorithm must be printed on screen. In combination with `-c`, also prints the execution time for *BLAS*.