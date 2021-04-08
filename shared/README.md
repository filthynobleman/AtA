# Shared Memory AtA
In this project a parallel version of the AtA algorithm in a shared memory model is provided, together with a sequential, highly optimized implementation of the Strassen's algorithm.

## Prerequisite
In order to build the shared memory parallel implementation of the AtA algorithm, you need a working installation of *Intel MKL*, which you can either find [here](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html) or install through the package manager with the command
```
sudo apt-get install intel-mkl
```
As an alternative, you can build the project using any *CBLAS* backend. For example, you can use [*OpenBlas*](https://www.openblas.net/), or the distribution provided from the *Linux* package manager, which you can easily install with the command
```
sudo apt-get install libblas-dev
```

Also, you will need a working installation of [*OpenMP*](https://www.openmp.org/), which you can download from the website or install through the package manager with the command
```
sudo apt-get install libomp-dev
```

In case you want to build the project with *Intel MKL*, you will also need a working installation of *Threading Building Blocks* (*TBB*), which you can find [here](https://github.com/oneapi-src/oneTBB) as a compilable source, [here](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html) as an intel distribution, or install through the package manager with the command
```
sudo apt-get install tbb
```
> :warning: **Read this if you choose a TBB version different from the *Intel* one:** the `CMakeLists.txt` file for this project relies on the *CMake* module `FindTBB.cmake`, located in the directory `cmake`. This module actually works well with *Intel TBB* distribution, but you could need to modify it in case you choose a different distribution.

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
 - `WHICH_BLAS`=**MKL**|BLAS, which determines if the BLAS backend for the project is the *Intel MKL* implementation or the system default *BLAS* implementation
 - `ATA_GEMM`=**STRASSEN**|GEMM, which determines the matrix multiplication algorithm used by sequential ATA to perform standard `A^T * B` multiplication.
 - `ATA_MT_GEMM`=**STRASSEN**|GEMM, which determines the matrix multiplication algorithm used by multi-threaded ATA by threads performing standard matrix multiplication.
 - `ATA_MT_SYRK`=**ATA**|SYRK, which determines the matrix multiplication algorithm used by multi-threaded ATA by threads performing `A^T * A` multiplication.

Compile
```
make
```

## Running the examples
The building process produces, in the building directory, two executables
 - `TestAtAShared`
 - `TestSYRK`


### Testing shared memory AtA algorithm
The program `TestAtAShared` generates a random matrix, `A`, with dimensions `k` by `n`, and computes, using `p` threads, another matrix `C`, with dimension `n` by `n`, result of the operation `A^T * A`.

The command line execution for `TestAtA` is the following
```
./TestAtAShared [-N n] [-K k] [-P p] [-c] [-t]
```
where the options are:
 - `-N` to select the number of columns of `A` and the dimensions of `C`. Default is `1000`;
 - `-K` to select the number of rows of `A`. Default is `1000`;
 - `-P` to select the number of threads. Default is `16`;
 - `-c` to determine if the result of the operation must be tested against the result from the *BLAS* routine `?syrk`;
 - `-t` to determine if the execution time of the algorithm must be printed on screen. In combination with `-c`, also prints the execution time for *BLAS*.

### Testing multi-threaded SYRK algorithm
The program `TestSYRK` generates a random matrix, `A`, with dimensions `k` by `n`, and computes, using `p` threads, another matrix `C`, with dimension `n` by `n`, result of the operation `A^T * A`.

The command line execution for `TestSYRK` is the following
```
./TestSYRK [-N n] [-K k] [-P p] [-t]
```
where the options are:
 - `-N` to select the number of columns of `A` and the dimensions of `C`. Default is `1000`;
 - `-K` to select the number of rows of `A`. Default is `1000`;
 - `-P` to select the number of threads. Default is `16`;
 - `-t` to determine if the execution time of the algorithm must be printed on screen.