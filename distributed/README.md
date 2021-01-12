# Distributed Memory AtA
In this project a parallel version of the AtA algorithm in a distributed memory model is provided, together with a sequential, highly optimized implementation of the Strassen's algorithm.

## Prerequisite
In order to build the distributed memory parallel implementation of the AtA algorithm, you need a working installation of *Intel MKL*, which you can either find [here](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html) or install through the package manager with the command
```
sudo apt-get install intel-mkl
```
We plan to extend the codebase to also support any other implementaion of the [*CBLAS*](http://www.netlib.org/blas/) and [*ScaLAPACK*](http://www.netlib.org/scalapack/) interfaces, but for now the code completely relies on *Intel MKL*.

Also, you will need a working *MPI* installation. The [*Intel oneAPI Base toolkit*](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html) provides one (together with *MKL*), but you can also use (**recommended**) the [*OpenMPI*](https://www.open-mpi.org/) distribution, which can be either downloaded from the website or installed through a package manager with a command like
```
sudo apt-get install openmpi
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


Compile
```
make
```

## Running the examples
The building process produces, in the building directory, the executable
 - `TestAtADistributed`


### Testing distributed memory AtA algorithm
The program `TestAtADistributed` generates a random matrix, `A`, with dimensions `k` by `n`, and computes, using `p` threads, another matrix `C`, with dimension `n` by `n`, result of the operation `A^T * A`.

The command line execution for `TestAtADistributed` is the following
```
mpirun -np p ./TestAtADistributed [-N n] [-K k] [-c] [-t]
```
where the options are:
 - `-N` to select the number of columns of `A` and the dimensions of `C`. Default is `1000`;
 - `-K` to select the number of rows of `A`. Default is `1000`;
 - `-c` to determine if the result of the operation must be tested against the result from the *ScaLAPACK* routine `p?syrk`;
 - `-t` to determine if the execution time of the algorithm must be printed on screen. In combination with `-c`, also prints the execution time for *ScaLAPACK*.