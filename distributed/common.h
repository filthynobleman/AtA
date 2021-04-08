#ifndef ATA_COMMON_H_
#define ATA_COMMON_H_

#define AtA_Distributed_VERSION_MAJOR 
#define AtA_Distributed_VERSION_MINOR 
#define DOUBLE_PRECISION
/* #undef SYNCHRONIZE_PROCESSES */
/* #undef ONLY_COMPUTE */


#include <stdlib.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#define TimerStart(beg)     gettimeofday(&(beg), NULL)
#define TimerStop(end)      gettimeofday(&(end), NULL);
#define TimerGet(beg, end)  (((end).tv_sec - (beg).tv_sec) + ((end).tv_usec - (beg).tv_usec) / 1.0e6)


#ifdef LONG_INTEGERS
typedef unsigned long integer;
#else
typedef unsigned int integer;
#endif

#ifdef DOUBLE_PRECISION
typedef double real;
#else
typedef float real;
#endif


// typedef void (*syrkfun)(const enum CBLAS_ORDER, const enum CBLAS_UPLO, const enum CBLAS_TRANSPOSE,
//                      const int, const int, const real, const real*, const int, const real, real*, const int);
// typedef void (*gemmfun)(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const enum CBLAS_TRANSPOSE,
//                      const int, const int, const int, const real, const real*, const int, const real*, const int,
//                      const real, real*, const int);


#ifdef DOUBLE_PRECISION
// const syrkfun cblas_syrk = &cblas_dsyrk;
// const gemmfun cblas_gemm = &cblas_dgemm;
#define cblas_syrk cblas_dsyrk
#define cblas_gemm cblas_dgemm
#define cblas_axpy cblas_daxpy
#define cblas_scal cblas_dscal
#define cblas_dot  cblas_ddot
#define cblas_gemv cblas_dgemv
#define MPI_RREAL MPI_DOUBLE
#else
// const syrkfun cblas_syrk = &cblas_ssyrk;
// const gemmfun cblas_gemm = &cblas_sgemm;
#define cblas_syrk cblas_ssyrk
#define cblas_gemm cblas_sgemm
#define cblas_axpy cblas_saxpy
#define cblas_scal cblas_sscal
#define cblas_dot  cblas_sdot
#define cblas_gemv cblas_sgemv
#define MPI_RREAL MPI_FLOAT
#endif




#endif // ATA_COMMON_H_
