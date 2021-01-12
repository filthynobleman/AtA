#ifndef ATA_COMMON_H_
#define ATA_COMMON_H_


#include <stdlib.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


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
// 					 	const int, const int, const real, const real*, const int, const real, real*, const int);
// typedef void (*gemmfun)(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE, const enum CBLAS_TRANSPOSE,
// 					 	const int, const int, const int, const real, const real*, const int, const real*, const int,
// 					 	const real, real*, const int);


#ifdef DOUBLE_PRECISION
// const syrkfun cblas_syrk = &cblas_dsyrk;
// const gemmfun cblas_gemm = &cblas_dgemm;
#define cblas_syrk cblas_dsyrk
#define cblas_gemm cblas_dgemm
#define cblas_axpy cblas_daxpy
#define cblas_scal cblas_dscal
#define MPI_RREAL MPI_DOUBLE
#else
// const syrkfun cblas_syrk = &cblas_ssyrk;
// const gemmfun cblas_gemm = &cblas_sgemm;
#define cblas_syrk cblas_ssyrk
#define cblas_gemm cblas_sgemm
#define cblas_axpy cblas_saxpy
#define cblas_scal cblas_sscal
#define MPI_RREAL MPI_FLOAT
#endif




#endif // ATA_COMMON_H_
