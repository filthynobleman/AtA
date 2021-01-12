#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <mpi.h>
#include <string.h>
#include <unistd.h>
#include <mkl.h>
#include <mkl_blas.h>
#include <mkl_cblas.h>
#include <mkl_scalapack.h>
#include <mkl_pblas.h>
#include <mkl_blacs.h>
#include <mkl_lapacke.h>
#include <math.h>
#include <ctype.h>
#include "matrices.h"
#include "common.h"

void GenerateRandomMatrix(real* M, size_t n, size_t m)
{
	register size_t i;
	for (i = 0; i < n * m; i++)
		M[i] = rand() / ((real)RAND_MAX);
}


int main(int argc, const char** argv){


	size_t NumRuns;
	size_t N = 1000;
	size_t M = 1000;
	srand(0);

	
	if (argc < 5)
	{
		fprintf(stderr, "Not enough input arguments.\n");
		printf("Usage: AtA [N M]\n");
		exit(1);
	}

	NumRuns = atol(argv[1]);
	N = atol(argv[2]);
	M = atol(argv[3]);
	srand(0);

	double *Times = (double*)malloc(NumRuns * sizeof(double));
	double *CompTimes = (double*)malloc(NumRuns * sizeof(double));
	MPI_Init(NULL, NULL);
	int id, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	double t1, t2, t3, diff;

	real* A = (real*)malloc(N * M * sizeof(real));
	if (A == NULL)
	{
		fprintf(stderr, "Cannot allocate A->");
		exit(1);
	}
	real *C = (real*)malloc(M * M * sizeof(real));

	for (int i = 0; i < NumRuns; ++i){

		if (id == 0){
			GenerateRandomMatrix(A, N, M);
			memset(C, 0, M * M * sizeof(real));
		} else if (i == 0){
			free(A);
			A = NULL;
			free(C);
			C = NULL;
		}
		

		// long long bhandle = sys2blacs_handle(MPI_COMM_WORLD);
	 //    long long icontxt = bhandle;
		long long icontxt;
		long long what = 0;
	    long long npcol;// 	qui npcol = 1;
	    long long nprow;//	qui nprow = size;
	    long long myrow, mycol;
	    int dims[] = {0, 0};

	    MPI_Dims_create(size, 2, dims);
	    nprow = dims[0];
	    npcol = dims[1];
	//    fprintf(deb, "grid is %d x %d\n", nprow, npcol);
	    long long negone = -1;
	    long long one = 1;
	    long long zero = 0;
	    long long n = (long long)N;
	    long long m = (long long)M;
	    long long k = (long long)M;
	    MKL_INT descA[9], descC[9], descSubA[9], descSubC[9];
	    


	    blacs_get(&negone, &what, &icontxt);
	    blacs_gridinit( &icontxt, "R", &nprow, &npcol );
	    blacs_gridinfo( &icontxt, &nprow, &npcol, &myrow, &mycol );
	    MKL_INT np, mp, kp, nb = 1, mb = 1, kb = 1;

	    np = numroc_( &n, &nb, &myrow, &zero, &nprow );
	    mp = numroc_( &m, &nb, &myrow, &zero, &nprow );
	    kp = numroc_( &k, &kb, &mycol, &zero, &npcol );
	    long long lldA = n;
	    long long lldC = m;
	    long long lldSubA = np;
	    long long lldSubC = kp;
	    MKL_INT info;
	    real dzero = 0;
	    real done = 1;
	    real *subA = (real*)calloc(np * mp, sizeof(real));
	    real *subC = (real*)malloc(mp * kp * sizeof(real));


		descinit_(descA, &n, &m, &n, &m, &zero, &zero, &icontxt, &lldA, &info);    
		descinit_(descC, &m, &m, &m, &m, &zero, &zero, &icontxt, &lldC, &info); 

        descinit_(descSubA, &n, &m, &nb, &mb, &zero, &zero, &icontxt, &lldSubA, &info);    
		descinit_(descSubC, &m, &m, &mb, &kb, &zero, &zero, &icontxt, &lldSubC, &info); 

		pdgeadd_("N", &n, &m, &done, A, &one, &one, descA, &dzero, subA, &one, &one, descSubA);
		t1 = MPI_Wtime();
		pdsyrk_("U", "N", &m, &n, &done, subA, &one, &one, descSubA, &dzero, subC, &one, &one, descSubC);
		t3 = MPI_Wtime();
		
		
		pdgeadd_("N", &m, &m, &done, subC, &one, &one, descSubC, &dzero, C, &one, &one, descC);
		
		t2 = MPI_Wtime();

		free(subA);
		free(subC);
		double diff1 = t3 - t1; // time for computation
		double diff2 = t2 - t3; // time for matrix retrievement 
		double diff3 = t2 - t1; // total time
		Times[i] = diff3;
		CompTimes[i] = diff1;

	    blacs_freebuff(&icontxt, &zero);
	    blacs_gridexit(&icontxt);
	}
	long long one = 1;
	blacs_exit(&one);
	if (id == 0){
		/*
		real *C2 = (real*)calloc(M * M, sizeof(real));
		cblas_dsyrk(CblasRowMajor, CblasLower, CblasTrans, N, M, 1, A, M, 0, C2, M);
		real error = 0;
		for (int i = 0; i  < M; i++){
			for (int j = 0; j <= i; ++j)
				error += (C2[i * M + j] - C[i * M + j]) * (C2[i * M + j] - C[i * M + j]);
		}
		printf("error = %f\n", sqrt(error));
		*/
		FILE *stream = fopen(argv[4], "w");
		if (stream == NULL){
			fprintf(stderr, "Stream is NULL\n");
			exit(1);
		}
		for (int i = 0; i < NumRuns; ++i){
			fprintf(stream, "%f,%f\n", Times[i], CompTimes[i]);
		}
		fclose(stream);
	}
	MPI_Finalize();
    return 0;
}
