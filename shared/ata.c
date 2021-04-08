/**
 * @file       ata.c
 *
 * @brief      Implementation of ATA algorithm.
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#include "ata.h"
#include "strassen.h"
#include "simulexec.h"

#include <omp.h>

#define ATA_BASE_CASE   512
#define ALPHA 0.5

// Computes C = C + alpha * A^T * A.
// A is K-by-N, C is N-by-N
void AtA(const MKL_INT N, const MKL_INT K, 
         const real alpha, const real* A, const MKL_INT lda,
         real* C, const MKL_INT ldc)
{
    if (N <= ATA_BASE_CASE)
    {
        cblas_syrk(CblasRowMajor, CblasLower, CblasTrans,
                   N, K, alpha, A, lda, 1, C, ldc);
        return;
    }


    // Half sizes
    const MKL_INT N2 = N / 2;
    const MKL_INT K2 = K / 2;

    // Create submatrices of A
    const real* A11 = A;
    const real* A12 = A + N2;
    const real* A21 = A + K2 * lda;
    const real* A22 = A + K2 * lda + N2;
    // Create submatrices of C
    real* C11 = C;
    real* C12 = C + N2;
    real* C21 = C + N2 * ldc;
    real* C22 = C + N2 * ldc + N2;


    // C11 = A11^T * A11 + A21^T * A21
    AtA(N2, K2, alpha, A11, lda, C11, ldc);
    AtA(N2, K - K2, alpha, A21, lda, C11, ldc);

    // C22 = A12^T * A12 + A22^T * A22
    AtA(N - N2, K2, alpha, A12, lda, C22, ldc);
    AtA(N - N2, K - K2, alpha, A22, lda, C22, ldc);

    // C21 = A12^T * A11 + A22^T * A21
#ifdef ATA_USE_STRASSEN
    StrassenATB(N - N2, N2, K2, alpha, A12, lda, A11, lda, C21, ldc);
    StrassenATB(N - N2, N2, K - K2, alpha, A22, lda, A21, lda, C21, ldc);
#else
    cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, N - N2, N2, K2, alpha, A12, lda, A11, lda, 1, C21, ldc);
    cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, N - N2, N2, K - K2, alpha, A22, lda, A21, lda, 1, C21, ldc);
#endif
}





typedef struct AtAParams
{
	SubMatrix* SA;
	SubMatrix* SC;
	integer    id;
	integer    num_threads;
} AtAParams;


void AtALocal_MT(AtAParams* params)
{
	SubMatrix* A = params->SA;
	SubMatrix* C = params->SC;
	integer myid = params->id;
	integer numth = params->num_threads;


	TaskTree tree = SimulateExecution(A->M, C->M, numth, ALPHA);
	TaskNode* Node;
	Task T;
	size_t Level = tree.Height;
	do
	{
		GetTaskNode(&tree, myid, Level--, &Node);
	} while(Node == NULL);
	T = Node->T;
	if (T.Type == ATA)
	{
		real* rA = T.A.M->A + SUBMATIDX(&(T.A), 0, 0);
		real* rC = T.C.M->A + SUBMATIDX(&(T.C), 0, 0);
		MKL_INT N = T.A.NumCols;
		MKL_INT K = T.A.NumRows;
		MKL_INT lda = T.A.M->NumCols;
		MKL_INT ldc = T.C.M->NumCols;
#ifdef ATA_MT_USE_ATA
		AtA(N, K, 1, rA, lda, rC, ldc);
#else
        cblas_syrk(CblasRowMajor, CblasLower, CblasTrans, N, K, 1, rA, lda, 1, rC, ldc);
#endif
	}
	else
	{
		real* rA = T.A.M->A + SUBMATIDX(&(T.A), 0, 0);
		real* rB = T.B.M->A + SUBMATIDX(&(T.B), 0, 0);
		real* rC = T.C.M->A + SUBMATIDX(&(T.C), 0, 0);
		MKL_INT M = T.A.NumCols;
		MKL_INT N = T.B.NumCols;
		MKL_INT K = T.A.NumRows;
		MKL_INT lda = T.A.M->NumCols;
		MKL_INT ldb = T.B.M->NumCols;
		MKL_INT ldc = T.C.M->NumCols;
#ifdef ATA_MT_USE_STRASSEN
		StrassenATB(M, N, K, 1, rA, lda, rB, ldb, rC, ldc);
#else
        cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1, rA, lda, rB, ldb, 1, rC, ldc);
#endif
	}
}


void AtA_MT(const Matrix* A, Matrix *C, integer NumThreads)
{
	if (C == NULL)
	{
		fprintf(stderr, "Ata: C is NULL.");
		return;
	}
	// if (C->A == NULL)
	// 	C->A = (real *)calloc(A->NumCols * A->NumCols, sizeof(real));
	// else
	// 	memset(C->A, 0, A->NumCols * A->NumCols * sizeof(real));
	if (C->A == NULL)
	{
		fprintf(stderr, "AtA: Cannot allocate C.");
		return;
	}
	C->NumRows = A->NumCols;
	C->NumCols = A->NumCols;

	const SubMatrix A0 = {A, 0, 0, A->NumRows, A->NumCols};
	SubMatrix C0 = {C, 0, 0, C->NumRows, C->NumCols};
	integer id;
	AtAParams* ataps = (AtAParams*)malloc(NumThreads * sizeof(AtAParams));
	
	omp_set_dynamic(0);
	omp_set_num_threads(NumThreads);
	#pragma omp parallel
	#pragma omp for
	for (id = 0; id < NumThreads; ++id)
	{
		ataps[id].SA = &A0;
		ataps[id].SC = &C0;
		ataps[id].id = id;
		ataps[id].num_threads = NumThreads;
		AtALocal_MT(ataps + id);
	}

	free(ataps);
}
