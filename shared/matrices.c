/**
 * @file       simulexec.h
 *
 * @brief      Implementation of matrices.h
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#include "matrices.h"
#include "strassen.h"



void mat_gemtm(const Matrix *A, const Matrix* B, Matrix* C, const real alpha, const real beta)
{
	//cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans,
	// strassen(CblasTrans, CblasNoTrans,
	// 		   A->NumCols, B->NumCols, A->NumRows,
	// 		   alpha, A->A, A->NumCols, B->A, B->NumCols,
	// 		   beta, C->A, C->NumCols);
	StrassenATB(A->NumCols, B->NumCols, B->NumRows, 
				alpha, A->A, A->LDD,
				B->A, B->LDD,
				C->A, C->LDD);
}


void submat_gemtm(const SubMatrix *A, const SubMatrix *B, SubMatrix *C, const real alpha, const real beta)
{
	const real* SA = A->M->A + SUBMATIDX(A, 0, 0);
	const real* SB = B->M->A + SUBMATIDX(B, 0, 0);
	real* SC = C->M->A + SUBMATIDX(C, 0, 0);
	// cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans,
	// 		   A->NumCols, B->NumCols, A->NumRows,
	// 		   alpha, SA, A->M->NumCols, SB, B->M->NumCols,
	// 		   beta, SC, C->M->NumCols);

	StrassenATB(A->NumCols, B->NumCols, B->NumRows, 
				alpha, SA, A->M->LDD,
				SB, B->M->LDD,
				SC, C->M->LDD);
}



void mat_syrk(const Matrix *A, Matrix *C, const real alpha, const real beta)
{
	cblas_syrk(CblasRowMajor, CblasLower, CblasTrans,
			   C->NumRows, A->NumRows, alpha, A->A, A->LDD, beta, C->A, C->LDD);
}


void submat_syrk(const SubMatrix *A, SubMatrix *C, const real alpha, const real beta)
{
	const real* SA = A->M->A + SUBMATIDX(A, 0, 0);
	real* SC = C->M->A + SUBMATIDX(C, 0, 0);
	cblas_syrk(CblasRowMajor, CblasLower, CblasTrans,
			   C->NumRows, A->NumRows, alpha, SA, A->M->LDD, beta, SC, C->M->LDD);
}


real mat_err(const Matrix* A, const Matrix *B)
{
	if (A->NumRows != B->NumRows || A->NumCols != B->NumCols)
	{
		fprintf(stderr, "AtA::mat_err::MATRIX_SIZE: Incompatible matrix sizes (%lu-by-%lu and %lu-by-%lu).\n", A->NumRows, A->NumCols, B->NumRows, B->NumCols);
		return -1.0;
	}
	register integer i, j;
	real Result = 0;
	for (i = 0; i < A->NumRows; i++)
	{
		for (j = 0; j <= i; ++j)
		{
			real dij = MATELM(A, i, j) - MATSYMELM(B, i, j);
			Result += dij * dij;	
		}
	}
	return sqrt(Result);
}