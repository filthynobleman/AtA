#include "matrices.h"
//#include "strassen.h"


void mat_gemtm(const Matrix *A, const Matrix* B, Matrix* C, const real alpha, const real beta)
{
	cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans,
//	strassen(CblasTrans, CblasNoTrans, 
			   A->NumCols, B->NumCols, A->NumRows,
			   alpha, A->A, A->NumCols, B->A, B->NumCols,
			   beta, C->A, C->NumCols);
}


void submat_gemtm(const SubMatrix *A, const SubMatrix *B, SubMatrix *C, const real alpha, const real beta)
{
	const real* SA = A->M->A + SUBMATIDX(A, 0, 0);
	const real* SB = B->M->A + SUBMATIDX(B, 0, 0);
	real* SC = C->M->A + SUBMATIDX(C, 0, 0);
	cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans,
			   A->NumCols, B->NumCols, A->NumRows,
			   alpha, SA, A->M->NumCols, SB, B->M->NumCols,
			   beta, SC, C->M->NumCols);
}



void mat_syrk(const Matrix *A, Matrix *C, const real alpha, const real beta)
{
	cblas_syrk(CblasRowMajor, CblasLower, CblasTrans,
			   C->NumRows, A->NumRows, alpha, A->A, A->NumCols, beta, C->A, C->NumCols);
}


void submat_syrk(const SubMatrix *A, SubMatrix *C, const real alpha, const real beta)
{
	const real* SA = A->M->A + SUBMATIDX(A, 0, 0);
	real* SC = C->M->A + SUBMATIDX(C, 0, 0);
	cblas_syrk(CblasRowMajor, CblasLower, CblasTrans,
			   C->NumRows, A->NumRows, alpha, SA, A->M->NumCols, beta, SC, C->M->NumCols);
}
void fprint_sym_mat(real *A, int n, FILE *file){
#ifdef DEBUG_MODE
	int count = 0;
	for (int i = 0; i < n; ++i){
		for (int j = 0; j <= i; ++j){
			fprintf(file, "%.3e\t", A[count++]);
		}
		fprintf(file, "\n");
	}
#endif
}
void print_mat(real *A, int noRowsToPrint, int noColsToPrint, int n, int m){
#ifdef DEBUG_MODE
	for (int i = 0; i < noRowsToPrint; ++i){
		for (int j = 0; j < noColsToPrint; ++j){
			printf("%.3e\t", A[i * m + j]);
		}
		printf("\n");
	}
#endif	
}

void fprint_mat(real *A, int noRowsToPrint, int noColsToPrint, int n, int m, FILE *deb){
#ifdef DEBUG_MODE
	for (int i = 0; i < noRowsToPrint; ++i){
		for (int j = 0; j < noColsToPrint; ++j){
			fprintf(deb, "%.3e\t", A[i * m + j]);
		}
		fprintf(deb, "\n");
	}	
#endif
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
