#ifndef ATA_MATRICES_H_
#define ATA_MATRICES_H_


#include "common.h"


struct Matrix
{
	real *A;
	integer NumRows;
	integer NumCols;
};
typedef struct Matrix Matrix;
typedef Matrix mat_t;


struct SubMatrix
{
	Matrix* M;
	integer BegRow;
	integer BegCol;
	integer NumRows;
	integer NumCols;
};
typedef struct SubMatrix SubMatrix;
typedef SubMatrix submat_t;




#define MATIDX(M, i, j)		        ((i) * (M)->NumCols + (j))
#define MATSYMIDX_(M, i, j)         ((i) * ((i) + 1) / 2 + (j))
#define MATSYMIDX(M, i, j)          ((i) >= (j) ? MATSYMIDX_(M, i, j) : MATSYMIDX_(M, j, i))
#define MATELM(M, i, j)		        (((M)->A)[MATIDX(M, i, j)])
#define MATSYMELM(M, i, j)          (((M)->A)[MATSYMIDX(M, i, j)])

#define SUBMATIDX(SM, i, j)	        MATIDX((SM)->M, (SM)->BegRow + (i), (SM)->BegCol + (j))
#define SUBMATSYMIDX(SM, i, j)      MATSYMIDX((SM)->M, (SM)->BegRow + (i), (SM)->BegCol + (j))
#define SUBMATELM(SM, i, j)         MATELM((SM)->M, (SM)->BegRow + (i), (SM)->BegCol + (j))
#define SUBMATSYMELM(SM, i, j)      MATSYMELM((SM)->M, (SM)->BegRow + (i), (SM)->BegCol + (j))




void mat_gemtm(const Matrix *A, const Matrix* B, Matrix* C, const real alpha, const real beta);
void submat_gemtm(const SubMatrix *A, const SubMatrix *B, SubMatrix *C, const real alpha, const real beta);

void mat_syrk(const Matrix *A, Matrix *C, const real alpha, const real beta);
void submat_syrk(const SubMatrix *A, SubMatrix *C, const real alpha, const real beta);

real mat_err(const Matrix* A, const Matrix *B);
void fprint_sym_mat(real *A, int n, FILE *file);
void print_mat(real *A, int noRowsToPrint, int noColsToPrint, int n, int m);

void fprint_mat(real *A, int noRowsToPrint, int noColsToPrint, int n, int m, FILE *deb);





#endif // ATA_MATRICES_H_