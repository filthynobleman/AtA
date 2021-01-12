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


#define ATA_BASE_CASE   2048

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
    cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, 
               N - N2, N2, K2, alpha, A12, lda, A11, lda, 1, C21, ldc);
    cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, 
               N - N2, N2, K - K2, alpha, A22, lda, A21, lda, 1, C21, ldc);
#endif
}