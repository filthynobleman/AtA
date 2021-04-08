#include "strassen.h"

#define ATB_BASE 512



void StrassenATBRec(const MKL_INT m, const MKL_INT n, const MKL_INT k,
                    const real alpha, const real* A, const MKL_INT lda,
                    const real* B, const MKL_INT ldb,
                    real* C, const MKL_INT ldc,
                    real* M, real* lh, real* rh)
{
    if (k <= ATB_BASE || n <= ATB_BASE || m <= ATB_BASE)
    {
        cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                   m, n, k, alpha, A, lda, B, ldb, 0, C, ldc);
        return;
    }


    // Recursive case
    // Get the ceiling of half the dimensions
    MKL_INT m2 = (m + 1) / 2;
    MKL_INT n2 = (n + 1) / 2;
    MKL_INT k2 = (k + 1) / 2;


    // Compute the submatrices to use as working matrices in the recursive calls
    real* sM = M + m2 * n2;
    real* slh = lh + k2 * m2;
    real* srh = rh + k2 * n2;

    // Compute the input submatrices for the recursive steps
    const real* A11 = A;
    const real* A12 = A + m2;
    const real* A21 = A + k2 * lda;
    const real* A22 = A + k2 * lda + m2;
    const real* B11 = B;
    const real* B12 = B + n2;
    const real* B21 = B + k2 * ldb;
    const real* B22 = B + k2 * ldb + n2;
    real* C11 = C;
    real* C12 = C + n2;
    real* C21 = C + m2 * ldc;
    real* C22 = C + m2 * ldc + n2;


    // C11 = M1 + M4 - M5 + M7
    // C12 = M3 + M5
    // C21 = M2 + M4
    // C22 = M1 - M2 + M3 + M6

    register MKL_INT i;
    MKL_INT lhlen;
    MKL_INT rhlen;
    // M1 = (A11 + A22)^T * (B11 + B22)
    // lh = A11 + A22 is k2-by-m2
    // rh = B11 + B22 is k2-by-n2
    lhlen = m2 * sizeof(real);
    rhlen = n2 * sizeof(real);
    for (i = 0; i < k - k2; ++i)
    {
        memcpy(lh + i * m2, A11 + i * lda, lhlen);
        cblas_axpy(m - m2, 1, A22 + i * lda, 1, lh + i * m2, 1);

        memcpy(rh + i * n2, B11 + i * ldb, rhlen);
        cblas_axpy(n - n2, 1, B22 + i * ldb, 1, rh + i * n2, 1);
    }
    if (k2 > k - k2)
    {
        memcpy(lh + (k2 - 1) * m2, A11 + (k2 - 1) * lda, lhlen);
        memcpy(rh + (k2 - 1) * n2, B11 + (k2 - 1) * ldb, rhlen);
    }
    // Recursive call
    memset(M, 0, m2 * n2 * sizeof(real));
    StrassenATBRec(m2, n2, k2, alpha, lh, m2, rh, n2, M, n2, sM, slh, srh);
    // print_mat(M, m2, n2);
    // Update C11 and C22
    for (i = 0; i < m - m2; ++i)
    {
        cblas_axpy(n2, 1, M + i * n2, 1, C11 + i * ldc, 1);
        cblas_axpy(n - n2, 1, M + i * n2, 1, C22 + i * ldc, 1);
    }
    if (m2 > m - m2)
        cblas_axpy(n2, 1, M + (m2 - 1) * n2, 1, C11 + (m2 - 1) * ldc, 1);

    // M2 = (A12 + A22)^T * B11
    // lh = A12 + A22 is k2-by-(m - m2)
    lhlen = (m - m2) * sizeof(real);
    for (i = 0; i < k - k2; ++i)
    {
        memcpy(lh + i * (m - m2), A12 + i * lda, lhlen);
        cblas_axpy(m - m2, 1, A22 + i * lda, 1, lh + i * (m - m2), 1);
    }
    if (k2 > k - k2)
        memcpy(lh + (k2 - 1) * (m - m2), A12 + (k2 - 1) * lda, lhlen);
    // Recurseve call
    memset(M, 0, m2 * n2 * sizeof(real));
    StrassenATBRec(m - m2, n2, k2, alpha, lh, m - m2, B11, ldb, M, n2, sM, slh, srh);
    // print_mat(M, m - m2, n2);
    // Update C21 and C22(-)
    for (i = 0; i < m - m2; ++i)
    {
        cblas_axpy(n2, 1, M + i * n2, 1, C21 + i * ldc, 1);
        cblas_axpy(n - n2, -1, M + i * n2, 1, C22 + i * ldc, 1);
    }

    // M3 = A11 * (B12 - B22)
    // rh = B12 - B22 is k2-by-(n - n2)
    rhlen = (n - n2) * sizeof(real);
    for (i = 0; i < k - k2; ++i)
    {
        memcpy(rh + i * (n - n2), B12 + i * ldb, rhlen);
        cblas_axpy(n - n2, -1, B22 + i * ldb, 1, rh + i * (n - n2), 1);
    }
    if (k2 > k - k2)
        memcpy(rh + (k2 - 1) * (n - n2), B12 + (k2 - 1) * ldb, rhlen);
    // Recursive call
    memset(M, 0, m2 * n2 * sizeof(real));
    StrassenATBRec(m2, n - n2, k2, alpha, A11, lda, rh, n - n2, M, n - n2, sM, slh, srh);
    // print_mat(M, m - m2, n2);
    // Update C12 and C22
    for (i = 0; i < m - m2; ++i)
    {
        cblas_axpy(n - n2, 1, M + i * (n - n2), 1, C12 + i * ldc, 1);
        cblas_axpy(n - n2, 1, M + i * (n - n2), 1, C22 + i * ldc, 1);
    }
    if (m2 > m - m2)
        cblas_axpy(n - n2, 1, M + (m2 - 1) * (n - n2), 1, C12 + (m2 - 1) * ldc, 1);

    // M4 = A22^T * (B21 - B11)
    // rh = B21 - B11 is (k - k2)-by-n2, since multiplication by A22 excludes the last row
    rhlen = n2 * sizeof(real);
    for (i = 0; i < k - k2; ++i)
    {
        memcpy(rh + i * n2, B21 + i * ldb, rhlen);
        cblas_axpy(n2, -1, B11 + i * ldb, 1, rh + i * n2, 1);
    }
    // Recursive call
    memset(M, 0, m2 * n2 * sizeof(real));
    StrassenATBRec(m - m2, n2, k - k2, alpha, A22, lda, rh, n2, M, n2, sM, slh, srh);
    // print_mat(M, m - m2, n2);
    // Update C11 and C21
    for (i = 0; i < m - m2; ++i)
    {
        cblas_axpy(n2, 1, M + i * n2, 1, C11 + i * ldc, 1);
        cblas_axpy(n2, 1, M + i * n2, 1, C21 + i * ldc, 1);
    }

    // M5 = (A11 + A21)^T * B22
    // lh = A11 + A21 is (k - k2)-by-m2, since multiplication by B22 excludes the last row
    lhlen = m2 * sizeof(real);
    for (i = 0; i < k - k2; ++i)
    {
        memcpy(lh + i * m2, A11 + i * lda, lhlen);
        cblas_axpy(m2, 1, A21 + i * lda, 1, lh + i * m2, 1);
    }
    // Recursive call
    memset(M, 0, m2 * n2 * sizeof(real));
    StrassenATBRec(m2, n - n2, k - k2, alpha, lh, m2, B22, ldb, M, n - n2, sM, slh, srh);
    // print_mat(M, m2, n - n2);
    // Update C11(-) and C12
    for (i = 0; i < m2; ++i)
    {
        cblas_axpy(n - n2, -1, M + i * (n - n2), 1, C11 + i * ldc, 1);
        cblas_axpy(n - n2, 1, M + i * (n - n2), 1, C12 + i * ldc, 1);
    }

    // M6 = (A12 - A11)^T * (B11 + B12)
    // lh = A12 - A11 is k2-by-m2
    // rh = B11 + B12 is k2-by-n2
    lhlen = (m - m2) * sizeof(real);
    rhlen = (n - n2) * sizeof(real);
    for (i = 0; i < k2; ++i)
    {
        memcpy(lh + i * m2, A12 + i * lda, lhlen);
        if (m2 > m - m2)
            *(lh + i * m2 + m - m2) = 0;
        cblas_axpy(m2, -1, A11 + i * lda, 1, lh + i * m2, 1);

        memcpy(rh + i * (n - n2), B11 + i * ldb, rhlen);
        cblas_axpy(n - n2, 1, B12 + i * ldb, 1, rh + i * (n - n2), 1);
    }
    // Recursive call
    memset(M, 0, m2 * n2 * sizeof(real));
    StrassenATBRec(m - m2, n - n2, k2, alpha, lh, m2, rh, n - n2, M, n - n2, sM, slh, srh);
    // print_mat(M, m - m2, n - n2);
    // Update C22
    for (i = 0; i < m - m2; ++i)
        cblas_axpy(n - n2, 1, M + i * (n - n2), 1, C22 + i * ldc, 1);

    // M7 = (A21 - A22) * (B21 + B22)
    // lh = A21 - A22 is (k - k2)-by-m2
    // rh = B21 + B22 is (k - k2)-by-n2
    lhlen = m2 * sizeof(real);
    rhlen = n2 * sizeof(real);
    for (i = 0; i < k - k2; ++i)
    {
        memcpy(lh + i * m2, A21 + i * lda, lhlen);
        cblas_axpy(m - m2, -1, A22 + i * lda, 1, lh + i * m2, 1);

        memcpy(rh + i * n2, B21 + i * ldb, rhlen);
        cblas_axpy(n - n2, 1, B22 + i * ldb, 1, rh + i * n2, 1);
    }
    // Recursive call
    memset(M, 0, m2 * n2 * sizeof(real));
    StrassenATBRec(m2, n2, k - k2, alpha, lh, m2, rh, n2, M, n2, sM, slh, srh);
    // print_mat(M, m - m2, n2);
    // Update C11
    for (i = 0; i < m2; ++i)
        cblas_axpy(n2, 1, M + i * n2, 1, C11 + i * ldc, 1);
}



void StrassenATB(const MKL_INT m, const MKL_INT n, const MKL_INT k,
                 const real alpha, const real* A, const MKL_INT lda,
                 const real* B, const MKL_INT ldb,
                 real* C, const MKL_INT ldc)
{
    if (k <= ATB_BASE || n <= ATB_BASE || m <= ATB_BASE)
    {
        cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                   m, n, k, alpha, A, lda, B, ldb, 1, C, ldc);
        return;
    }

    real* M = (real*)malloc((m * n) / 2 * sizeof(real));
    real* lh = (real*)malloc((m * k) / 2 * sizeof(real));
    real* rh = (real*)malloc((k * n) / 2 * sizeof(real));
    if (M == NULL || lh == NULL || rh == NULL)
    {
        fprintf(stderr, "Cannot allocate the support matrices.\n");
        exit(1);
    }

    StrassenATBRec(m, n, k, alpha, A, lda, B, ldb, C, ldc, M, lh, rh);

    free(M);
    free(lh);
    free(rh);
}


void StrassenATBC0(const MKL_INT m, const MKL_INT n, const MKL_INT k,
                   const real alpha, const real* A, const MKL_INT lda,
                   const real* B, const MKL_INT ldb,
                   real* C, const MKL_INT ldc)
{
    register MKL_INT i;
    for (i = 0; i < m; ++i)
        memset(C + i * ldc, 0, n * sizeof(real));

    StrassenATB(m, n, k, alpha, A, lda, B, ldb, C, ldc);
}