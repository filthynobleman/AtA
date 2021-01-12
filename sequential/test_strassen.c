/**
 * @file       test_strassen.c
 *
 * @brief      Test Strassen's algorithm.
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#include "strassen.h"
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>


void Usage(const char* ProgName)
{
    printf("Usage: %s [-M m] [-N n] [-K k] [-t] [-c]\n", ProgName);
    printf("\n");
    printf("Generates two random matrices A (k-by-m) and B (k-by-n),\n");
    printf("and computes the product C = A^T * B (m-by-k).\n");
    printf("If the optional argument -c is given, also computes the\n");
    printf("product with cblas_?gemm.\n");
    printf("If the optional argument -t is given, the execution time\n");
    printf("is printed on stdout. Combined with -c, also gives the time\n");
    printf("for cblas_?gemm\n");
    printf("\n");
    printf("    -M default is 1000.\n");
    printf("    -N default is 1000.\n");
    printf("    -K default is 1000.\n");
    printf("    -t default is false.\n");
    printf("    -c default is false.\n");
}


int main(int argc, char * const *argv)
{
    printf("Configuration option: ");
#ifdef DOUBLE_PRECISION
    printf("Using double precision for floating point values.\n");
#else
    printf("Using single precision for floating point values.\n");
#endif

    MKL_INT M = 1000;
    MKL_INT N = 1000;
    MKL_INT K = 1000;
    int GetTime = 0;
    int Check = 0;

    int opt;
    while ((opt = getopt(argc, argv, "ctN:M:K:")) != -1)
    {
        switch(opt)
        {
            case 'c': Check = 1; break;
            case 't': GetTime = 1; break;
            case 'M': M = atol(optarg); break;
            case 'N': N = atol(optarg); break;
            case 'K': K = atol(optarg); break;
            default: Usage(argv[0]); exit(EXIT_FAILURE);
        }
    }


    struct timeval Beg, End;
    double Time;

    printf("Allocating matrices... ");
    fflush(stdout);
    TimerStart(Beg);
    real* A = (real*)malloc(K * M * sizeof(real));
    if (A == NULL)
    {
        fprintf(stderr, "Failed to allocate A\n");
        exit(EXIT_FAILURE);
    }
    real* B = (real*)malloc(K * N * sizeof(real));
    if (B == NULL)
    {
        fprintf(stderr, "Failed to allocate B\n");
        exit(EXIT_FAILURE);
    }
    real* C = (real*)calloc(M * N, sizeof(real));
    if (C == NULL)
    {
        fprintf(stderr, "Failed to allocate C\n");
        exit(EXIT_FAILURE);
    }
    real* CCheck = NULL;
    if (Check)
    {
        CCheck = (real*)calloc(M * N, sizeof(real));
        if (CCheck == NULL)
        {
            fprintf(stderr, "Failed to allocate CCheck\n");
            exit(EXIT_FAILURE);
        }
    }

    TimerStop(End);
    Time = TimerGet(Beg, End);
    if (GetTime)
        printf("Elapsed time is: %.3f seconds.\n", Time);
    else
        printf("Done.\n");

    printf("Randomly filling matrices.\n");
    fflush(stdout);
    TimerStart(Beg);
    register MKL_INT i;
    for (i = 0; i < K * M; ++i)
        A[i] = (real)rand() / (real)RAND_MAX;
    for (i = 0; i < K * N; ++i)
        B[i] = (real)rand() / (real)RAND_MAX;
    TimerStop(End);
    Time = TimerGet(Beg, End);
    if (GetTime)
        printf("Elapsed time is: %.3f seconds.\n", Time);
    else
        printf("Done.\n");


    printf("Computing C = A^T * B with Strassen's algorithm.\n");
    fflush(stdout);
    TimerStart(Beg);
    StrassenATB(M, N, K, 1, A, M, B, N, C, N);
    TimerStop(End);
    Time = TimerGet(Beg, End);
    if (GetTime)
        printf("Elapsed time is: %.3f seconds.\n", Time);
    else
        printf("Done.\n");


    if (!Check)
        return EXIT_SUCCESS;


    printf("Computing C = A^T * B with cblas_?gemm routine.\n");
    fflush(stdout);
    TimerStart(Beg);
    cblas_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1, A, M, B, N, 0, CCheck, N);
    TimerStop(End);
    Time = TimerGet(Beg, End);
    if (GetTime)
        printf("Elapsed time is: %.3f seconds.\n", Time);
    else
        printf("Done.\n");


    real Error = 0;
    real Norm = 0;

    for (i = 0; i < M * N; ++i)
    {
        real Dij = CCheck[i] - C[i];
        Error += Dij * Dij;
        Norm += CCheck[i] * CCheck[i];
    }

    printf("Absolute error: %.3e\n", sqrt(Error));
    printf("Relative error: %.3e\n", sqrt(Error) / sqrt(Norm));





    return 0;
}