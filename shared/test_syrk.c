/**
 * @file       test_syrk.c
 *
 * @brief      Test multi-threaded CBLAS routine ?syrk.
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#include "common.h"
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>


void Usage(const char* ProgName)
{
    printf("Usage: %s [-N n] [-K k] [-P p] [-t] [-c]\n", ProgName);
    printf("\n");
    printf("Generates a random matrix A (k-by-n) and computes the\n");
    printf("product C = A^T * A (n-by-n).\n");
    printf("The optional argument -P specifies the number of threads.\n");
    printf("If the optional argument -t is given, the execution time\n");
    printf("is printed on stdout.\n");
    printf("\n");
    printf("    -N default is 1000.\n");
    printf("    -K default is 1000.\n");
    printf("    -P default is 16.\n");
    printf("    -t default is false.\n");
}


int main(int argc, char* const *argv)
{
    printf("Configuration option: ");
#ifdef DOUBLE_PRECISION
    printf("Using double precision for floating point values.\n");
#else
    printf("Using single precision for floating point values.\n");
#endif

    MKL_INT N = 1000;
    MKL_INT K = 1000;
    MKL_INT NumThreads = 16;
    int GetTime = 0;

    int opt;
    while ((opt = getopt(argc, argv, "tN:K:P:")) != -1)
    {
        switch(opt)
        {
            case 'P': NumThreads = atol(optarg); break;
            case 't': GetTime = 1; break;
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
    real* A = (real*)malloc(K * N * sizeof(real));
    if (A == NULL)
    {
        fprintf(stderr, "Failed to allocate A\n");
        exit(EXIT_FAILURE);
    }
    real* C = (real*)calloc(N * N, sizeof(real));
    if (C == NULL)
    {
        fprintf(stderr, "Failed to allocate C\n");
        exit(EXIT_FAILURE);
    }
TimerStop(End);
    Time = TimerGet(Beg, End);
    if (GetTime)
        printf("Elapsed time is: %.3f seconds.\n", Time);
    else
        printf("Done.\n");

    printf("Randomly filling matrix A... ");
    fflush(stdout);
    TimerStart(Beg);
    register MKL_INT i;
    for (i = 0; i < K * N; ++i)
        A[i] = (real)rand() / (real)RAND_MAX;
    TimerStop(End);
    Time = TimerGet(Beg, End);
    if (GetTime)
        printf("Elapsed time is: %.3f seconds.\n", Time);
    else
        printf("Done.\n");


    printf("Computing C = A^T * A with cblas_?syrk routine... ");
    fflush(stdout);
    TimerStart(Beg);
#ifdef USE_MKL
    mkl_set_num_threads(NumThreads);
#else
    omp_set_num_threads(NumThreads);
#endif
    cblas_syrk(CblasRowMajor, CblasLower, CblasTrans, N, K, 1, A, N, 0, C, N);
    TimerStop(End);
    Time = TimerGet(Beg, End);
    if (GetTime)
        printf("Elapsed time is: %.3f seconds.\n", Time);
    else
        printf("Done.\n");





    return 0;
}