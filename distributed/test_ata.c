/**
 * @file       test_ata.c
 *
 * @brief      Test distributed memory ATA-D algorithm
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */

#include "ATA.h"
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>
#include <mpi.h>

void Usage(const char* ProgName)
{
    printf("Usage: %s [-N n] [-K k] [-t]\n", ProgName);
    printf("\n");
    printf("Generates a random matrix A (k-by-n) and computes the\n");
    printf("product C = A^T * A (n-by-n).\n");
    printf("If the optional argument -t is given, the execution time\n");
    printf("is printed on stdout. Combined with -c, also gives the time\n");
    printf("for cblas_?syrk\n");
    printf("\n");
    printf("    -N default is 1000.\n");
    printf("    -K default is 1000.\n");
    printf("    -t default is false.\n");
    printf("    -c default is false.\n");
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
    int GetTime = 0;
    int Check = 0;

    int opt;
    while ((opt = getopt(argc, argv, "ctN:K:")) != -1)
    {
        switch(opt)
        {
            case 'c': Check = 1; break;
            case 't': GetTime = 1; break;
            case 'N': N = atol(optarg); break;
            case 'K': K = atol(optarg); break;
            default: Usage(argv[0]); exit(EXIT_FAILURE);
        }
    }

    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    struct timeval Beg, End;
    double Time;

    real* A;
    real* CCheck;
    if (rank == 0)
    {
        printf("Allocating matrices... ");
        fflush(stdout);
        TimerStart(Beg);
        A = (real*)malloc(K * N * sizeof(real));
        if (A == NULL)
        {
            fprintf(stderr, "Failed to allocate A\n");
            exit(EXIT_FAILURE);
        }
        CCheck = NULL;
        if (Check)
        {
            CCheck = (real*)calloc(N * N, sizeof(real));
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

        printf("Randomly filling matrix A.\n");
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
    }
    else
    {
        A = NULL;
        CCheck = NULL;
    }
    Matrix mA = { A, K, N };
    Matrix mC = { NULL, N, N };


    if (rank == 0)
    {
        printf("Computing C = A^T * A with ATA-D algorithm.\n");
        fflush(stdout);
    }
    Time = AtAn(&mA, &mC);
    if (rank == 0)
    {
        if (GetTime)
            printf("Elapsed time is: %.3f seconds.\n", Time);
        else
            printf("Done.\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (!Check)
    {
        MPI_Finalize();
        return EXIT_SUCCESS;
    }


    if (rank == 0)
    {
        printf("Computing C = A^T * A with cblas_?syrk routine.\n");
        fflush(stdout);
        TimerStart(Beg);
        cblas_syrk(CblasRowMajor, CblasLower, CblasTrans, N, K, 1, A, N, 0, CCheck, N);
        TimerStop(End);
        Time = TimerGet(Beg, End);
        if (GetTime)
            printf("Elapsed time is: %.3f seconds.\n", Time);
        else
            printf("Done.\n");


        real Error = 0;
        real Norm = 0;

        real* C = mC.A;

        register size_t i, j;
        register size_t cnt = 0;
        for (i = 0; i < N; ++i)
        {
            for (j = 0; j <= i; ++j)
            {
                real Dij = CCheck[i * N + j] - C[cnt++];
                Error += Dij * Dij;
                Norm += CCheck[i * N + j] * CCheck[i * N + j];
            }
        }

        printf("Absolute error: %.3e\n", sqrt(Error));
        printf("Relative error: %.3e\n", sqrt(Error) / sqrt(Norm));
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();


    return 0;
}