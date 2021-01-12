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
#include <assert.h>
#include "matrices.h"
#include "simulexec.h"


// The best choice for ALPHA seems to be 0.5 according to theory
#define ALPHA 0.5


void GenerateRandomMatrix(real* M, size_t n, size_t m)
{
	register size_t i;
	for (i = 0; i < n * m; i++)
		M[i] = rand() / ((real)RAND_MAX);
}

void AtAn(Matrix* A, Matrix* C){


	int size, id;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if (C == NULL)
	{
		fprintf(stderr, "Ata: C is NULL.");
		exit(1);
	}
	if (id == 0){
		C->A = (real *)calloc(A->NumCols * (A->NumCols + 1) / 2, sizeof(real));
		if (C->A == NULL)
		{
			fprintf(stderr, "AtA: Cannot allocate C->");
			exit(1);
		}
	}
	integer N = A->NumRows;
	integer M = A->NumCols;
	C->NumRows = M;
	C->NumCols = M;
    TaskTree Tree = SimulateExecution(A, C, size, ALPHA);


    size_t Level;
    TaskNode* Node = NULL;
    int IsInnerNode = 0;
    real* C0 = NULL;
    real* A0 = NULL;
    real* B0 = NULL;
    Matrix sC0;
    MPI_Status status; 
    // Start from the leaf task and proceed up
    for (Level = Tree.Height; Level > 0; --Level)
    {
        // Retrieve my task at current level
        GetTaskNode(&Tree, id, Level, &Node);
        // If task is null, I did all my jobs and from now on I don't have anything to do
        if (Node == NULL)
            continue;

        // If task is not null, perform my task
        Task T = Node->T;
        // Do either A^T * A or A^T * B, depending on the task
        // WARNING: After the leaf task, all the subsequent tasks are only about recombining
        //          resuls, hence I only have to pass the result above
        if (! IsInnerNode)
        {// I am a leaf I do task T
            double ExecStart = MPI_Wtime();
            if (T.Type == ATA){
                integer n = T.A.NumRows;
                integer m = T.A.NumCols;
                integer ia = T.A.BegRow;
                integer ja = T.A.BegCol;

                A0 = (real*)malloc(n * m * sizeof(real));
                
                for (int i = 0; i < n; ++i)
                    memcpy(A0 + i * m, A->A + (ia + i) * M + ja, m * sizeof(real));

                C0 = (real*)malloc(m * m * sizeof(real)); //QUI
                
                Matrix sA0 = {A0, n, m};
                sC0.A = C0;
                sC0.NumRows = m;
                sC0.NumCols = m; // va spostato fuori dall'if
                mat_syrk(&sA0, &sC0, 1, 0); //UNCOMMENTING THIS -> UNCOMMENTING C0 MALLOC QUI
                for (register int i = 1; i < sC0.NumRows; ++i)
                    memcpy(C0 + i * (i + 1) / 2, C0 + i * m, (i + 1) * sizeof(real));
                // AtA(&sA0, &sC0);
                // C0 = sC0.A;
                free(A0);
            } else {
                integer n = T.A.NumRows;
                integer m = T.A.NumCols;
                integer k = T.B.NumCols;

                A0 = (real*)malloc(n * m * sizeof(real));
                B0 = (real*)malloc(n * k * sizeof(real));
                
                for (int i = 0; i < n; ++i)
                    memcpy(A0 + i * m, A->A + SUBMATIDX(&(T.A), i, 0), m * sizeof(real));
                
                for (int i = 0; i < n; ++i)
                    memcpy(B0 + i * k, A->A + SUBMATIDX(&(T.B), i, 0), k * sizeof(real));
                
                C0 = (real*)malloc(m * k * sizeof(real)); // QUI
                
                Matrix sA0 = {A0, n, m};
                Matrix sB0 = {B0, n, k};
                sC0.A = C0;
                sC0.NumRows = m;
                sC0.NumCols = k;// va spostato fuori dall'if
                mat_gemtm(&sA0, &sB0, &sC0, 1, 0); //UNCOMMENTING THIS -> UNCOMMENTING C0 MALLOC  QUI
                //MatMult(&sA0, &sB0, &sC0);
                
                free(A0);
                free(B0);
            }
#if defined(SYNCHRONIZE_PROCESSES) || defined(ONLY_COMPUTE)
            MPI_Barrier(MPI_COMM_WORLD);
#ifdef ONLY_COMPUTE
            return;
#endif
#endif
        } else {
            sC0.A = C0;
            sC0.NumRows = T.C.NumRows;
            sC0.NumCols = T.C.NumCols;
        }

        // Now I could have been generated either from a father or from myself
        // In both cases, I need to generate a communicator
        TaskNode* Father = Node->Father;
        size_t FatherID = Father->T.ID;
        int vCommLen;
        int* vComm;
        GetChildrenCommunicator(&Tree, FatherID, Level - 1, &vComm, &vCommLen);

        // Create the MPI communicator
        // code...
        if (FatherID != id)
        {
            // If I have been generated by another process, send everything to him
            // Use MPI to send result to FatherID
            if (T.Type == ATA)
                MPI_Send(sC0.A, sC0.NumRows * (sC0.NumCols + 1) / 2, MPI_RREAL, FatherID, 0, MPI_COMM_WORLD);
            else
                MPI_Send(sC0.A, sC0.NumRows * sC0.NumCols, MPI_RREAL, FatherID, 0, MPI_COMM_WORLD);
            // After that, I don't need anymore C0
            free(C0);
            C0 = NULL;
        }
        else
        {
            // If I have been generated by myself, receive the results and combine them
            // WARNING: Because of the structure of the tree and because of how the
            //          communicator is created, the father process is always the first
            //          in the communicator
            // WARNING: Instruction on how to combine the results are not given by
            //          the current node Node, but by node Father and its list of
            //          children Father->Children, whose length is Father->NumChildren.
            integer CN = Father->T.C.NumRows;
            integer CM = Father->T.C.NumCols;
            real* CRes; 
            if (Father->T.Type == ATA)
                CRes = (real*)calloc(CN * (CM + 1) / 2, sizeof(real));
            else
                CRes = (real*)calloc(CN * CM, sizeof(real));
            size_t CPSize = 0;
            for (int i = 0; i < vCommLen; ++i)
            {
                Task ChildTask = Father->Children[i].T;
                size_t CPSizeLoc = ChildTask.C.NumRows * ChildTask.C.NumCols;
                CPSize = CPSize > CPSizeLoc ? CPSize : CPSizeLoc;
            }
            real* CP = (real*)malloc(CPSize * sizeof(real));
            for (int i = 0; i < vCommLen; ++i)
            {
                int P = vComm[i];
                Task ChildTask = Father->Children[i].T;
                if (P != id)
                {
                    // CP = (real*)malloc(ChildTask.C.NumRows * ChildTask.C.NumCols * sizeof(real));
                    // Use MPI to receive from process P the result CP
                    if (ChildTask.Type == ATA){
                        MPI_Recv(CP, ChildTask.C.NumRows * (ChildTask.C.NumCols + 1) / 2, MPI_RREAL, P, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    else
                        MPI_Recv(CP, ChildTask.C.NumRows * ChildTask.C.NumCols, MPI_RREAL, P, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                }
                else{
                    // CP = C0;
                    if (T.Type == ATA)
                        memcpy(CP, C0, T.C.NumCols * (T.C.NumRows + 1) / 2 * sizeof(real));
                    else
                        memcpy(CP, C0, T.C.NumCols * T.C.NumRows * sizeof(real));
                    free(C0);
                }
                SubMatrix CF = Father->T.C;
                // Combine the result from CP in matrix CRes
                double CombineStart = MPI_Wtime();
                if (ChildTask.Type == ATA)
                {
                    if (Father->T.Type == ATA)
                    {
                        size_t BegRowIdx;
                        for (register int j = 0; j < ChildTask.C.NumRows; ++j)
                        {
                            BegRowIdx = (ChildTask.C.BegRow - CF.BegRow + j) * (ChildTask.C.BegRow - CF.BegRow + j + 1) / 2;
                            cblas_axpy(j + 1, 1, CP + j * (j + 1) / 2,
                                       1, CRes + BegRowIdx + ChildTask.C.BegCol - CF.BegCol, 1);
                        }
                    }
                    else
                    {
                        for (register int j = 0; j < ChildTask.C.NumRows; ++j)
                            cblas_axpy(j + 1, 1, CP + j * (j + 1) / 2,
                                       1, CRes + (ChildTask.C.BegRow - CF.BegRow + j) * CM + ChildTask.C.BegCol - CF.BegCol, 1);
                    }
                }
                else
                {
                    if (Father->T.Type == ATA)
                    {
                        size_t BegRowIdx;
                        for (register int j = 0; j < ChildTask.C.NumRows; ++j)
                        {
                            BegRowIdx = (ChildTask.C.BegRow - CF.BegRow + j) * (ChildTask.C.BegRow - CF.BegRow + j + 1) / 2;
                            if (ChildTask.Type == ATB)
                                cblas_axpy(ChildTask.C.NumCols, 1, CP + j * ChildTask.C.NumCols,
                                       1, CRes + BegRowIdx + ChildTask.C.BegCol - CF.BegCol, 1);
                            else if (ChildTask.Type == ATB_SPECIAL)
                                cblas_axpy(ChildTask.C.BegRow - CF.BegRow + j + 1, 1, CP + j * ChildTask.C.NumCols,
                                       1, CRes + BegRowIdx + ChildTask.C.BegCol - CF.BegCol, 1);
                        }
                    }
                    else
                    {
                        for (register int j = 0; j < ChildTask.C.NumRows; ++j)
                            cblas_axpy(ChildTask.C.NumCols, 1, CP + j * ChildTask.C.NumCols,
                                       1, CRes + (ChildTask.C.BegRow - CF.BegRow + j) * CM + ChildTask.C.BegCol - CF.BegCol, 1);
                    }
                }
            }
            free(CP);
            C0 = CRes;
        }
        if (vCommLen > 0){
            free(vComm);
        }
        
        // After the first iteration, I am no more in a leaf task
        IsInnerNode = 1;
    }
    if (id == 0)
    	C->A = C0;
	return;
}



double RunBench(real* A, size_t N, size_t M)
{
	int id, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	double t1, t2, diff;

#ifdef DEBUG_MODE
	printf("id %d\n", id);
	printf("Generating random %lu-by-%lu matrix... ", N, M);
	t1 = MPI_Wtime();
#endif
	if (A == NULL)
	{
		fprintf(stderr, "Cannot allocate A->");
		exit(1);
	}
	if (id == 0){
		GenerateRandomMatrix(A, N, M);
	}
	MPI_Bcast(A, N * M, MPI_RREAL, 0, MPI_COMM_WORLD);


	Matrix A0 = {A, N, M};
#ifdef DEBUG_MODE
	t2 = MPI_Wtime();
	diff = t2 - t1;
	printf("%1.2f seconds.\n", diff);
#endif
	
	MPI_Barrier(MPI_COMM_WORLD);
	Matrix AtAC0;
	if (size > 1){
#ifdef DEBUG_MODE
		printf("Computing product with AtA->.. ");
#endif
        t1 = MPI_Wtime();
		AtAn(&A0, &AtAC0);
        t2 = MPI_Wtime();
		diff = t2 - t1;
#ifdef DEBUG_MODE
		printf("%1.2f seconds. id %d\n", diff, id);
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	} else {
#ifdef DEBUG_MODE
		printf("Computing product with AtA->.. ");
#endif
        t1 = MPI_Wtime();
        AtAC0.A = (real*)malloc(M * M * sizeof(real));
        AtAC0.NumRows = M;
        AtAC0.NumCols = M;
		mat_syrk(&A0, &AtAC0, 1, 1);
        t2 = MPI_Wtime();
        for (register int i = 1; i < AtAC0.NumRows; ++i)
            memcpy(AtAC0.A + i * (i + 1) / 2, AtAC0.A + i * M, (i + 1) * sizeof(real));
		diff = t2 - t1;
#ifdef DEBUG_MODE
		printf("%1.2f seconds.\n", diff);
#endif
	}

#ifdef DEBUG_MODE
	real* CblasC;
	
	if (id == 0){

		printf("Computing product with CBLAS... ");
		t1 = MPI_Wtime();
		CblasC = (real*)calloc(M * M, sizeof(real));
		if (CblasC == NULL)
		{
			fprintf(stderr, "Cannot allocate C->");
			free(A);
			exit(1);
		}
	    Matrix CblasC0 = {CblasC, M, M};
		mat_syrk(&A0, &CblasC0, 1., 0.);
		t2 = MPI_Wtime();
		diff = t2 - t1;
		printf("%1.2f seconds.\n", diff);
        FILE* deb = fopen("cblas.txt", "w");
        fprint_mat(CblasC0.A, M, M, M, M, deb);
        fclose(deb);

		printf("|| C - C2 || = %.3e\n", mat_err(&CblasC0, &AtAC0));
		free(CblasC);
		free(AtAC0.A);
	}
#else
    if (id == 0)
        free(AtAC0.A);
#endif


	return diff;
}


int main(int argc, const char const **argv)
{
    MPI_Init(&argc, &argv);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);  
    printf("numProcs %d\n", size);
    fflush(stdout); 
   if (argc < 5)
    {
        fprintf(stderr, "Not enough input arguments.\n");
        fprintf(stderr, "Usage: %s NumRuns NumRows NumCols Output\n", argv[0]);
        exit(1);
    }

    size_t NumRuns = atol(argv[1]);
    size_t NumRows = atol(argv[2]);
    size_t NumCols = atol(argv[3]);
    real* A = (real*)malloc(NumRows * NumCols * sizeof(real));
    if (A == NULL)
    {
        fprintf(stderr, "Cannot allocate A->");
        exit(1);
    }

    srand(0);

    double* Times = (double*)malloc(NumRuns * sizeof(double));

    size_t i;
    for (i = 0; i < NumRuns; ++i)
    {
        Times[i] = RunBench(A, NumRows, NumCols);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(A);

    int ID;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    if (ID == 0)
    {
        FILE *stream = fopen(argv[4], "w");
        for (i = 0; i < NumRuns; ++i)
            fprintf(stream, "%.6f\n", Times[i]);
        fclose(stream);
    }

    MPI_Finalize();



    
    return 0;
}