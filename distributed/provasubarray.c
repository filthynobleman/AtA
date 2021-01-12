#include <mpi.h>
#include <stdio.h>
 
int main(int argc, char *argv[])
{
    int myrank;
    MPI_Status status;
    MPI_Datatype subarray;
    int array[9] = { -1, 1, 2, 3, -2, -3, -4, -5, -6 };
    int array_size[] = {9};
    int array_subsize[] = {3};
    int array_start[] = {1};
    int i;
 
    MPI_Init(&argc, &argv);
 
    /* Create a subarray datatype */
    MPI_Type_create_subarray(1, array_size, array_subsize, array_start, MPI_ORDER_C, MPI_INT, &subarray);
    MPI_Type_commit(&subarray);
 
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 
    if (myrank == 0)
    {
        MPI_Send(array, 1, subarray, 1, 123, MPI_COMM_WORLD);
    }
    else if (myrank == 1)
    {
        for (i=0; i<9; i++)
            array[i] = 0;
        MPI_Recv(array, 1, subarray, 0, 123, MPI_COMM_WORLD, &status);
        for (i=0; i<9; i++)
            printf("array[%d] = %d\n", i, array[i]);
        fflush(stdout);
    }
 
    MPI_Finalize();
    return 0;
}
