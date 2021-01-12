#include "common.h"
#include "simulexec.h"


int main(const int const argc, char const **argv)
{
    size_t N = 100;
    size_t M = 100;
    size_t NumProcs = 1;
    double alpha = 0.5;
    if (argc == 2)
        NumProcs = atol(argv[1]);
    else if (argc == 3)
    {
        N = atol(argv[1]);
        M = atol(argv[2]);
    }
    else if (argc == 4)
    {
        N = atol(argv[1]);
        M = atol(argv[2]);
        NumProcs = atol(argv[3]);
    }
    else if (argc == 5)
    {
        N = atol(argv[1]);
        M = atol(argv[2]);
        NumProcs = atol(argv[3]);
        alpha = atof(argv[4]);
    }


    Matrix A = { NULL, N, M };
    Matrix C = { NULL, M, M };
    TaskTree Tree = SimulateExecution(&A, &C, NumProcs, alpha);

    FILE* stream = fopen("tree.tex", "w");
    LaTeXCodeTree(&Tree, NumProcs, stream);
    fclose(stream);


    return 0;
}