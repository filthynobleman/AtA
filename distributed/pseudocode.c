void AtA(Matrix* A, Matrix* C)
{
    TaskTree Tree = SimulateExecution(A, C, NumProcs, alpha);


    size_t Level;
    TaskNode* Node = NULL;
    int IsInnerNode = 0;
    real* C0 = NULL;
    // Start from the leaf task and proceed up
    for (Level = Tree.Height; Level > 0; --Level)
    {
        // Retrieve my task at current level
        GetTaskNode(&Tree, MyID, Level, &Node);
        // If task is null, I did all my jobs and from now on I don't have anything to do
        if (Node == NULL)
            continue;

        // If task is not null, perform my task
        Task T = Node->T;
        // Do either A^T * A or A^T * B, depending on the task
        // WARNING: After the leaf task, all the subsequent tasks are only about recombining
        //          resuls, hence I only have to pass the result above
        if (! IsInnerNode)
        {
            // code...
            // code...
            // code...
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
        if (FatherID != MyID)
        {
            // If I have been generated by another process, send everything to him
            // Use MPI to send result to FatherID
            // code...
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
            real* CRes = (real*)calloc(CN * CM, sizeof(real));
            real* CP = NULL;
            int i;
            for (i = 0; i < vCommLen; ++i)
            {
                int P = vComm[i];
                Task ChildTask = Father->Children[i].T;
                if (P != MyID)
                {
                    CP = (real*)malloc(ChildTask.C.NumRows * ChildTask.T.NumCols * sizeof(real));
                    // Use MPI to receive from process P the result CP
                    // code...
                }
                else
                    CP = C0;
                // Combine the result from CP in matrix CRes
                // code... 
            }
        }

        // After the first iteration, I am no more in a leaf task
        IsInnerNode = 1;
    }
}