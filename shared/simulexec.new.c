#include "simulexec.h"
#include "queue.h"


TaskTree SimulateExecution(Matrix* Matrix_A, Matrix* Matrix_C, size_t NumProcs, double alpha)
{
	Queue q;
	QueueInit(&q, sizeof(Task), NumProcs, NumProcs / 2);

	// Task* Tasks = (Task*)calloc(2 * NumProcs, sizeof(Task));
	TaskTree Tree;
	TaskTreeInit(&Tree);


	size_t LastUsed = 0;


	// Main task
	SubMatrix A = { Matrix_A, 0, 0, Matrix_A->NumRows, Matrix_A->NumCols };
	SubMatrix C = { Matrix_C, 0, 0, Matrix_A->NumCols, Matrix_A->NumCols };
	Task t = { ATA, A, A, C , NumProcs, 0 };
	size_t CurLevel;

	QueueEnqueue(&q, &t);
	TaskTreeAddTask(&Tree, t, 0);
	while (q.QSize > 0)
	{
		QueueDequeue(&q, &t);
		if (t.Budget == 1)
		{
			// if (Tasks[2 * t.ID].Type == NONE)
			// 	Tasks[2 * t.ID] = t;
			// else
			// 	Tasks[2 * t.ID + 1] = t;
			continue;
		}

		TaskTreeLastProcessLevel(&Tree, t.ID, &CurLevel);

		if (t.Type == ATA)
		{

			integer N2 = t.A.NumRows / 2;
			integer M2 = t.A.NumCols / 2;


			if (t.Budget == 2)
			{
				SubMatrix Ax2 = { t.A.M, t.A.BegRow, t.A.BegCol + M2, t.A.NumRows, t.A.NumCols - M2 };
				SubMatrix Ax1 = { t.A.M, t.A.BegRow, t.A.BegCol, t.A.NumRows, M2 };
				SubMatrix C11 = { t.C.M, t.C.BegRow, t.C.BegCol, M2, M2 };
				SubMatrix C2x = { t.C.M, t.C.BegRow + M2, t.C.BegCol, t.C.NumRows - M2, t.C.NumCols };
				// C2x = Ax2^T * A
				Task st1 = { ATB_SPECIAL, Ax2, t.A, C2x, 1, t.ID };
				QueueEnqueue(&q, &st1);
				TaskTreeAddTask(&Tree, st1, t.ID);
				// C11 = Ax1^T * Ax1
				Task st2 = { ATA, Ax1, Ax1, C11, 1, ++LastUsed };
				QueueEnqueue(&q, &st2);
				TaskTreeAddTaskLevel(&Tree, st2, t.ID, CurLevel);

				continue;
			}
			size_t AtBProcs = (size_t)(alpha * t.Budget);
			size_t AtAProcs = t.Budget - AtBProcs;
			if (AtBProcs == 0)
			{
				AtBProcs = 1;
				AtAProcs -= 1;
			}
			else if (AtAProcs < 2)
			{
				AtAProcs = 2;
				AtBProcs = t.Budget - AtAProcs;
			}
			// else if (AtAProcs == 0)
			// {
			// 	AtBProcs -= 1;
			// 	AtAProcs = 1;
			// }
			else if (AtAProcs % 2 == 1)
			{
				AtBProcs += 1;
				AtAProcs -= 1;
			}

			// AtB: C21 = A12^T * A11 + A22^T * A21 = Ax2^T * Ax1
			{
				SubMatrix Ax2 = { t.A.M, t.A.BegRow, t.A.BegCol + M2, t.A.NumRows, t.A.NumCols - M2 };
				SubMatrix Ax1 = { t.A.M, t.A.BegRow, t.A.BegCol, t.A.NumRows, M2 };
				SubMatrix C21 = { t.C.M, t.C.BegRow + M2, t.C.BegCol, t.A.NumCols - M2, M2 };
				Task st = { ATB, Ax2, Ax1, C21, AtBProcs, t.ID };
				QueueEnqueue(&q, &st);
				TaskTreeAddTask(&Tree, st, t.ID);
			}

			// AtA
			// C11 = A11^T * A11 + A21^T * A21 = Ax1^T * Ax1
			// C22 = A12^T * A12 + A22^T * A22 = Ax2^T * Ax2
			{
				SubMatrix Ax1 = { t.A.M, t.A.BegRow, t.A.BegCol, t.A.NumRows, M2 };
				SubMatrix C11 = { t.C.M, t.C.BegRow, t.C.BegCol, M2, M2 };
				Task st1 = { ATA, Ax1, Ax1, C11, AtAProcs / 2, ++LastUsed };
				QueueEnqueue(&q, &st1);
				TaskTreeAddTaskLevel(&Tree, st1, t.ID, CurLevel);


				SubMatrix Ax2 = { t.A.M, t.A.BegRow, t.A.BegCol + M2, t.A.NumRows, t.A.NumCols - M2 };
				SubMatrix C22 = { t.C.M, t.C.BegRow + M2, t.C.BegCol + M2, t.A.NumCols - M2, t.A.NumCols - M2 };
				Task st2 = { ATA, Ax2, Ax2, C22, AtAProcs / 2, ++LastUsed };
				QueueEnqueue(&q, &st2);
				TaskTreeAddTaskLevel(&Tree, st2, t.ID, CurLevel);
			}
		}
		else
		{
			integer NA2 = t.A.NumRows / 2;
			integer MA2 = t.A.NumCols / 2;
			integer NB2 = t.B.NumRows / 2;
			integer MB2 = t.B.NumCols / 2;

			// C = A^T * B
			// if budg is multiple of 4
			// 		C11 = A11^T * B11 + A21^T * B21
			// 		C21 = A12^T * B11 + A22^T * B21
			// 		C12 = A11^T * B12 + A21^T * B22
			// 		C22 = A12^T * B12 + A22^T * B22
			// else
			//		Ci = A^T * B_i
			if (t.Budget % 4 == 0)
			{
				// C11 = A11^T * B11 + A21^T * B21 = Ax1^T * Bx1
				SubMatrix Ax1 = { t.A.M, t.A.BegRow, t.A.BegCol, t.A.NumRows, MA2 };
				SubMatrix Bx1 = { t.B.M, t.B.BegRow, t.B.BegCol, t.B.NumRows, MB2 };
				SubMatrix C11 = { t.C.M, t.C.BegRow, t.C.BegCol, MA2, MB2 };
				Task st1 = { ATB, Ax1, Bx1, C11, t.Budget / 4, t.ID };
				QueueEnqueue(&q, &st1);
				TaskTreeAddTask(&Tree, st1, t.ID);
				// C21 = A12^T * B11 + A22^T * B21 = Ax2^T * Bx1
				SubMatrix Ax2 = { t.A.M, t.A.BegRow, t.A.BegCol + MA2, t.A.NumRows, t.A.NumCols - MA2 };
				SubMatrix C21 = { t.C.M, t.C.BegRow + MA2, t.C.BegCol, t.A.NumCols - MA2, MB2 };
				Task st2 = { ATB, Ax2, Bx1, C21, t.Budget / 4, ++LastUsed };
				QueueEnqueue(&q, &st2);
				TaskTreeAddTaskLevel(&Tree, st2, t.ID, CurLevel);
				// C12 = A11^T * B12 + A21^T * B22 = Ax1^T * Bx2
				SubMatrix Bx2 = { t.B.M, t.B.BegRow, t.B.BegCol + MB2, t.B.NumRows, t.B.NumCols - MB2 };
				SubMatrix C12 = { t.C.M, t.C.BegRow, t.C.BegCol + MB2, MA2, t.B.NumCols - MB2 };
				Task st3 = { ATB, Ax1, Bx2, C12, t.Budget / 4, ++LastUsed };
				QueueEnqueue(&q, &st3);
				TaskTreeAddTaskLevel(&Tree, st3, t.ID, CurLevel);
				// C22 = A12^T * B12 + A22^T * B22 = Ax2^T * Bx2
				SubMatrix C22 = { t.C.M, t.C.BegRow + MA2, t.C.BegCol + MB2, t.A.NumCols - MA2, t.B.NumCols - MB2 };
				Task st4 = { ATB, Ax2, Bx2, C22, t.Budget / 4, ++LastUsed };
				QueueEnqueue(&q, &st4);
				TaskTreeAddTaskLevel(&Tree, st4, t.ID, CurLevel);
			}
			else 
			{
				size_t MSB = t.B.NumCols / t.Budget;
				size_t i;
				for (i = 0; i < t.Budget; ++i)
				{
					size_t NumCols = (i == t.Budget - 1) ? MSB : (t.B.NumCols - MSB * (t.Budget - 1));
					SubMatrix Bxi = { t.B.M, t.B.BegRow, t.B.BegCol + i * MSB, t.B.NumRows, NumCols };
					SubMatrix Cxi = { t.C.M, t.C.BegRow, t.C.BegCol + i * MSB, t.C.NumRows, NumCols };
					Task st = { ATB, t.A, Bxi, Cxi, 1, i == 0 ? t.ID : ++LastUsed };
					QueueEnqueue(&q, &st);
					TaskTreeAddTaskLevel(&Tree, st, t.ID, CurLevel);
				}
			}
		}
	}

	QueueFree(&q);

	// return Tasks;
	return Tree;
}


TaskTree SimulateStrassen(Matrix* Matrix_A, Matrix* Matrix_B, Matrix* Matrix_C, size_t NumProcs)
{
	Queue q;
	QueueInit(&q, sizeof(Task), NumProcs, NumProcs / 2);

	// Task* Tasks = (Task*)calloc(2 * NumProcs, sizeof(Task));
	TaskTree Tree;
	TaskTreeInit(&Tree);


	size_t LastUsed = 0;


	// Main task
	SubMatrix A = { Matrix_A, 0, 0, Matrix_A->NumRows, Matrix_A->NumCols };
	SubMatrix B = { Matrix_B, 0, 0, Matrix_B->NumRows, Matrix_B->NumCols };
	SubMatrix C = { Matrix_C, 0, 0, Matrix_A->NumCols, Matrix_B->NumCols };
	Task t = { ATB, A, B, C , NumProcs, 0 };
	size_t CurLevel;

	QueueEnqueue(&q, &t);
	TaskTreeAddTask(&Tree, t, 0);
	while (q.QSize > 0)
	{
		QueueDequeue(&q, &t);
		if (t.Budget == 1)
		{
			// if (Tasks[2 * t.ID].Type == NONE)
			// 	Tasks[2 * t.ID] = t;
			// else
			// 	Tasks[2 * t.ID + 1] = t;
			continue;
		}

		TaskTreeLastProcessLevel(&Tree, t.ID, &CurLevel);

		
		integer NA2 = t.A.NumRows / 2;
		integer MA2 = t.A.NumCols / 2;
		integer NB2 = t.B.NumRows / 2;
		integer MB2 = t.B.NumCols / 2;

		// C = A^T * B
		// if budg is multiple of 4
		// 		C11 = A11^T * B11 + A21^T * B21
		// 		C21 = A12^T * B11 + A22^T * B21
		// 		C12 = A11^T * B12 + A21^T * B22
		// 		C22 = A12^T * B12 + A22^T * B22
		// else
		//		Ci = A^T * B_i
		if (t.Budget % 4 == 0)
		{
			// C11 = A11^T * B11 + A21^T * B21 = Ax1^T * Bx1
			SubMatrix Ax1 = { t.A.M, t.A.BegRow, t.A.BegCol, t.A.NumRows, MA2 };
			SubMatrix Bx1 = { t.B.M, t.B.BegRow, t.B.BegCol, t.B.NumRows, MB2 };
			SubMatrix C11 = { t.C.M, t.C.BegRow, t.C.BegCol, MA2, MB2 };
			Task st1 = { ATB, Ax1, Bx1, C11, t.Budget / 4, t.ID };
			QueueEnqueue(&q, &st1);
			TaskTreeAddTask(&Tree, st1, t.ID);
			// C21 = A12^T * B11 + A22^T * B21 = Ax2^T * Bx1
			SubMatrix Ax2 = { t.A.M, t.A.BegRow, t.A.BegCol + MA2, t.A.NumRows, t.A.NumCols - MA2 };
			SubMatrix C21 = { t.C.M, t.C.BegRow + MA2, t.C.BegCol, t.A.NumCols - MA2, MB2 };
			Task st2 = { ATB, Ax2, Bx1, C21, t.Budget / 4, ++LastUsed };
			QueueEnqueue(&q, &st2);
			TaskTreeAddTaskLevel(&Tree, st2, t.ID, CurLevel);
			// C12 = A11^T * B12 + A21^T * B22 = Ax1^T * Bx2
			SubMatrix Bx2 = { t.B.M, t.B.BegRow, t.B.BegCol + MB2, t.B.NumRows, t.B.NumCols - MB2 };
			SubMatrix C12 = { t.C.M, t.C.BegRow, t.C.BegCol + MB2, MA2, t.B.NumCols - MB2 };
			Task st3 = { ATB, Ax1, Bx2, C12, t.Budget / 4, ++LastUsed };
			QueueEnqueue(&q, &st3);
			TaskTreeAddTaskLevel(&Tree, st3, t.ID, CurLevel);
			// C22 = A12^T * B12 + A22^T * B22 = Ax2^T * Bx2
			SubMatrix C22 = { t.C.M, t.C.BegRow + MA2, t.C.BegCol + MB2, t.A.NumCols - MA2, t.B.NumCols - MB2 };
			Task st4 = { ATB, Ax2, Bx2, C22, t.Budget / 4, ++LastUsed };
			QueueEnqueue(&q, &st4);
			TaskTreeAddTaskLevel(&Tree, st4, t.ID, CurLevel);
		}
		else 
		{
			size_t MSB = t.B.NumCols / t.Budget;
			size_t i;
			for (i = 0; i < t.Budget; ++i)
			{
				size_t NumCols = (i == t.Budget - 1) ? MSB : (t.B.NumCols - MSB * (t.Budget - 1));
				SubMatrix Bxi = { t.B.M, t.B.BegRow, t.B.BegCol + i * MSB, t.B.NumRows, NumCols };
				SubMatrix Cxi = { t.C.M, t.C.BegRow, t.C.BegCol + i * MSB, t.C.NumRows, NumCols };
				Task st = { ATB, t.A, Bxi, Cxi, 1, i == 0 ? t.ID : ++LastUsed };
				QueueEnqueue(&q, &st);
				TaskTreeAddTaskLevel(&Tree, st, t.ID, CurLevel);
			}
		}
	}

	QueueFree(&q);

	// return Tasks;
	return Tree;
}


void MatlabCode(size_t N, size_t M, Task* Tasks, size_t NumProcs, FILE* stream)
{
	fprintf(stream, "A = rand(%lu, %lu);\n", N, M);
	fprintf(stream, "C = zeros(%lu, %lu);\n", M, M);
	fprintf(stream, "CBench = A\' * A;\n");

	size_t i;
	Task t;
	for (i = 0; i < NumProcs; ++i)
	{
		t = Tasks[2 * i];
		fprintf(stream, "%% Process %lu, Task 1\n", i);
		PrintMatHeader(&(t.C), 'C', 0, 0, stream);
		fprintf(stream, " = ");
		PrintMatHeader(&(t.C), 'C', 0, 0, stream);
		fprintf(stream, " + ");
		PrintMatHeader(&(t.A), 'A', 0, 0, stream);
		fprintf(stream, "\' * ");
		if (t.Type == ATA)
			PrintMatHeader(&(t.A), 'A', 0, 0, stream);
		else
			PrintMatHeader(&(t.B), 'A', 0, 0, stream);
		fprintf(stream, ";\n");


		t = Tasks[2 * i + 1];
		if (t.Type == NONE)
			continue;
		fprintf(stream, "%% Process %lu, Task 2\n", i);
		PrintMatHeader(&(t.C), 'C', 0, 0, stream);
		fprintf(stream, " = ");
		PrintMatHeader(&(t.C), 'C', 0, 0, stream);
		fprintf(stream, " + ");
		PrintMatHeader(&(t.A), 'A', 0, 0, stream);
		fprintf(stream, "\' * ");
		if (t.Type == ATA)
			PrintMatHeader(&(t.A), 'A', 0, 0, stream);
		else
			PrintMatHeader(&(t.B), 'A', 0, 0, stream);
		fprintf(stream, ";\n");
	}


	fprintf(stream, "disp(norm(tril(C) - tril(CBench), \'fro\'));\n");
}


void PrintMatHeader(SubMatrix* M, char MName, size_t IndentLevel, int NewLine, FILE* stream)
{
	int i;
	for (i = 0; i < IndentLevel; ++i)
		fprintf(stream, "    ");
	fprintf(stream, "%c(%d:%d, %d:%d)", MName, 
			M->BegRow + 1, M->BegRow + M->NumRows,
			M->BegCol + 1, M->BegCol + M->NumCols);
	if (NewLine)
		fprintf(stream, "\n");
}



void PrintTask(Task t, size_t IndentLevel)
{
	if (t.Type == NONE)
		return;
	size_t i;
	for (i = 0; i < IndentLevel; ++i)
		fprintf(stdout, "    ");
	printf("Task is %s.\n", t.Type == ATA ? "AtA" : "AtB");
	PrintMatHeader(&(t.C), 'C', IndentLevel + 1, 0, stdout);
	printf(" = ");
	PrintMatHeader(&(t.A), 'A', 0, 0, stdout);
	printf("^T * ");
	if (t.Type == ATA)
		PrintMatHeader(&(t.A), 'A', 0, 1, stdout);
	else
		PrintMatHeader(&(t.B), 'A', 0, 1, stdout);
}



void TaskNodeFree(TaskNode* node)
{
	if (node->Children == NULL)
		free(node);
	else
	{
		size_t i;
		TaskNodeFree(node->Children);
	}
}

void TaskTreeFree(TaskTree* tree)
{
	if (tree == NULL)
		return;
	TaskNodeFree(tree->Root);
}



void TaskTreeInit(TaskTree* Tree)
{
	if (Tree == NULL)
	{
		fprintf(stderr, "AtA::TaskTreeInit::NullPtr: Cannot initialize a NULL tree.\n");
		exit(1);
	}
	Tree->Height = 0;
	Tree->Root = NULL;
}

void TaskTreeAddTaskLevel(TaskTree* Tree, Task NewTask, size_t Father, size_t Level)
{
	if (Tree == NULL)
	{
		fprintf(stderr, "AtA::TaskTreeAddTaskLevel::NullPtr: Given tree is NULL.\n");
		exit(1);
	}

	if (Level > Tree->Height)
	{
		fprintf(stderr, "AtA::TaskTreeAddTaskLevel::IllegalValue: Level is higher than height (%lu > %lu).\n", 
				Level, Tree->Height);
		exit(1);
	}

	if (Tree->Root == NULL)
	{
		Tree->Root = (TaskNode*)malloc(sizeof(TaskNode));
		if (Tree->Root == NULL)
		{
			fprintf(stderr, "AtA::TaskTreeAddTaskLevel::NullPtr: Failed root allocation.\n");
			exit(1);
		}
		Tree->Root->T = NewTask;
		Tree->Root->Father = NULL;
		Tree->Root->Children = NULL;
		Tree->Root->NumChildren = 0;
		Tree->Root->Tree = Tree;
		Tree->Root->Level = 0;
		return;
	}

	TaskNode* Node;
	GetTaskNode(Tree, Father, Level, &Node);
	if (Node == NULL)
	{
		fprintf(stderr, "AtA::TaskTreeAddTaskLevel::IDNotFound: Cannot find any node with id %lu at level %lu.\n", 
				Father, Level);
		exit(1);
	}

	size_t Len = Node->NumChildren;
	Node->Children = (TaskNode*)realloc(Node->Children, (Len + 1) * sizeof(TaskNode));
	if (Node->Children == NULL)
	{
		fprintf(stderr, "AtA::TaskTreeAddTaskLevel::NullPtr: Cannot reallocate children.\n");
		exit(1);
	}
	TaskNode NewNode = { NewTask, Node, NULL, 0, Tree, Level + 1 };
	Node->Children[Len] = NewNode;
	Node->NumChildren += 1;
	if (Level == Tree->Height)
		Tree->Height += 1;
}

void TaskTreeAddTask(TaskTree* Tree, Task NewTask, size_t Father)
{
	if (Tree == NULL)
	{
		fprintf(stderr, "AtA::TaskTreeAddTask::NullPtr: Given tree is NULL.\n");
		exit(1);
	}

	if (Tree->Root == NULL)
	{
		Tree->Root = (TaskNode*)malloc(sizeof(TaskNode));
		if (Tree->Root == NULL)
		{
			fprintf(stderr, "AtA::TaskTreeAddTask::NullPtr: Failed root allocation.\n");
			exit(1);
		}
		Tree->Root->T = NewTask;
		Tree->Root->Father = NULL;
		Tree->Root->Children = NULL;
		Tree->Root->NumChildren = 0;
		Tree->Root->Tree = Tree;
		Tree->Root->Level = 0;
		return;
	}

	size_t i;
	TaskNode* Node;
	for (i = Tree->Height; i >= 0; --i)
	{
		// TaskTreeAddTaskLevel(Tree, NewTask, Father, i);
		GetTaskNode(Tree, Father, i, &Node);
		if (Node == NULL)
			continue;

		size_t Len = Node->NumChildren;
		Node->Children = (TaskNode*)realloc(Node->Children, (Len + 1) * sizeof(TaskNode));
		if (Node->Children == NULL)
		{
			fprintf(stderr, "AtA::TaskTreeAddTask::NullPtr: Cannot reallocate children.\n");
			exit(1);
		}
		TaskNode NewNode = { NewTask, Node, NULL, 0, Tree, i + 1 };
		Node->Children[Len] = NewNode;
		Node->NumChildren += 1;
		if (i == Tree->Height)
			Tree->Height += 1;
		return;
	}
	fprintf(stderr, "AtA::TaskTreeAddTask::IDNotFound: Cannot find any node with id %lu.\n", Father);
	exit(1);
}


void GetChildrenCommunicator(TaskTree* Tree, size_t ID, size_t Level, int** vComm, int* vCommLen)
{
	if (Tree == NULL)
	{
		fprintf(stderr, "AtA::GetChildrenCommunicator::NullPtr: Given tree is NULL.\n");
		exit(1);
	}
	if (vComm == NULL)
	{
		fprintf(stderr, "AtA::GetChildrenCommunicator::NullPtr: Given communicator is NULL.\n");
		exit(1);
	}
	if (Tree == NULL)
	{
		fprintf(stderr, "AtA::GetChildrenCommunicator::NullPtr: Given communicator length is NULL.\n");
		exit(1);
	}


	TaskNode* Node;
	GetTaskNode(Tree, ID, Level, &Node);
	if (Node->NumChildren == 0)
	{
		*vCommLen = 0;
		*vComm = NULL;
		return;
	}

	*vCommLen = Node->NumChildren;
	*vComm = (int*)malloc((*vCommLen) * sizeof(int));
	size_t i;
	for (i = 0; i < *vCommLen; ++i)
		(*vComm)[i] = Node->Children[i].T.ID;
}


void GetTask(TaskTree* Tree, size_t ID, size_t Level, Task** T)
{
	if (Tree == NULL)
	{
		fprintf(stderr, "AtA::GetTask::NullPtr: Given tree is NULL.\n");
		exit(1);
	}
	if (T == NULL)
	{
		fprintf(stderr, "AtA::GetTask::NullPtr: Given task is NULL.\n");
		exit(1);
	}

	TaskNode* Node;
	GetTaskNode(Tree, ID, Level, &Node);
	if (Node == NULL)
	{
		fprintf(stderr, "AtA::GetTask::IDNotFound: Cannot find any node with id %lu.\n", ID);
		exit(1);
	}

	*T = &(Node->T);
}



void GetTaskNode(TaskTree* Tree, size_t ID, size_t Level, TaskNode** Node)
{
	if (Tree == NULL)
	{
		fprintf(stderr, "AtA::GetTaskNode::NullPtr: Given tree is NULL.\n");
		exit(1);
	}
	if (Node == NULL)
	{
		fprintf(stderr, "AtA::GetTaskNode::NullPtr: Given node is NULL.\n");
		exit(1);
	}

	Queue q;
	QueueInit(&q, sizeof(TaskNode*), 100, 100);
	QueueEnqueue(&q, &(Tree->Root));

	TaskNode* Node_;
	while (q.QSize > 0)
	{
		QueueDequeue(&q, &Node_);
		if (Node_->T.ID == ID && Node_->Level == Level)
		{
			*Node = Node_;
			QueueFree(&q);
			return;
		}

		size_t i;
		for (i = 0; i < Node_->NumChildren; ++i)
		{
			TaskNode* Child = Node_->Children + i;
			QueueEnqueue(&q, &Child);
		}
	}
	QueueFree(&q);
	*Node = NULL;
}


void TaskTreeLastProcessLevel(TaskTree* Tree, size_t ID, size_t* Level)
{
	if (Level == NULL)
	{
		fprintf(stderr, "AtA::TaskTreeLastProcessLevel::NullPtr: Given level is NULL.\n");
		exit(1);
	}
	*Level = Tree->Height;
	TaskNode* Node;
	for (; (*Level) >= 0; (*Level) -= 1)
	{
		GetTaskNode(Tree, ID, *Level, &Node);
		if (Node == NULL)
			continue;
		return;
	}
	*Level = Tree->Height + 1;
}


void LaTeXCodeNode(TaskNode* Node, FILE* stream)
{
	Task T = Node->T;
	size_t ID = T.ID;
	TaskType Type = T.Type;
	char* TypeName = Type == ATA ? "$A^\\top A$" : "$A^\\top B$";
	size_t ia = T.A.BegRow;
	size_t ja = T.A.BegCol;
	size_t na = T.A.NumRows;
	size_t ma = T.A.NumCols;

	size_t ib = Type == ATA ? T.A.BegRow : T.B.BegRow;
	size_t jb = Type == ATA ? T.A.BegCol : T.B.BegCol;
	size_t nb = Type == ATA ? T.A.NumRows : T.B.NumRows;
	size_t mb = Type == ATA ? T.A.NumCols : T.B.NumCols;

	size_t ic = T.C.BegRow;
	size_t jc = T.C.BegCol;
	size_t nc = T.C.NumRows;
	size_t mc = T.C.NumCols;
	fprintf(stream, ".{%lu: %s ($C_{%lu:%lu, %lu:%lu} = A_{%lu:%lu, %lu:%lu}^\\top A_{%lu:%lu, %lu:%lu}$)}\n", 
			ID, TypeName, 
			ic, ic + nc, jc, jc + mc,
			ia, ia + na, ja, ja + ma, 
			ib, ib + nb, jb, jb + mb);
	// fprintf(stream, ".%lu\n", Node->T.ID);
	// 
	// printf("Node %lu has %lu children.\n", Node->T.ID, Node->NumChildren);

	if (Node->Children == NULL)
		return;

	size_t i;
	for (i = 0; i < Node->NumChildren; ++i)
	{
		fprintf(stream, "[\n");
		LaTeXCodeNode(Node->Children + i, stream);
		fprintf(stream, "]\n");
	}
}



void LaTeXCodeTree(TaskTree* Tree, size_t NumProcs, FILE* stream)
{
	fprintf(stream, "\\documentclass{article}\n");
	fprintf(stream, "\\usepackage{amsfonts}\n");
	fprintf(stream, "\\usepackage{tikz}\n");
	fprintf(stream, "\\usepackage{tikz-qtree}\n");
	fprintf(stream, "\\usepackage[left=0cm,right=0cm,top=0cm,bottom=0cm,paperwidth=%lucm,paperheight=%lucm]{geometry}\n", 
			(Tree->Height + 1) * 15, 4 * NumProcs / 5 );
	fprintf(stream, "\n");
	fprintf(stream, "\\pagenumbering{gobble}\n");
	fprintf(stream, "\n");
	fprintf(stream, "\\begin{document}\n");
	fprintf(stream, "\\begin{figure}\n");
	fprintf(stream, "\\centering\n");
	fprintf(stream, "\\begin{tikzpicture}[grow\'=right,level distance=5in]\n");

	fprintf(stream, "\\Tree[\n");
	LaTeXCodeNode(Tree->Root, stream);
	fprintf(stream, "]\n");

	fprintf(stream, "\\end{tikzpicture}\n");
	fprintf(stream, "\\caption{A tree of %lu processes distributing a matrix $A \\in \\mathbb{R}^{%d\\times%d}$. The height of the tree is %lu}\n",
			NumProcs, Tree->Root->T.A.NumRows, Tree->Root->T.A.NumCols, Tree->Height);
	fprintf(stream, "\\end{figure}\n");
	fprintf(stream, "\\end{document}\n");
}