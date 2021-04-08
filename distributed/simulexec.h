/**
 * @file       simulexec.h
 *
 * @brief      Interface exposing functions and data structures for generating the
 *             ATA recursion tree.
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#ifndef ATA_SIMULEXEC_H_
#define ATA_SIMULEXEC_H_ 


#include "matrices.h"


enum TaskType
{
	NONE = 0,
	ATA = 1,
	ATB = 2,
    ATB_SPECIAL = 3
};
typedef enum TaskType TaskType;


struct Task
{
	TaskType  Type;
	SubMatrix A;
	SubMatrix B;
	SubMatrix C;
	size_t 	  Budget;
	size_t    ID;
};
typedef struct Task Task;


struct TaskNode;
typedef struct TaskNode TaskNode;

struct TaskTree
{
    size_t Height;
    TaskNode* Root;
};
typedef struct TaskTree TaskTree;

struct TaskNode
{
    Task T;
    TaskNode* Father;
    TaskNode* Children;
    size_t NumChildren;
    TaskTree* Tree;
    size_t Level;
};



void TaskTreeInit(TaskTree* Tree);
void TaskTreeAddTask(TaskTree* Tree, Task NewTask, size_t Father);
void TaskTreeAddTaskLevel(TaskTree* Tree, Task NewTask, size_t Father, size_t Level);
void GetChildrenCommunicator(TaskTree* Tree, size_t ID, size_t Level, int** vComm, int* vCommLen);
void GetTask(TaskTree* Tree, size_t ID, size_t Level, Task** T);
void GetTaskNode(TaskTree* Tree, size_t ID, size_t Level, TaskNode** Node);
void TaskTreeLastProcessLevel(TaskTree* Tree, size_t ID, size_t* Level);
void TaskTreeFree(TaskTree* Tree);




TaskTree SimulateExecution(Matrix* Matrix_A, Matrix* Matrix_C, size_t NumProcs, double alpha);

void MatlabCode(size_t N, size_t M, Task* Tasks, size_t NumProcs, FILE* stream);
void LaTeXCodeTree(TaskTree* Tree, size_t NumProcs, FILE* stream);
void PrintMatHeader(SubMatrix* M, char MName, size_t IndentLevel, int NewLine, FILE* stream);
void PrintTask(Task t, size_t IndentLevel);



#endif // ATA_SIMULEXEC_H_