/**
 * @file       queue.h
 *
 * @brief      Interface for a generic, simple queue structure.
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#ifndef ATA_QUEUE_H_
#define ATA_QUEUE_H_


#include "common.h"


struct Queue
{
	size_t ElementSize;
	size_t First;
	size_t Last;
	size_t QSize;
	size_t QCapacity;
	void*  Content;
	size_t QCapStep;
};
typedef struct Queue Queue;
typedef Queue queue_t;



void QueueInit(Queue* q, size_t ElementSize, size_t QCapacity, size_t QCapStep);
void QueueEnqueue(Queue *q, void* Element);
void QueueDequeue(Queue* q, void* Element);
void QueueFree(Queue *q);





#endif // ATA_QUEUE_H_