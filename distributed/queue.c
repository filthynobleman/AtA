/**
 * @file       queue.h
 *
 * @brief      Implementation of a generic, simple queue structure.
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#include "queue.h"


void QueueInit(Queue* q, size_t ElementSize, size_t QCapacity, size_t QCapStep)
{
	if (q == NULL)
	{
		fprintf(stderr, "AtA::QueueInit::NullPtr: Given queue is NULL.\n");
		exit(1);
	}

	q->ElementSize = ElementSize;
	q->QCapacity = QCapacity;
	q->QCapStep = QCapStep;
	q->First = 0;
	q->Last = 0;
	q->QSize = 0;

	if (q->QCapacity == 0)
		q->Content == NULL;
	else
	{
		q->Content = malloc(q->QCapacity * q->ElementSize);
		if (q->Content == NULL)
		{
			fprintf(stderr, "AtA::QueueInit::NullPtr: Cannot allocate the queue.\n");
			exit(1);
		}
	}
}


void QueueEnqueue(Queue *q, void* Element)
{
	if (q == NULL)
	{
		fprintf(stderr, "AtA::QueueEnqueue::NullPtr: Given queue is NULL.\n");
		exit(1);
	}

	if (Element == NULL)
	{
		fprintf(stderr, "AtA::QueueEnqueue::NullPtr: Given element is NULL.\n");
		exit(1);
	}

	if (q->QSize == q->QCapacity)
	{
		q->Content = realloc(q->Content, (q->QCapacity + q->QCapStep) * q->ElementSize);
		if (q->Content == NULL)
		{
			fprintf(stderr, "AtA::QueueEnqueue::NullPtr: Cannot reallocate queue content.\n");
			exit(1);
		}

		if (q->Last <= q->First)
		{
			if (q->Last > 0)
			{
				void* NewContent = malloc((q->QCapacity + q->QCapStep) * q->ElementSize);
				memcpy(NewContent, q->Content + q->First * q->ElementSize, (q->QCapacity - q->First) * q->ElementSize);
				memcpy(NewContent + (q->QCapacity - q->First) * q->ElementSize, q->Content, q->Last * q->ElementSize);
				free(q->Content);
				q->Content = NewContent;
				q->First = 0;
				q->Last = q->QCapacity;
				//memcpy(q->Content, q->Content + q->QCapacity * q->ElementSize, q->Last * q->ElementSize);
			}
			else
				q->Last += q->QCapacity;
		}
		q->QCapacity += q->QCapStep;
	}

	memcpy(q->Content + q->Last * q->ElementSize, Element, q->ElementSize);
	q->Last = (q->Last + 1) % q->QCapacity;
	q->QSize += 1;
}


void QueueDequeue(Queue* q, void* Element)
{
	if (q == NULL)
	{
		fprintf(stderr, "AtA::QueueDequeue::NullPtr: Given queue is NULL.\n");
		exit(1);
	}

	if (Element == NULL)
	{
		fprintf(stderr, "AtA::QueueDequeue::NullPtr: Given element is NULL.\n");
		exit(1);
	}

	memcpy(Element, q->Content + q->First * q->ElementSize, q->ElementSize);
	q->First = (q->First + 1) % q->QCapacity;
	q->QSize -= 1;
}


void QueueFree(Queue *q)
{
	if (q == NULL)
	{
		fprintf(stderr, "AtA::QueueFree::NullPtr: Given queue is NULL.\n");
		exit(1);
	}

	if (q->Content == NULL)
	{
		fprintf(stderr, "AtA::QueueFree::NullPtr: Queue content is NULL.\n");
		exit(1);
	}

	free(q->Content);
}