/**
 * @file       ata.h
 *
 * @brief      Interface for ATA algorithm.
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#ifndef ATA_ATA_H_
#define ATA_ATA_H_


#include "common.h"
#include "matrices.h"



void AtA(const MKL_INT N, const MKL_INT K, 
         const real alpha, const real* A, const MKL_INT lda,
         real* C, const MKL_INT ldc);
void AtA_MT(const Matrix* A, Matrix *C, integer NumThreads);





#endif // ATA_ATA_H_