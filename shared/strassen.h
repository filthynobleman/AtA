/**
 * @file       strassen.h
 *
 * @brief      Interface for Strassen's algorithm performing A^T * B.
 *
 * @author     Filippo Maggioli, Viviana Arrigoni
 *             {maggioli,arrigoni}@di.uniroma1.it
 *             Sapienza, University of Rome - Department of Computer Science
 *             
 * @date       12 Jan 2021
 */
#ifndef ATA_STRASSEN_H_
#define ATA_STRASSEN_H_

#include "common.h"


void StrassenATB(const MKL_INT m, const MKL_INT n, const MKL_INT k,
                 const real alpha, const real* A, const MKL_INT lda,
                 const real* B, const MKL_INT ldb,
                 real* C, const MKL_INT ldc);

void StrassenATBC0(const MKL_INT m, const MKL_INT n, const MKL_INT k,
                   const real alpha, const real* A, const MKL_INT lda,
                   const real* B, const MKL_INT ldb,
                   real* C, const MKL_INT ldc);


// void StrassenATB_MT(const MKL_INT m, const MKL_INT n, const MKL_INT k,
//                     const real alpha, const real* A, const MKL_INT lda,
//                     const real* B, const MKL_INT ldb,
//                     real* C, const MKL_INT ldc,
//                     const MKL_INT NumThreads);

// void StrassenATBC0_MT(const MKL_INT m, const MKL_INT n, const MKL_INT k,
//                       const real alpha, const real* A, const MKL_INT lda,
//                       const real* B, const MKL_INT ldb,
//                       real* C, const MKL_INT ldc,
//                       const MKL_INT NumThreads);





#endif // ATA_STRASSEN_H_