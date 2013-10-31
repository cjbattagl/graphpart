#ifndef INC_MATVEC_H
#define INC_MATVEC_H

#include "csr.h"

void csr_matvec__sequential (double* y, const csr_t* A, const double* x);
int matvec (const csr_t* Adata, const double* b, double* x, int n, double* rhist, int maxiter);
#endif