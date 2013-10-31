#ifndef INC_MATVEC_H
#define INC_MATVEC_H

#include "csr.h"

void axpy (double* dest, double alpha, const double* x, const double* y, int n);
double dot (const double *x, const double *y, int n);
int matvec (const csr_t* Adata, const double* b, double* x, int n, double* rhist, int maxiter);
#endif