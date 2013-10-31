#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matvec.h"

int matvec (const csr_t* Adata, const double* b, double* x, int n, double* rhist, int maxiter)
{
  csr_matvec__sequential (b, Adata, x);
}

void
csr_matvec__sequential (double* y, const csr_t* A, const double* x)
{
  assert (A);
  assert ((x && y) || !A->m);

  for (int i = 0; i < A->m; ++i) {
    int k;
    double y_i = 0;
    for (int k = A->ptr[i]; k < A->ptr[i+1]; ++k)
      y_i += A->val[k] * x[A->ind[k]];
    y[i] = y_i;
  }
}
