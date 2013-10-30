#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>
#include <bebop/smc/csr_matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "csr.h"
#include "timer.h"

void csr_instantiate (int m, int nnz, int* ptr, int* ind, double* val, csr_t* A)
{
  assert (A);
  A->m = m;
  A->nnz = nnz;
  A->ptr = ptr;
  A->ind = ind;
  A->val = val;
}

void csr_free (csr_t* A)
{
  if (A) {
    if (A->ptr) { free (A->ptr); A->ptr = 0; }
    if (A->ind) { free (A->ind); A->ind = 0; }
    if (A->val) { free (A->val); A->val = 0; }
  }
}

void csr_readHB (const char* file, csr_t* B) {
  struct sparse_matrix_t* A = NULL;
  struct stopwatch_t* timer = NULL;
  int err;

  assert (file);

  timer = stopwatch_create ();
  assert (timer);

  fprintf (stderr, "Loading '%s'...", file); fflush (stderr);
  stopwatch_start (timer);
  A = load_sparse_matrix (HARWELL_BOEING, file); assert (A);
  stopwatch_stop (timer);
  fprintf (stderr, "done [%Lg secs].\n", stopwatch_elapsed (timer));

  fprintf (stderr, "Expanding to unsymmetric (if symmetric)..."); fflush (stderr);
  stopwatch_start (timer);
  err = sparse_matrix_expand_symmetric_storage (A); assert (!err);
  stopwatch_stop (timer);
  fprintf (stderr, "done [%Lg secs].\n", stopwatch_elapsed (timer));

  fprintf (stderr, "Converting to CSR..."); fflush (stderr);
  stopwatch_start (timer);
  sparse_matrix_convert (A, CSR);
  stopwatch_stop (timer);
  fprintf (stderr, "done [%Lg secs].\n", stopwatch_elapsed (timer));

  stopwatch_destroy (timer);

  assert (A->format == CSR);
  struct csr_matrix_t* A_csr = (struct csr_matrix_t *)A->repr;
  assert (A_csr);
  assert (A_csr->nnz == (A_csr->rowptr[A_csr->m] - A_csr->rowptr[0]));

  fprintf (stderr, "--> Matrix is %d x %d with %d non-zeros.\n", A_csr->m, A_csr->n, A_csr->nnz);
}
