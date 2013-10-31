#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/bcoo_matrix.h>
#include <bebop/smc/bcsr_matrix.h>
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/jad_matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <bebop/util/config.h>

#include "fennel.h"
#include "timer.h"
#include "matvec.h"
#include "csr.h"

struct stopwatch_t* g_timer = 0;

int main (int argc, char *argv[])
{
  /* Process arguments */
  if (argc != 2) {
    fprintf (stderr, "usage: %s <input-matrix>\n", argv[0]);
    return -1;
  }

  const char* matrix_filename = argv[1];

  /* Open and read matrix */
  csr_t A;
  read_matrix_market_real_sparse(matrix_filename, &A);

  stopwatch_init ();
  g_timer = stopwatch_create ();

  /* Run CG */
  fprintf (stdout, "\n===== Running matvec =====\n");
  run_driver (&A, "rhist__sequential.out");

  /* Done -- clean-up */
  stopwatch_destroy (g_timer);
  csr_free (&A);
  return 0;
}

static void run_driver (const csr_t* A, const char* rhist_filename)
{
  const int maxiter = getenv ("MAXITER") ? atoi (getenv ("MAXITER")) : A->m;
  const int numtrials = getenv ("NUMTRIALS") ? atoi (getenv ("NUMTRIALS")) : 40;

  double* b = NULL; /* right-hand side */
  double* x = NULL; /* solution vector */
  setup_mv (A->m, &b, &x);

  testing_driver (A, maxiter, b, x, rhist_filename);
  timing_driver (A, maxiter, numtrials, b, x);

  /* clean-up */
  free (x);
  free (b);
}

static void setup_mv (int n, double** p_b, double** p_x)
{
  if (n <= 0) return;

  if (p_b) {
    if (!(*p_b)) { // not yet allocated
      *p_b = (double *)malloc (n * sizeof (double));
      assert (*p_b);
    }
    for (int i = 0; i < n; ++i)
      (*p_b)[i] = 1.0;
  }
  if (p_x) {
    if (!(*p_x)) { // not yet allocated
      *p_x = (double *)malloc (n * sizeof (double));
      assert (*p_x);
    }
  }
}


static void timing_driver (const csr_t* A, const int maxiter, const int numtrials, double* b, double* x)
{
  if (numtrials <= 0) return;

  long double* times = (long double *)malloc (numtrials * sizeof (long double));
  assert (times);

  int num_iters = 0;
  g_malloc = 0;
  g_init = 0;
  g_scan = 0;


  for (int trial = 0; trial < numtrials; ++trial) {
    setup_mv (A->m, &b, &x);
    stopwatch_start (g_timer);
    num_iters += matvec(A, b, x, A->m, NULL, maxiter);
    times[trial] = stopwatch_stop (g_timer);
  }

  fprintf (stderr, "Raw execution times:");
  for (int trial = 0; trial < numtrials; ++trial)
    fprintf (stderr, " %Lg", times[trial]);
  fprintf (stderr, "\n\n");

  /* Summary statistics */
  long double t_min = times[0], t_max = times[0], t_sum = times[0];
  for (int trial = 1; trial < numtrials; ++trial) {
    if (times[trial] < t_min) t_min = times[trial];
    if (times[trial] > t_max) t_max = times[trial];
    t_sum += times[trial];
  }
  long double t_med = median__long_double (numtrials, times);

  const long double avg_iters = (long double)num_iters / numtrials;
  const long double num_flops = /* total # of floating-point operations */
    (long double)4.0 * A->m /* CG start-up */
    + avg_iters * (2.0 * A->nnz /* SpMV */
                   + 5.0 * A->m /* 2 dots + 3 axpys */)
    ;
  printf ("Minimum time (best performance): %Lg secs (~ %.2Lg Gflop/s)\n",
          t_min, 1.0e-9 * num_flops / t_min);
  printf ("Maximum time (worst performance): %Lg secs (~ %.2Lg Gflop/s)\n",
          t_max, 1.0e-9 * num_flops / t_max);
  printf ("Average time (performance): %Lg (~ %.2Lg Gflop/s)\n",
          t_sum / numtrials, 1.0e-9 * num_flops / t_sum * numtrials);
  printf ("Average malloc time: %Lg\n", g_malloc / numtrials);
  printf ("Average scan time: %Lg\n", g_scan / numtrials);
  printf ("Average init time: %Lg\n", g_init / numtrials);
  printf ("Median time (performance): %Lg secs (~ %.2Lg Gflop/s)\n",
          t_med, 1.0e-9 * num_flops / t_min);
}

static void testing_driver (const csr_t* A, const int maxiter, double* b, double* x, const char* rhist_filename)
{
  double* rhist = NULL;
  if (maxiter && rhist_filename) {
    rhist = (double *)malloc (maxiter * sizeof (double));
    assert (rhist);
  }

  assert ((b && x) || (!A->m));
  setup_mv (A->m, &b, &x);
  int retval = matvec(A, b, x, A->m, rhist, maxiter);

  int num_iters; /* Number of iterations executed */
  if (retval < 0) {
    fprintf (stderr, "Iteration failed to converge! (maxiter = %d)\n", maxiter);
    num_iters = maxiter + 1;
  } else {
    fprintf (stderr, "Converged after %d iterations.\n",
             retval);
    num_iters = retval;
  }

  /* Dump residuals */
  if (rhist && rhist_filename) {
    FILE* rhist_fp = fopen (rhist_filename, "w");
    assert (rhist_fp);
    for (int i = 0; i < num_iters; ++i)
      fprintf (rhist_fp, "%g\n", rhist[i]);
    fclose (rhist_fp);
    fprintf (stderr, "(See '%s' for convergence history)\n", rhist_filename);
  }

  if (rhist) free (rhist);

  assert (num_iters != (maxiter+1)); // Abort if no convergence
}

static int cmp__long_double (const void* pa, const void* pb)
{
  assert (pa && pb);
  long double a = ((long double *)pa)[0];
  long double b = ((long double *)pb)[0];
  if (a < b) return -1;
  if (a > b) return 1;
  return 0;
}

static long double median__long_double (int n, const long double* X)
{
  /* Makes a sorted copy */
  long double* X_copy = (long double *)malloc (n * sizeof (long double));
  assert (X_copy);
  memcpy (X_copy, X, n * sizeof (long double));
  qsort (X_copy, n, sizeof (long double), cmp__long_double);

  long double x_median;
  if ((n % 2) == 0) // even number of elements
    x_median = (long double)0.5 * (X_copy[n/2-1] + X_copy[n/2]);
  else // odd number of elements
    x_median = X_copy[n/2];

  free (X_copy);
  return x_median;
}
