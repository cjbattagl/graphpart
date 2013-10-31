#include "timer.h"
#include "csr.h"

static void run_driver (const csr_t* A, const char* rhist_filename);
static void setup_mv (int n, double** p_b, double** p_x);

static void timing_driver (const csr_t* A, const int maxiter, const int numtrials, double* b, double* x);
static void testing_driver (const csr_t* A, const int maxiter, double* b, double* x, const char* rhist_filename);

static int cmp__long_double (const void* pa, const void* pb);
static long double median__long_double (int n, const long double* X);
