static void usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args);
static int do_validate_matrix (int argc, char *argv[], struct arginfo* arglist);
//static void run_driver(const csr_t* A, const char* rhist_filename);
static void setup_mv (int n, double** p_b, double** p_x);

static int cmp__long_double (const void* pa, const void* pb);
static long double median__long_double (int n, const long double* X);
