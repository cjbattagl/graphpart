static void usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args);
static int do_validate_matrix (int argc, char *argv[], struct arginfo* arglist);

static int run_fennel (const struct csr_matrix_t* A, int nparts);