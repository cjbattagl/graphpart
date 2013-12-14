static void mpi_usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args);

// Calculations
int run_mpi_fennel (const struct csr_matrix_t* A, int nparts, float gamma);
int mpi_fennel_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx, 
    bool **parts, float alpha, float gamma, int *emptyverts); 
    
static int compute_cut(int *emptyparts, int *redparts, int *rowptr, int *colidx, bool **parts, int nparts, int n, FILE* out);
static float calc_dc(float alpha, float gamma, int len);
