static void mpi_usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args);

// Calculations
int mpi_run_fennel(int* rowptr, int* colidx, int n, int n_local, 
    int nnz, int nnz_local, int v_offset, float gamma);

int mpi_fennel_kernel(int n, int n_local, int offset, int nparts, int *partsize, int *rowptr, int *colidx, 
    int **parts, float alpha, float gamma, int *emptyverts); 
    
static int mpi_compute_cut(int *emptyparts, int *redparts, int *rowptr, int *colidx, int **parts, 
    int nparts, int n_local, int offset, FILE* out);
    
static float calc_dc(float alpha, float gamma, int len);
