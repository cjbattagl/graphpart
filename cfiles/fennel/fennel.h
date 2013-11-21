static void usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args);

// CSR Conversion
void csr_to_metis (int n, int nnz, int *rowptr, int *colidx, idx_t **xadj, idx_t **adjncy, idx_t **vwgt, idx_t **adjwgt);
struct csr_matrix_t* mat_to_csr (struct coo_matrix_t* A);
void sort_coo (void* coord_array, const int length, enum value_type_t value_type);
static int coo_to_csr_convert(struct sparse_matrix_t* A);

// Calculations
static int run_fennel (const struct csr_matrix_t* A, int nparts, float gamma);
static int fennel_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx, 
    bool **parts, float alpha, float gamma, int *emptyverts);
static float calc_dc(float alpha, float gamma, int len);
