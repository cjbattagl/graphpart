static void usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args);
static int do_validate_matrix (int argc, char *argv[], struct arginfo* arglist);

//util
struct csr_matrix_t* mat_to_csr (struct coo_matrix_t* A);
void sort_coo (void* coord_array, const int length, enum value_type_t value_type);
int coord_elem_by_row_real (const void* a, const void* b);
int coord_elem_by_row_pattern (const void* a, const void* b);
int coord_elem_by_col_real (const void* a, const void* b);
int coord_elem_by_col_pattern (const void* a, const void* b);
static int coo_to_csr_convert(struct sparse_matrix_t* A);
float calc_dc(float alpha, float gamma, int len);

static int fennel_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx, 
    bool **parts, float alpha, float gamma, int *emptyverts);
static int run_fennel (const struct csr_matrix_t* A, int nparts, float gamma);
