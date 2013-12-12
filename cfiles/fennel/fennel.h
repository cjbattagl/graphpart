static void usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args);
static struct csr_matrix_t* mat_to_csr (struct coo_matrix_t* A);
static void sort_coo (void* coord_array, const int length, enum value_type_t value_type);
static int coo_to_csr_convert(struct sparse_matrix_t* A);
