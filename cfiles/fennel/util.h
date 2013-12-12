void checkCall(int rval);

const char* format_to_string(enum sparse_matrix_storage_format_t type);
const char* value_type_to_string(enum value_type_t type);
const char* symmetry_type_to_string(enum symmetry_type_t type);

int coord_elem_by_row_real (const void* a, const void* b);
int coord_elem_by_row_pattern (const void* a, const void* b);
int coord_elem_by_col_real (const void* a, const void* b);
int coord_elem_by_col_pattern (const void* a, const void* b);

struct csr_matrix_t* mat_to_csr (struct coo_matrix_t* A);
void sort_coo (void* coord_array, const int length, enum value_type_t value_type);
int coo_to_csr_convert(struct sparse_matrix_t* A);
