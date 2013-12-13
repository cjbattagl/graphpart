// CSR Conversion
static void csr_to_metis (int n, int nnz, int *rowptr, int *colidx, idx_t **xadj, idx_t **adjncy, idx_t **vwgt, idx_t **adjwgt);
static struct csr_matrix_t* mat_to_csr (struct coo_matrix_t* A);
static void sort_coo (void* coord_array, const int length, enum value_type_t value_type);
static int coo_to_csr_convert(struct sparse_matrix_t* A);

// Calculations
int run_fennel (const struct csr_matrix_t* A, int nparts, float gamma);
int fennel_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx, 
    bool **parts, float alpha, float gamma, int *emptyverts); 
int sample_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx, 
    bool **parts, float alpha, float gamma, int *emptyverts, float prob);
    
static int compute_cut(int *emptyparts, int *redparts, int *rowptr, int *colidx, bool **parts, int nparts, int n, FILE* out);
static float calc_dc(float alpha, float gamma, int len);

struct parameter_data {
#define METISMEX_OPTION_WGTFLAG -1
    int wgtflag;
#define METISMEX_OPTION_ADJWGT -2
    int adjwgt;
#define METISMEX_OPTION_VSIZE -3
    idx_t *vsize;
};

#define FUNCNAMELEN 25