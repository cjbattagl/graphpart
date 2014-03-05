// CSR Conversion
static void csr_to_metis (int n, int nnz, int *rowptr, int *colidx, idx_t **xadj, idx_t **adjncy, idx_t **vwgt, idx_t **adjwgt);

// Calculations
int run_fennel (const struct csr_matrix_t* A, int nparts, float gamma);
int fennel_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx, 
    bool **parts, float alpha, float gamma, int *emptyverts); 
int deg_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx, 
    bool **parts, float alpha, float gamma, int *emptyverts, int cutoff);
int sample_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx, 
    bool **parts, float alpha, float gamma, int *emptyverts, float prob);
    
static int compute_cut(int *emptyparts, int *redparts, int *rowptr, int *colidx, bool **parts, int nparts, int n, FILE* out, int cutoff);
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