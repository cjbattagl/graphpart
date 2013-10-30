#if !defined (INC_CSR_H)
#define INC_CSR_H /*!< csr.h included. */

long double g_malloc;
long double g_init;
long double g_scan;

typedef struct
{
  int m;       /*! no. of rows */
  int nnz;     /*!< no. of stored non-zeros */
  int* ptr;    /*! row pointers, 0-based indices */
  int* ind;    /*! column indices, 0-based */
  double* val; /*! stored values */
} csr_t;

/**
 *  \brief Instantiates a CSR matrix using the data given.
 *  \note This routine makes a __shallow-copy__ of the inputs.
 */
void csr_instantiate (int m, int nnz, int* ptr, int* ind, double* val,
                      csr_t* A);

/**
 *  \brief Reads the matrix from the specified file, which must be in
 *  Harwell-Boeing format, and returns it as a handle to a CSR
 *  matrix. If the matrix is stored symmetrically, it is expanded to
 *  full format, i.e., the non-zeros in both the lower and upper
 *  triangles are represented explicitly.
 */
void csr_readHB (const char* filename, csr_t* B);
void csr_matvec__sequential (double* y, const csr_t* A, const double* x);
void csr_free (csr_t* A);

#endif
