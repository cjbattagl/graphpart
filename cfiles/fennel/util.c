#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/read_mm.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/interface.h>
#include <bebop/smc/sparse_matrix_ops.h>

#include <bebop/util/merge_sort.h>
#include <bebop/util/config.h>
#include <bebop/util/get_options.h>
#include <bebop/util/init.h>
#include <bebop/util/log.h>
#include <bebop/util/malloc.h>
#include <bebop/util/timer.h>
#include <bebop/util/util.h>

// #include <metis.h>
#include <assert.h>

#include "util.h"

// void checkCall(int rval) {
//   switch (rval) {
//     case METIS_OK:
//       return;
//     case METIS_ERROR_INPUT:
//       fprintf(stdout,"metismex:metisError / metis input error");
//       break;
//     case METIS_ERROR_MEMORY:
//       fprintf(stdout,"metismex:metisError / metis memory error");
//       break;
//     default:
//       fprintf(stdout,"metismex:metisError / unknown metis error");
//       break;
//   }
// }

const char* format_to_string(enum sparse_matrix_storage_format_t type ) {
  switch(type) {
    case CSR:
      return("CSR\n");
    case COO:
      return("COO\n");
    case CSC:
      return("CSC\n");
    default:
      return(stdout, "ERROR\n");
  }
}

const char* value_type_to_string(enum value_type_t type) {
  switch(type) {
    case REAL:
      return("REAL ");
    case PATTERN:
      return("PATTERN ");
    case COMPLEX:
      return("COMPLEX ");
    default:
      return("ERROR ");
  }
}

const char* symmetry_type_to_string(enum symmetry_type_t type) {
  switch(type) {
    case SYMMETRIC:
      return("SYMMETRIC");
    case UNSYMMETRIC:
      return("UNSYMMETRIC\n");
    case SKEW_SYMMETRIC:
      return("SKEW_SYMMETRIC\n");
    case HERMITIAN:
      return("HERMITIAN\n");
    default:
      return("ERROR\n");
  }
}

int coord_elem_by_row_real (const void* a, const void* b) {
  struct coord_elem_t
  {
    int r;
    int c;
    double val;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->r < y->r)
    return -1;
  else if (x->r > y->r)
    return +1;

  /* else */
  return 0;
}

int coord_elem_by_row_pattern (const void* a, const void* b) {
  struct coord_elem_t
  {
    int r;
    int c;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->r < y->r)
    return -1;
  else if (x->r > y->r)
    return +1;

  /* else */
  return 0;
}


/**
* Compares two coordinate-array elements by their column indices, and
* returns values like strcmp does. Comparison function for sorting.
*/
int coord_elem_by_col_real (const void* a, const void* b) {
  struct coord_elem_t
  {
    int r;
    int c;
    double val;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->c < y->c)
    return -1;
  else if (x->c > y->c)
    return +1;

  /* else */
  return 0;
}

int coord_elem_by_col_pattern (const void* a, const void* b) {
  struct coord_elem_t
  {
    int r;
    int c;
  };

  const struct coord_elem_t* x = ((struct coord_elem_t*) a);
  const struct coord_elem_t* y = ((struct coord_elem_t*) b);

  if (x->c < y->c)
    return -1;
  else if (x->c > y->c)
    return +1;

  /* else */
  return 0;
}


struct csr_matrix_t* mat_to_csr (struct coo_matrix_t* A) {
  int m = A->m;
  int n = A->n;
  int nnz = A->nnz;
  enum symmetry_type_t symmetry_type = A->symmetry_type;
  enum symmetric_storage_location_t symmetric_storage_location = A->symmetric_storage_location;
  enum value_type_t value_type = A->value_type;
  int *rowptr = NULL;
  int *colidx = NULL;
  void *__values = NULL;
  void *__coord_array = NULL;
  int i, j, currow = 0;
  struct csr_matrix_t* B = bebop_calloc (1, sizeof (struct csr_matrix_t));
  int index_base = A->index_base;

  /* bebop_set_debug_level (2); */

  fprintf (stdout, "=== coo_to_csr ===\n");
  fprintf (stdout, "index_base=%d\n",index_base);
  fprintf (stdout, "\tm = %d, n = %d, nnz = %d, value_type = %d\n", m, n, nnz, value_type);

  if (A->value_type == REAL)
    __values = bebop_calloc (nnz, sizeof (double));
  else if (A->value_type == PATTERN)
    __values = NULL;
  else {
      fprintf (stdout, "*** coo_to_csr: input matrix in COO format has"
         " invalid value type %d ***\n", A->value_type);
      bebop_exit (EXIT_FAILURE);
  }

  rowptr = bebop_calloc ((m+1), sizeof (int));
  colidx = bebop_calloc (nnz, sizeof (int));

  if (nnz == 0) {
      /* calloc fills in rowptr with zeros, which is all
we need if the matrix is empty */
      fprintf (stdout, "\tMatrix is empty, done\n");
      init_csr_matrix (B, m, n, nnz, __values, colidx, rowptr, UNSYMMETRIC,
                 UPPER_TRIANGLE, value_type, LIBRARY_DEALLOCATES,
                 &free, NO_COPY);
      fprintf (stdout, "=== Done with coo_to_csr_matrix ===\n");
      return B;
    }

  /* Intermediate conversion to coordinate array, so we can sort the entries */
  fprintf (stdout, "\tIntermediate conversion to coordinate array format\n");
  coo_matrix_to_coord_elem_array (&__coord_array, &nnz, A);

  fprintf (stdout, "\tSorting elements of coordinate array\n");
  sort_coo (__coord_array, nnz, value_type);

  if (value_type == REAL) {
      struct coord_array_t { int r; int c; double val; };
      struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;
      double* values = (double*) __values;

      fprintf (stdout, "\tReal value specialization\n");

      /*
      * Having sorted the elements of the coordinate array, first by initial
      * row, then within each row by initial column, now the first row with
      * nonzeros is coord_array[0].r.
      */
      fprintf (stdout, "\tInitializing start of rowptr\n");
      currow = coord_array[0].r - index_base;
      for (i = 0; i <= currow; i++)
        {
         /*
         * Until we get to first row with an entry in it, all the rowptrs
         * before then are zero. The rowptr for that first column is also
         * zero.
         */
         rowptr[i] = 0;
        }
      /* For each coord_array entry i, add it to the CSR matrix */
      fprintf (stdout, "\tAdding entries to CSR matrix\n");
      for (i = 0; i < nnz; i++)
        {
         //fprintf (stdout, "\t\ti = %d of %d\n", i, nnz);
         if (coord_array[i].r - index_base > currow)
         {
         /*
         * We may jump more than one row at a time, so set the rowptr
         * entries for the empty rows in between.
         */
         for (j = currow+1; j <= coord_array[i].r - index_base; j++)
                {
                 if (j - index_base < 0 || j - index_base > m)
                 {
                 bebop_log (0, "*** At entry %d of %d, "
                                "j = %d is out of the valid range "
                                "[%d,%d] ***\n",
                                i, nnz, j - index_base, 0, m);
                 bebop_free (values);
                 bebop_free (rowptr);
                 bebop_free (colidx);
                 bebop_free (B);
                 return NULL;
                 }
                 rowptr[j] = i;
                }

         currow = coord_array[i].r - index_base;
         }

         values[i] = coord_array[i].val;
         colidx[i] = coord_array[i].c - index_base;
        }

      /* Set the last entries in rowptr appropriately */
      for (j = currow+1; j <= m; j++)
        rowptr[j] = nnz;

      init_csr_matrix (B, m, n, nnz, __values, colidx, rowptr, symmetry_type,
                 symmetric_storage_location, value_type,
                 LIBRARY_DEALLOCATES, &free, NO_COPY);
    }
  else if (value_type == PATTERN)
    {
      struct coord_array_t { int r; int c; };
      struct coord_array_t* coord_array = (struct coord_array_t*) __coord_array;

      /*
* Having sorted the elements of the coordinate array, first by initial
* row, then within each row by initial column, now the first row with nonzeros is
* coord_array[0].r.
*/
      currow = coord_array[0].r - index_base;
      for (i = 0; i <= currow; i++) {
         /*
         * Until we get to first row with an entry in it, all the rowptrs
         * before then are zero. The rowptr for that first row is also zero.
         */
         rowptr[i] = 0;
         //fprintf (stdout, "\tset rowptr[%d] to %d\n", i, 0);
         }
      /* For each coord_array entry i, add it to the CSR matrix */
      for (i = 0; i < nnz; i++) {
      //fprintf (stdout, "\ti = %d, c.r[i] = %d, c.c[i] = %d\n", i, coord_array[i].r, coord_array[i].c);
         if (coord_array[i].r - index_base > currow) {
         /*
         * We may jump more than one row at a time, so set the rowptr
         * entries for the empty rows in between.
         */
         for (j = currow+1; j <= coord_array[i].r - index_base; j++) {
                 rowptr[j] = i;
                 //fprintf (stdout, "\tset rowptr[%d] to %d\n", j, i);
                 }
         currow = coord_array[i].r - index_base;
         }
                //fprintf (stdout, "\tset colidx[%d] to %d\n", i, coord_array[i].c - index_base);
        colidx[i] = coord_array[i].c - index_base;
         }

      /* Set the last entries in rowptr appropriately */
      for (j = currow+1; j <= m; j++) { rowptr[j] = nnz;
      //fprintf (stdout, "\tset rowptr[%d] to %d\n", j, nnz);
}

      init_csr_matrix (B, m, n, nnz, NULL, colidx, rowptr, symmetry_type,
                 symmetric_storage_location, value_type,
                 LIBRARY_DEALLOCATES, &free, NO_COPY);
    }

  bebop_free (__coord_array);
  fprintf (stdout, "=== Done with coo_to_csr ===\n");
  return B;
}

void sort_coo (void* coord_array, const int length, enum value_type_t value_type) {
  fprintf (stdout, "=== sort_coord_elem_array_for_csr_conversion ===\n");

  if (value_type == REAL) {
    struct coord_elem_t { int r; int c; double val; };
    merge_sort (coord_array, length, sizeof (struct coord_elem_t), coord_elem_by_col_real);
    merge_sort (coord_array, length, sizeof (struct coord_elem_t), coord_elem_by_row_real);
  }
  else if (value_type == PATTERN) {
    struct coord_elem_t { int r; int c; };
    merge_sort (coord_array, length, sizeof (struct coord_elem_t), coord_elem_by_col_pattern);
    merge_sort (coord_array, length, sizeof (struct coord_elem_t), coord_elem_by_row_pattern);
  }
  fprintf (stdout, "=== Done with sort_coord_elem_array_for_csr_conversion ===\n");
}

int coo_to_csr_convert(struct sparse_matrix_t* A) {
  assert(A->format == COO);
  struct csr_matrix_t* B = mat_to_csr (A->repr);
  if (B == NULL) { return -1; }
  destroy_coo_matrix (A->repr);
  A->repr = B; A->format = CSR;
  return 0;
}
