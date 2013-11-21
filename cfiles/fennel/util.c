#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/sparse_matrix.h>
#include <metis.h>
#include "util.h"

void checkCall(int rval) {
  switch (rval) {
    case METIS_OK:
      return;
    case METIS_ERROR_INPUT:
      fprintf(stdout,"metismex:metisError / metis input error");
      break;
    case METIS_ERROR_MEMORY:
      fprintf(stdout,"metismex:metisError / metis memory error");
      break;
    default:
      fprintf(stdout,"metismex:metisError / unknown metis error");
      break;
  }
}

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