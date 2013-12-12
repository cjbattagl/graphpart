/**
* Driver program for our FENNEL implementation
* Note that portions of this code are adapted from the Sparse Matrix Converter,
* provided with the following LICENSE:
*
* Copyright (c) 2008, Regents of the University of California
* All rights reserved.
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the
* following conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in
* the documentation and/or other materials provided with the
* distribution.
*
* * Neither the name of the University of California, Berkeley, nor
* the names of its contributors may be used to endorse or promote
* products derived from this software without specific prior
* written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
* STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
* OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/read_mm.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/interface.h>
#include <bebop/smc/sparse_matrix_ops.h>

#include <bebop/util/merge_sort.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <bebop/util/config.h>
#include <bebop/util/get_options.h>
#include <bebop/util/init.h>
#include <bebop/util/log.h>
#include <bebop/util/malloc.h>
#include <bebop/util/timer.h>
#include <bebop/util/util.h>

#include <metis.h>

#include "fennel.h"
#include "fennel_kernel.h"

#include "randperm.h"
#include "util.h"

struct
cmdlineopts_t {
  char* input_filename;
  char* input_file_format;
  char* output_filename;
  char* output_file_format;
} opts;

static void
usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args) {
  fprintf (out, "Usage:\n");
  fprintf (out, "fennel [options] <in-filename> <in-format>\n");
  fprintf (out, "<in-filename>: name of file containing the sparse matrix to read in\n");
  fprintf (out, "<in-format>: format of the input file (\"HB\", \"ML\", \"MM\")\n");
  fprintf (out, "<out-filename>: name of file to which to output results\n");
  fprintf (out, "[options]: -e -- expand symmetric into unsymmetric storage\n");
  fprintf (out, " -v -- verbose mode\n");
  fprintf (out, "EX: ./fennel -v -e 'test.mtx' 'MM' 4\n");
}

int main (int argc, char *argv[]) {
  extern int optind;
  enum sparse_matrix_file_format_t informat = 0;
  enum sparse_matrix_file_format_t outformat = 0;
  struct sparse_matrix_t* A = NULL;
  struct arginfo *arglist = NULL;
  double seconds = 0.0;
  int errcode = 0;

  bebop_default_initialize (argc, argv, &errcode);
  if (errcode != 0) {
      fprintf (stderr, "*** Failed to initialize BeBOP Utility Library "
         "(error code %d) ***\n", errcode);
      bebop_exit (EXIT_FAILURE);
  }

  const char* matrix_filename = argv[1];
  register_usage_function (usage);
  arglist = register_arginfo (arglist, 'v', NULLARG, NULL, "If specified, ac"
                         "tivate verbose mode");
                         

  arglist = register_arginfo (arglist, 'e', NULLARG, NULL, "If specified, ex"
                         "pand the input matrix from symmetric storage "
                         "into unsymmetric storage");
  arglist = register_arginfo (arglist, 'a', NULLARG, NULL, "If specified, va"
                         "lidate the input matrix, without outputting a"
                         "nything");

  get_options (argc, argv, arglist, NULL);

  if (argc - optind != 3) {
      fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
         "e specified, but there should be %d ***\n", argc - optind, 3);
      dump_usage (stderr, argv[0], arglist, NULL);
      bebop_exit (EXIT_FAILURE); /* stops logging */
  }

  opts.input_filename = argv[optind];
  opts.input_file_format = argv[optind+1];
  int parts = atoi(argv[optind+2]);
  //opts.output_filename = argv[optind+2];
  
  if (strcmp (opts.input_file_format, "HB") == 0 ||
      strcmp (opts.input_file_format, "hb") == 0) { informat = HARWELL_BOEING; }
  else if (strcmp (opts.input_file_format, "MM") == 0 ||
         strcmp (opts.input_file_format, "mm") == 0) { informat = MATRIX_MARKET; }
  else if (strcmp (opts.input_file_format, "ML") == 0 ||
           strcmp (opts.input_file_format, "ml") == 0) { informat = MATLAB; }
  else {
      fprintf (stderr, "*** Unsupported input file format \"%s\" ***\n",
         opts.input_file_format);
      dump_usage (stderr, argv[0], arglist, NULL);
      destroy_arginfo_list (arglist);
      bebop_exit (EXIT_FAILURE); /* stops logging */
  }

  if (got_arg_p (arglist, 'v')) {
      printf ("Loading sparse matrix...");
      fflush (stdout); /* Otherwise the message may not be printed until
                         after loading is complete */
  }
  
  seconds = get_seconds();
  A = load_sparse_matrix (informat, opts.input_filename);
  seconds = get_seconds() - seconds;
  
  if (A == NULL) {
      fprintf (stderr, "*** Failed to load input matrix file \"%s\" ***\n",
         opts.input_filename);
      destroy_arginfo_list (arglist);
      bebop_exit (1); /* stops logging */
  }
  
  if (got_arg_p (arglist, 'v')) { printf ("done, in %g seconds\n", seconds); }
  if (got_arg_p (arglist, 'e')) {
    if (got_arg_p (arglist, 'v')) {
      printf ("Expanding sparse matrix into unsymmetric storage...");
      fflush (stdout);
    }
    seconds = get_seconds();
    errcode = sparse_matrix_expand_symmetric_storage (A);
    seconds = get_seconds() - seconds;
    if (errcode != 0) {
      fprintf (stderr, "*** Failed to expand matrix into symmetric storage ***\n");
      destroy_sparse_matrix (A);
      destroy_arginfo_list (arglist);
      bebop_exit (2);
    }
    if (got_arg_p (arglist, 'v')) { printf ("done, in %g seconds\n", seconds); }
  }
  const char* format = format_to_string(A->format);
  fprintf (stdout, "Current format: %s",format);

  fprintf (stdout, "\n===== Converting to CSR =====\n");
  if(A->format == COO) {
    coo_to_csr_convert(A);
  }
  else {
    sp_convert(A, "CSR");
  }

  struct csr_matrix_t* repr = A->repr;
  const char* type = value_type_to_string(repr->value_type);
  const char* sym_type = symmetry_type_to_string(repr->symmetry_type);

  fprintf (stdout, "m = %d, n = %d, nnz = %d, density = %f, val type = %s%s",
    repr->m,repr->n,repr->nnz,
    (float)repr->nnz/(float)(repr->m * repr->n),
    type,sym_type);
    
  FILE *LambdaFile;
  LambdaFile = fopen("lambda.txt", "a");
  assert(LambdaFile != NULL);
  fprintf(LambdaFile, "%s",opts.input_filename);
  fprintf(LambdaFile, " ,");
  fclose(LambdaFile);
  
  // ********** Run FENNEL ***************************************
  fprintf (stdout, "\n===== Running fennel =====\n");
  run_fennel(repr, parts, 1.5); //todo: nparts, gamma as inputs
  // *************************************************************
  //errcode = save_sparse_matrix ("out.mtx", A, MATRIX_MARKET);
  destroy_sparse_matrix (A);
  return 0;
}


static struct csr_matrix_t* mat_to_csr (struct coo_matrix_t* A) {
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

static void sort_coo (void* coord_array, const int length, enum value_type_t value_type) {
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

static int coo_to_csr_convert(struct sparse_matrix_t* A) {
  assert(A->format == COO);
  struct csr_matrix_t* B = mat_to_csr (A->repr);
  if (B == NULL) { return -1; }
  destroy_coo_matrix (A->repr);
  A->repr = B; A->format = CSR;
  return 0;
}
