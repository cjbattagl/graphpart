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
  fprintf (out, "EX: ./fennel -v -e 'test.mtx' 'MM' 4 10000\n");
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

  if (argc - optind != 4) {
      fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
         "e specified, but there should be %d ***\n", argc - optind, 4);
      dump_usage (stderr, argv[0], arglist, NULL);
      bebop_exit (EXIT_FAILURE); /* stops logging */
  }

  opts.input_filename = argv[optind];
  opts.input_file_format = argv[optind+1];
  int parts = atoi(argv[optind+2]);
  int cutoff = atoi(argv[optind+3]);
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
  run_fennel(repr, parts, 1.5, cutoff); //todo: nparts, gamma as inputs
  // *************************************************************
  //errcode = save_sparse_matrix ("out.mtx", A, MATRIX_MARKET);
  destroy_sparse_matrix (A);
  return 0;
}

