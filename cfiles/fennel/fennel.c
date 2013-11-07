/**
 * Driver program for our FENNEL implementation
 * Note that this code is a modification of driver code from Bebop,
 * provided with the following LICENSE:
 *
 * Copyright (c) 2008, Regents of the University of California 
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright 
 *   notice, this list of conditions and the following disclaimer in 
 *   the documentation and/or other materials provided with the 
 *   distribution.
 *
 * * Neither the name of the University of California, Berkeley, nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior
 *   written permission.  
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

#include "fennel.h"
#include "randperm.h"

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
  fprintf (out, "fennel <in-filename> <in-format> <out-filename> <out-format> [options]\n");
  fprintf (out, "<in-filename>:   name of file containing the sparse matrix to read in\n");
  fprintf (out, "<in-format>:     format of the input file (\"HB\", \"ML\", \"MM\")\n");
  fprintf (out, "<out-filename>:  name of file to which to output results\n");
  fprintf (out, "[options]: -a -- validate the input matrix only, without outputting anything\n");
  fprintf (out, "           -e -- expand symmetric into unsymmetric storage\n");
  fprintf (out, "           -v  -- verbose mode\n");
  fprintf (out, "EX: ./fennel -a 'test.mtx' 'MM'\n");
  fprintf (out, "EX: ./fennel -v 'test.mtx' 'MM'\n");
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
  /* Set the get_options usage function to "usage", instead of using the default
   * usage function.  This is necessary because the command-line arguments include
   * things that are not "options" in the strict sense, because they do not follow
   * a "-[a-z]".
   */
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

  if (got_arg_p (arglist, 'a')) {
      int errcode = do_validate_matrix (argc, argv, arglist);
      destroy_arginfo_list (arglist);
      deinit_timer();
      bebop_exit (errcode); /* stops logging */
  }

  if (argc - optind != 3) {
      fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
	       "e specified, but there should be %d ***\n", argc - optind, 3);
      dump_usage (stderr, argv[0], arglist, NULL);
      bebop_exit (EXIT_FAILURE); /* stops logging */
  }

  opts.input_filename = argv[optind];
  opts.input_file_format = argv[optind+1];
  opts.output_filename = argv[optind+2];
  //opts.output_file_format = argv[optind+3];

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


  fprintf (stdout, "\n===== Converting to CSR =====\n");
  sp_convert(A, "CSR");
  assert(A->format == CSR);
  
  struct csr_matrix_t* repr = A->repr;
  
  fprintf (stdout, "m = %d, n = %d, nnz = %d, density = %f\n",repr->m,repr->n,repr->nnz,
    (float)repr->nnz/(float)(repr->m * repr->n));
  
  // ********** Run FENNEL ***************************************
  fprintf (stdout, "\n===== Running fennel =====\n");
  run_fennel(repr, 2, 1.4); //todo: nparts, gamma as inputs
  // *************************************************************
  
  destroy_sparse_matrix (A);
  return 0;
}


static int run_fennel(const struct csr_matrix_t* A, int nparts, float gamma) {
  int m, n, nnz;
  void* values;
  int* colidx;
  int* rowptr;
  enum symmetry_type_t symmetry_type;
  enum symmetric_storage_location_t symmetric_storage_location;
  enum value_type_t value_type;
  
  double seconds = 0.0;

  unpack_csr_matrix (A, &m, &n, &nnz, &values, &colidx, &rowptr, &symmetry_type,
		   &symmetric_storage_location, &value_type);
		  
  // Set alpha
  float alpha = sqrt(2) * (nnz/pow(n,gamma));
  // float alpha = nnz * pow(2,(gamma-1))/pow(m,gamma);
  fprintf (stdout, "----> gamma = %f, alpha = %f\n",gamma,alpha);
  	   
  // Generate vorder
  fprintf (stdout, "----> Gen rand perm of size %d\n",n);
  int *vorder = genRandPerm(n);
  
  // Allocate partition vectors
  fprintf (stdout, "----> Gen %d partition vectors of size %d\n",nparts,n);
  bool** parts = (bool**)malloc(nparts * sizeof(bool*));
  for(int i = 0; i < nparts; i++) { 
    parts[i] = (bool*)malloc(n * sizeof(bool));
    for( int j=0; j<n; j++ ) { parts[i][j] = 0; } //fill with zeros 
    assert(parts[i]);
  }
  assert(parts);
  
  // Iterate over vorder
  fprintf (stdout, "----> Begin outer loop . . . \n");
  int v, k, s, node;
  int *row;
  int nnz_row;
  int *partscore = (int*)malloc(nparts * sizeof(int));
  int *partsize = (int*)malloc(nparts * sizeof(int));
  
  //initialize part sizes to 0
  for (s = 0; s < nparts; s++) { partsize[s] = 0; }

  int best_part;
  float curr_score, best_score;
  
  float c1, c2;
  float dc;
  
  seconds = get_seconds();
  // iterate through nodes in vorder
  for (int i = 0; i < n; i++) {
    for (s = 0; s < nparts; s++) { partscore[s]=0; }
    v = vorder[i];
    row = &rowptr[v];
    nnz_row = *(row+1) - *(row); //nnz in row
    
    // generate partition scores for each partition
    for (k = *row; k < ((*row)+nnz_row); k++) {
      node = colidx[k];
      for (s = 0; s < nparts; s++) { if (parts[s][node] == 1) { partscore[s]++; break; }}
    }
    
    // choose optimal partition (initializing first)
    best_score = partscore[0] - ((alpha*pow(partsize[0]+1,gamma)) - (alpha*pow(partsize[0],gamma)));
    best_part = 0;
    
    for (s = 1; s < nparts; s++) {
      curr_score = partscore[s] - ((alpha*pow(partsize[s]+1,gamma)) - (alpha*pow(partsize[s],gamma)));
      if (curr_score > best_score) { best_score = curr_score; best_part = s; }
    }
    
    parts[best_part][v] = 1; partsize[best_part]++;
  }
  seconds = get_seconds() - seconds;

 
  // Compute load balance
  int max_load = partsize[0];
  int min_load = partsize[0];
  for (s = 1; s < nparts; s++) {
    if (partsize[s] > max_load) {max_load = partsize[s];}
    if (partsize[s] < min_load) {min_load = partsize[s];}
  }

  fprintf (stdout, "\n===== Fennel Complete in %g seconds =====\n", seconds);
  fprintf (stdout, "----> Partition sizes: ");
  for (s = 0; s < nparts; s++) {
    fprintf (stdout, "| %d |", partsize[s]);
  }
  fprintf (stdout, "\n----> Load balance: %d / %d = %f\n",max_load,min_load,(float)max_load/min_load);

  // Compute cut quality
  int cutedges = 0;
  int v_part = 0;
  
  for (int i = 0; i < n; i++) {
    v = i;
    row = &rowptr[v];
    nnz_row = *(row+1) - *(row); //nnz in row
    
    // find v's partition
    for (s = 0; s < nparts; s++) {
      if (parts[s][node] == 1) { v_part = s; break; }
      //assert(0); //We shouldn't reach this if every vertex has a partition
    }
    
    // count edges to other partitions
    for (k = *row; k < ((*row)+nnz_row); k++) {
      if (parts[v_part][colidx[k]] != 1) { cutedges++; }
    }
  }
  
  fprintf (stdout, "\n----> Percent edges cut = %d / %d = %f\n",cutedges,nnz,(float)cutedges/nnz);
}

/**
 * Perform the matrix validation operation specified by the "-a" command-line option.
 */
static int
do_validate_matrix (int argc, char *argv[], struct arginfo* arglist)
{
  extern int optind;

  enum sparse_matrix_file_format_t informat = 0;
  struct sparse_matrix_t* A = NULL;
  double seconds = 0.0;
  int errcode = 0;

  if (argc - optind != 2) {
      fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
	       "e specified, but there should be %d ***\n", argc - optind, 2);
      dump_usage (stderr, argv[0], arglist, NULL);
      return -1;
  }

  opts.input_filename = argv[optind];
  opts.input_file_format = argv[optind+1];

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
      return -1;
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
      return 1;
  }
  if (got_arg_p (arglist, 'v')) { printf ("done, in %g seconds\n", seconds); }
  
  if (got_arg_p (arglist, 'v')) {
      printf ("Validating sparse matrix...");
      fflush (stdout); 
  }
  seconds = get_seconds();
  errcode = valid_sparse_matrix (A);
  seconds = get_seconds() - seconds;
  if (got_arg_p (arglist, 'v')) { printf ("done, in %g seconds\n", seconds); }

  if (valid_sparse_matrix (A)) { printf ("\n### Sparse matrix is valid ###\n\n"); }
  else { printf ("\n### Invalid sparse matrix! ###\n\n"); }

  destroy_sparse_matrix (A);
  return 0;
}