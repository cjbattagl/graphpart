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

int fennel_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx,
    bool **parts, float alpha, float gamma, int *emptyverts) {
      
  int *partscore = (int*)malloc(nparts * sizeof(int));
  int *row;
  int vert, k, s, nnz_row, best_part, randidx, nededges = 0, node = 0;
  float curr_score, best_score;
  int *vorder = genRandPerm(n);
  int oldpart;

  emptyverts = 0;

  for (int i = 0; i < n; i++) {
    for (s = 0; s < nparts; s++) { partscore[s]=0; }
    vert = vorder[i];
    row = &rowptr[vert];
    nnz_row = *(row+1) - *row;
    oldpart = -1;
   if(nnz_row != 0) {
      // generate partition scores for each partition
      for (k = *row; k < ((*row)+nnz_row); k++) {
        node = colidx[k];
        for (s = 0; s < nparts; s++) { if (parts[s][node] == 1) { partscore[s]++; /*break;*/ }}
      }
        
      // choose optimal partition (initializing first)
      best_score = (partscore[0]-nnz_row) - calc_dc(alpha,gamma,partsize[0]);
      best_part = 0;
      for (s = 1; s < nparts; s++) {
        curr_score = (partscore[s]-nnz_row) - calc_dc(alpha,gamma,partsize[s]);
        if (curr_score > best_score) { best_score = curr_score; best_part = s; }
      }
      for (s = 0; s < nparts; s++) { 
        if (parts[s][vert] == 1) {
          oldpart = s;
        }
        parts[s][vert] = 0; 
      }
      parts[best_part][vert] = 1;
      //int sum=0;
      //for (s = 0; s < nparts; s++) { sum += parts[s][vert]; }
      //assert(sum==1);
      partsize[best_part]++;
      if (oldpart >= 0) {
        partsize[oldpart]--;
      }
      
    } else { // empty vertex for some reason... assign it to random permutation
      emptyverts++;
      randidx = irand(nparts);
      for (s = 1; s < nparts; s++) {
        if (parts[s][vert] == 1) {
          oldpart = s;
        }
        parts[s][vert] = 0; 
      }
      parts[randidx][vert] = 1;
      partsize[randidx]++;
      if (oldpart >= 0) {
        partsize[oldpart]--;
      }
    }
  }
  
  free(vorder);
}

int sample_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx,
    bool **parts, float alpha, float gamma, int *emptyverts, float prob) {
      
  int *partscore = (int*)malloc(nparts * sizeof(int));
  int *row;
  int vert, k, s, nnz_row, best_part, randidx, nededges = 0, node = 0;
  float curr_score, best_score;
  int *vorder = genRandPerm(n);
  int oldpart;

  emptyverts = 0;

  for (int i = 0; i < n; i++) {
    for (s = 0; s < nparts; s++) { partscore[s]=0; }
    vert = vorder[i];
    row = &rowptr[vert];
    nnz_row = *(row+1) - *row;
    int *randOrder = genRandPerm(nnz_row);

    oldpart = -1;
   if(nnz_row != 0) {
      // generate partition scores for each partition
      k = *row;
      
      int end = floor(prob*nnz_row);
      int end_idx = (nnz_row > 200) ? end : nnz_row;
      
      for (int j=0; j<end_idx; j++) {
        node = colidx[k + randOrder[j]];
        for (s = 0; s < nparts; s++) { if (parts[s][node] == 1) { partscore[s]++; /*break;*/ }}
      }
      
      // choose optimal partition (initializing first)
      best_score = (partscore[0]-nnz_row) - calc_dc(alpha,gamma,partsize[0]);
      best_part = 0;
      for (s = 1; s < nparts; s++) {
        curr_score = (partscore[s]-nnz_row) - calc_dc(alpha,gamma,partsize[s]);
        if (curr_score > best_score) { best_score = curr_score; best_part = s; }
      }
      for (s = 0; s < nparts; s++) { 
        if (parts[s][vert] == 1) {
          oldpart = s;
        }
        parts[s][vert] = 0; 
      }
      parts[best_part][vert] = 1;
      //int sum=0;
      //for (s = 0; s < nparts; s++) { sum += parts[s][vert]; }
      //assert(sum==1);
      partsize[best_part]++;
      if (oldpart >= 0) {
        partsize[oldpart]--;
      }
      
    } else { // empty vertex for some reason... assign it to random permutation
      emptyverts++;
      randidx = irand(nparts);
      for (s = 1; s < nparts; s++) {
        if (parts[s][vert] == 1) {
          oldpart = s;
        }
        parts[s][vert] = 0; 
      }
      parts[randidx][vert] = 1;
      partsize[randidx]++;
      if (oldpart >= 0) {
        partsize[oldpart]--;
      }
    }
    free(randOrder);
  }
  
  free(vorder);
}

int run_fennel(const struct csr_matrix_t* A, int nparts, float gamma) {
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

  float alpha = sqrt(2) * (nnz/pow(n,gamma));
  //float alpha = nnz * pow(2,(gamma-1))/pow(m,gamma);
  fprintf (stdout, "----> gamma = %f, alpha = %f\n",gamma,alpha);
          
  // Allocate partition vectors
  fprintf (stdout, "----> Gen %d partition vectors\n",nparts);
  bool** parts = (bool**)malloc(nparts * sizeof(bool*));
  for(int i = 0; i < nparts; i++) {
    parts[i] = (bool*)malloc(n * sizeof(bool));
    for( int j=0; j<n; j++ ) { parts[i][j] = 0; } //fill with zeros
    assert(parts[i]);
  }
  assert(parts);

  // Iterate over vorder
  fprintf (stdout, "----> Begin outer loop... \n");
  int vert, k, s, node = 0;
  int *row;
  int nnz_row = 0;
  int emptyverts;
  
  //initialize part sizes
  int *partsize = (int*)malloc(nparts * sizeof(int));
  for (s = 0; s < nparts; s++) { partsize[s] = 0; }

  seconds = get_seconds();
  fennel_kernel(n, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts);
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
  for (s = 0; s < nparts; s++) { fprintf (stdout, "| %d |", partsize[s]); }
  fprintf (stdout, "\n----> Load balance: %d / %d = %1.3f\n",max_load,min_load,(float)max_load/min_load);

  // Compute cut quality
  int cutedges = 0;
  int emptyparts = 0; //empty assignments
  int redparts = 0; //redundant assignments
  
  cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, NULL);
    
  fprintf (stdout, "----> Percent edges cut = %d / %d = %1.3f\n",cutedges,nnz,(float)cutedges/nnz);
  fprintf (stdout, "----> Percent of random: %1.3f\n\n",((float)cutedges/nnz)/((float)(nparts-1)/nparts));
  fprintf (stdout, "===== Sanity Check =====\n");
  fprintf (stdout, "----> Unassigned vertices (error): %d\n",emptyparts);
  fprintf (stdout, "----> Overassigned vertices (error): %d\n", redparts);
  fprintf (stdout, "----> Empty vertices: %d\n", emptyverts);
  
  FILE *PartFile;
  PartFile = fopen("parts.mat", "w");
  assert(PartFile != NULL);
  
  FILE *LambdaFile;
  LambdaFile = fopen("lambda.txt", "a");
  assert(LambdaFile != NULL);

  
  /////////////////////////////////////////////////////////////////////////////
  ///////////// EXPERIMENTAL: DO ADDITIONAL RUNS ON THE NEW PARTITION /////////
  /////////////////////////////////////////////////////////////////////////////

  for (s = 0; s < nparts; s++) { partsize[s] = 0; }
  for(int i = 0; i < nparts; i++) {
      for( int j=0; j<n; j++ ) { parts[i][j] = 0; } //fill with zeros
  }
  int numruns = 2;
  for (int run=0; run<numruns; run++) {
  
    seconds = get_seconds();
    fennel_kernel(n, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts);
    seconds = get_seconds() - seconds;
   //fprintf (stdout, "\n %g seconds =====\n", seconds);


    // Compute load balance  
    max_load = partsize[0];
    min_load = partsize[0];
    for (s = 1; s < nparts; s++) {
      if (partsize[s] > max_load) {max_load = partsize[s];}
      if (partsize[s] < min_load) {min_load = partsize[s];}
    }

    //fprintf (stdout, "----> Run |%d|: Load balance: %d / %d = %1.2f\t",run,max_load,min_load,(float)max_load/min_load);
    cutedges = 0;
    emptyparts = 0; //empty assignments
    redparts = 0; //redundant assignments
    
    if (run < numruns-1) {
      cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, NULL);
    }
    else {
      cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, PartFile);
    }
    fprintf (stdout, "Percent edges cut = %d / %d = %1.3f\n", cutedges,nnz,(float)cutedges/nnz);
    fprintf(LambdaFile, "%1.3f, ",(float)cutedges/nnz);

    if (run == numruns-1) {
    //fprintf (stdout, "----> Unassigned vertices (error): %d\n",emptyparts);
    //fprintf (stdout, "----> Overassigned vertices (error): %d\n", redparts);
    //fprintf (stdout, "----> Empty vertices: %d\n", emptyverts);
    //fprintf (stdout, "----> Percent of random: %1.3f\n\n",((float)cutedges/nnz)/((float)(nparts-1)/nparts));
    }
  }
  fclose(PartFile);
      fprintf(LambdaFile, "\n");

  fclose(LambdaFile);

  //////////////////////////////////
  //////// METIS RUN ///////////////
  //////////////////////////////////

  idx_t *xadj, *adjncy, *vwgt, *adjwgt, *part, *perm, *iperm, *sep;
  idx_t options[METIS_NOPTIONS] = {0}, edgecut, sepsize;
  
  struct parameter_data params = {
    0, // wgtflag = 0
    0, // adjwgt = 0
    NULL, // vsize = NULL
  };
    
  double *optarray, *partpr, *permpr, *ipermpr, *seppr;
  METIS_SetDefaultOptions(options);
  //csr_to_metis (n, nnz, rowptr, colidx, &xadj, &adjncy, &vwgt, &adjwgt);

  // Allocate vsize
  params.vsize = (idx_t*)malloc(sizeof(idx_t)*n);

  // Set which METIS function to do here, or take in command line
  const char* funcname = "PartGraphRecursive";

  if (strcasecmp(funcname,"PartGraphRecursive")==0 ||
    strcasecmp(funcname,"PartGraphKway")==0 ) {

/*  // Figure out addl options:
    parseOptions(opt, options, &params);
    if (params.wgtflag == 0) {
      for (i=0; i<n; ++i) { vwgt[i] = 1; }
    }
        
    if (params.adjwgt == 0) {
      for (i=0; i<xadj[n]; ++i) { adjwgt[i] = 1; }
    }

    if (nparts < 2) {
      mexErrMsgTxt("nparts must be at least 2");
    }
*/
/*
    // Allocate memory for result of call
    part = (idx_t*) malloc (n*sizeof(idx_t));
        
    idx_t ncon = 1;
    for (int i=0; i<n; ++i) {
      params.vsize[i] = 1;
    }
    
    idx_t n_i = (idx_t)n;
    idx_t nparts_i = (idx_t)nparts;
    
    // Do the call
    seconds = get_seconds();
    if (strcasecmp(funcname,"PartGraphRecursive") == 0) {
      checkCall(METIS_PartGraphRecursive(&n_i, &ncon, xadj, adjncy, vwgt, params.vsize, adjwgt, &nparts_i, NULL, NULL, options, &edgecut, part)); 
    } 
    else if (strcasecmp(funcname, "PartGraphKway") == 0) {
      checkCall(METIS_PartGraphKway (&n_i, &ncon, xadj, adjncy, vwgt, params.vsize, adjwgt, &nparts_i, NULL, NULL, options, &edgecut, part));
    }
    else {
      fprintf(stdout,"METIS: unhandled case\n");
      return -1;
    }
    seconds = get_seconds() - seconds;
    fprintf (stdout, "\n===== METIS Complete in %g seconds =====\n", seconds);
    fprintf (stdout, "\tMETIS edges cut = %d / %d = %1.3f\n",(int)edgecut,nnz,(float)edgecut/nnz);*/
  }
}

static void csr_to_metis (int n, int nnz, int *rowptr, int *colidx, idx_t **xadj, idx_t **adjncy, idx_t **vwgt, idx_t **adjwgt) {
    int i, j, jbar;
    /* Allocate room for METIS's structure */
    *xadj = (idx_t*) malloc (n+1 * sizeof(idx_t));
    *adjncy = (idx_t*) malloc (nnz * sizeof(idx_t));
    *vwgt = (idx_t*) malloc (n * sizeof(idx_t));
    *adjwgt = (idx_t*) malloc (nnz * sizeof(idx_t));

    (*xadj)[0] = 0;
    jbar = 0;
    int nnz_row;
    for (i = 1; i <= n; i++) {
        nnz_row = rowptr[i] - rowptr[i-1];

        for (j = rowptr[i-1]; j < rowptr[i]; j++) {
            if (colidx[j] != i-1) {
                //fprintf(stdout,"%d ",jbar);
                (*adjncy)[jbar] = colidx[j];
                (*adjwgt)[jbar] = 1;
                jbar++;
            } else {
                (*vwgt)[i-1] = 1;
            }
        }
        (*xadj)[i] = jbar;
    }
}
  
int compute_cut(int *emptyparts, int *redparts, int *rowptr, int *colidx, bool **parts, int nparts, int n, FILE *out) {
  int vert, nnz_row, v_part;
  int cutedges = 0;
  int *row;
  for (int i = 0; i < n; i++) {
    vert = i;
    row = &rowptr[vert];
    nnz_row = *(row+1) - *(row); //nnz in row
    v_part = -1;
    
    // find v's partition
    for (int s = 0; s < nparts; s++) {
      if (parts[s][vert] == 1) {
        if(v_part != -1) { redparts++; }
        v_part = s;
      }
    }
    if (v_part == -1) {
      v_part = 0;
      emptyparts++;
    }
    
    if (out != NULL) {
      fprintf (out, "%d %d\n",i+1,v_part+1);
    }
    
    // count edges to other partitions
    for (int k = *row; k < ((*row)+nnz_row); k++) {
      if (parts[v_part][colidx[k]] != 1) { cutedges++; }
    }
  }
  
  return cutedges;
}

static float calc_dc(float alpha, float gamma, int len) {
  return (alpha*pow(len+0.5,gamma)) - (alpha*pow(len,gamma));
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

