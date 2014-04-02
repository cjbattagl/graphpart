#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>

// Bebop utils
#include <bebop/util/merge_sort.h>
#include <bebop/util/config.h>
#include <bebop/util/get_options.h>
#include <bebop/util/init.h>
#include <bebop/util/log.h>
#include <bebop/util/malloc.h>
#include <bebop/util/timer.h>
#include <bebop/util/util.h>

// Bebop smc utils
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/read_mm.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/interface.h>
#include <bebop/smc/sparse_matrix_ops.h>

#include <metis.h>

#include "fennel_kernel.h"
#include "randperm.h"
#include "util.h"

// Given the CSR vectors of a sparse matrix, and a partition vector, creates
// a streaming-partitioning assignment of nodes to partitions
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

int deg_kernel(int n, int nparts, int *partsize, int *rowptr, int *colidx,
    bool **parts, float alpha, float gamma, int *emptyverts, int cutoff) {
      
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
      if(nnz_row < cutoff) {  //only partition low degrees
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
        if (parts[s][vert] == 1) { oldpart = s; }
        parts[s][vert] = 0; 
      }
      parts[best_part][vert] = 1;
      partsize[best_part]++;
      if (oldpart >= 0) { partsize[oldpart]--; }  
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

// Modified FENNEL algorithm that samples 'prob' percent of the edges
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

// Prepare matrix A for partitioning, and call FENNEL
// then evaluate results
int run_fennel(const struct csr_matrix_t* A, int nparts, float gamma, int cutoff) {
  int m, n, nnz;
  void* values;
  int* colidx;
  int* rowptr;
  enum symmetry_type_t symmetry_type;
  enum symmetric_storage_location_t symmetric_storage_location;
  enum value_type_t value_type;
  double seconds = 0.0;
  //int cutoff = 30;

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
  //fennel_kernel(n, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts);
  deg_kernel(n, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts, cutoff);
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
  
  cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, NULL, cutoff);
    
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
  int numruns = 1;
  for (int run=0; run<numruns; run++) {
  
    seconds = get_seconds();
    //fennel_kernel(n, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts);
    deg_kernel(n, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts, cutoff);
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
      cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, NULL, cutoff);
    }
    else {
      cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, PartFile, cutoff);
    }
    fprintf (stdout, "Percent edges cut = %d / %d = %1.3f", cutedges,nnz,(float)cutedges/nnz);
    max_load = partsize[0];
    min_load = partsize[0];
    for (s = 1; s < nparts; s++) {
      if (partsize[s] > max_load) {max_load = partsize[s];}
      if (partsize[s] < min_load) {min_load = partsize[s];}
    }
    fprintf (stdout, "       ----> Load balance: %d / %d = %1.3f\n",max_load,min_load,(float)max_load/min_load);
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
  
// Compute number of edges cut by a given partitioning
static int compute_cut(int *emptyparts, int *redparts, int *rowptr, int *colidx, bool **parts, int nparts, int n, FILE *out, int cutoff) {
  int vert, nnz_row, v_part;
  int cutedges = 0;
  int numnodes = 0;
  int totnodes = 0;
  int *row;
  for (int i = 0; i < n; i++) {
    vert = i;
    row = &rowptr[vert];
    nnz_row = *(row+1) - *(row); //nnz in row
    if (nnz_row < cutoff) { numnodes++; }
    totnodes++;
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
      if (nnz_row < cutoff && rowptr[colidx[k]+1] - rowptr[colidx[k]] < cutoff && parts[v_part][colidx[k]] != 1) { cutedges++; }
    }
  }
  fprintf (stdout, "Pct nodes below cutoff: %f\n",(float)numnodes/totnodes);
  return cutedges;
}

static float calc_dc(float alpha, float gamma, int len) {
  return (alpha*pow(len+0.5,gamma)) - (alpha*pow(len,gamma));
}
