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

//#include <metis.h>

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

// Prepare matrix A for partitioning, and call FENNEL
// then evaluate results
int run_fennel(const struct csr_matrix_t* A, int nparts, float gamma) {
  int m, n, nnz;
  void* values;
  int* colidx;
  int* rowptr;
  enum symmetry_type_t symmetry_type;
  enum symmetric_storage_location_t symmetric_storage_location;
  enum value_type_t value_type;
  double seconds = 0.0;
  int cutoff = 300000;

  unpack_csr_matrix (A, &m, &n, &nnz, &values, &colidx, &rowptr, &symmetry_type,
                 &symmetric_storage_location, &value_type);

  float alpha = sqrt(2) * (nnz/pow(n,gamma));
  //float alpha = nnz * pow(2,(gamma-1))/pow(m,gamma);
  fprintf (stdout, "----> gamma = %f, alpha = %f\n",gamma,alpha);
          
  // Allocate partition vectors
  bool** parts = (bool**)malloc(nparts * sizeof(bool*));
  for(int i = 0; i < nparts; i++) {
    parts[i] = (bool*)malloc(n * sizeof(bool));
    for( int j=0; j<n; j++ ) { parts[i][j] = 0; } //fill with zeros
    assert(parts[i]);
  }
  assert(parts);

  // Iterate over vorder
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

  fprintf (stdout, "\n===== First Fennel Run Complete in %g seconds =====\n", seconds);
  // fprintf (stdout, "----> Partition sizes: ");
  // for (s = 0; s < nparts; s++) { fprintf (stdout, "| %d |", partsize[s]); }
  // fprintf (stdout, "\n----> Load balance: %d / %d = %1.3f\n",max_load,min_load,(float)max_load/min_load);

  // Compute cut quality
  int cutedges = 0;
  int emptyparts = 0; //empty assignments
  int redparts = 0; //redundant assignments
  
  cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, NULL);
    
  // fprintf (stdout, "----> Percent edges cut = %d / %d = %1.3f\n",cutedges,nnz,(float)cutedges/nnz);
  // fprintf (stdout, "----> Percent of random: %1.3f\n\n",((float)cutedges/nnz)/((float)(nparts-1)/nparts));
  // fprintf (stdout, "===== Sanity Check =====\n");
  // fprintf (stdout, "----> Unassigned vertices (error): %d\n",emptyparts);
  // fprintf (stdout, "----> Overassigned vertices (error): %d\n", redparts);
  // fprintf (stdout, "----> Empty vertices: %d\n", emptyverts);
  
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
  int numruns = 20;
  for (int run=0; run<numruns; run++) {
  
    seconds = get_seconds();
    fennel_kernel(n, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts);
    // deg_kernel(n, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts, cutoff);
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
    
    // if (run < numruns-1) {
    //   //cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, NULL, cutoff);
    // }
    // else {
      cutedges = compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, PartFile);
    // }
    fprintf (stdout, "Run [%d]: Percent edges cut = %d / %d = %1.3f \n", run,cutedges,nnz,(float)cutedges/nnz);
    max_load = partsize[0];
    min_load = partsize[0];
    for (s = 1; s < nparts; s++) {
      if (partsize[s] > max_load) {max_load = partsize[s];}
      if (partsize[s] < min_load) {min_load = partsize[s];}
    }
    fprintf(LambdaFile, "%1.3f, ",(float)cutedges/nnz);

    if (run == numruns-1) {
    fprintf (stdout, "----> Load balance: %d / %d = %1.3f\n",max_load,min_load,(float)max_load/min_load);
    fprintf (stdout, "----> Unassigned vertices (error): %d\n",emptyparts);
    fprintf (stdout, "----> Overassigned vertices (error): %d\n", redparts);
    fprintf (stdout, "----> Empty vertices: %d\n", emptyverts);
    fprintf (stdout, "----> Percent of random: %1.3f\n\n",((float)cutedges/nnz)/((float)(nparts-1)/nparts));
    }
  }
  fclose(PartFile);
      fprintf(LambdaFile, "\n");

  fclose(LambdaFile);
}

// Compute number of edges cut by a given partitioning
static int compute_cut(int *emptyparts, int *redparts, int *rowptr, int *colidx, bool **parts, int nparts, int n, FILE *out) {
  int vert, nnz_row, v_part;
  int nnz = 0;
  int cutedges = 0;
  int tot_lo_deg_edges = 0;
  int numnodes = 0;
  int totnodes = 0;
  int *row;

  size_t *cuts_per_part = (size_t*)malloc(nparts * sizeof(size_t));

  for (int i=0; i<nparts; i++) {
    cuts_per_part[i] = 0;
  }

  for (int i = 0; i < n; i++) {
    vert = i;
    row = &rowptr[vert];
    nnz_row = *(row+1) - *(row); //nnz in row
    numnodes++;
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
      v_part = nparts;
      emptyparts++;
    }
    
    if (out != NULL) { fprintf (out, "%d %d\n",i+1,v_part+1); }
    
    // count edges to other partitions
    for (int k = *row; k < ((*row)+nnz_row); k++) {
      nnz++;
      if (parts[v_part][colidx[k]] != 1) { 
        cuts_per_part[v_part]++;
        cutedges++; 
      }
    }
  }
  //fprintf (stdout, "%d / %d = \t",cutedges, tot_lo_deg_edges);
  // fprintf (stdout, "%f\t",(float)cutedges/tot_lo_deg_edges);
  //fprintf (stdout, "Pct nnz below cutoff: %f\n",(float)tot_lo_deg_edges/nnz);
  // fprintf (stdout, "%f\n",(float)numnodes/totnodes);
  for (int i=0; i<nparts; i++) {
      // fprintf (stdout, "%d\n",cuts_per_part[i]);
  }
  return cutedges;
}

static float calc_dc(float alpha, float gamma, int len) {
  return (alpha*pow(len+10,gamma)) - (alpha*pow(len,gamma));
}
