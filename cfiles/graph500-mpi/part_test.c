/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <mpi.h>
#include "common.h"
#include "../generator/make_graph.h"
#include "../generator/apply_permutation_mpi.h"

#define isPowerOfTwo(x) ((x != 0) && !(x & (x - 1)))
#define TO_NEW_IDX(x) g_perm[x]
#define NEW_PART_OF_IDX(x) parts[TO_NEW_IDX(x)]
#define VERTEX_TO_GLOBAL_ALT(r, i) ((int64_t)(MUL_SIZE((uint64_t)i) + (int)(r))) //Todo...
#define MUL_SIZE(x) ((x) * size)

static int* g_perm;
static int* parts;
/*static int num_total_lo_deg_nodes = 0;
static int num_local_lo_deg_nodes = 0;
static int num_total_hi_deg_nodes = 0;
static int num_local_hi_deg_nodes = 0;
static int local_lo_deg_offset = 0;*/
void apply_permutation_mpi(MPI_Comm comm, const int64_t local_perm_size, const int64_t* const local_vertex_perm, const int64_t N, const int64_t nedges, int64_t* result);

int print_parts(FILE* out, int* parts, int n, int n_local);
int print_graph_tuple(FILE* out, int64_t* result, int64_t n_local, int rank);
int print_graph_csr(FILE* out, size_t *rowptr, int64_t *colidx, int n_local);
int mpi_compute_cut(size_t *rowptr, int64_t *colidx, int *parts, int nparts, int n_local, int offset, int cutoff, csr_graph* const g);
void partition_graph_data_structure(csr_graph* const g);
int* genRandPerm(int* orderList, int size);
void shuffle_int(int *list, int len);
int irand(int n);
float calc_dc(float alpha, float gamma, int len);

enum {s_minimum, s_firstquartile, s_median, s_thirdquartile, s_maximum, s_mean, s_std, s_LAST};
static void get_statistics(const double x[], int n, double r[s_LAST]);
  static int compare_doubles(const void* a, const void* b);



int main(int argc, char* argv[]) {
  //int size, rank;
  double start, stop;
  size_t i;
//////
  MPI_Init(&argc, &argv);
  setup_globals();

  int SCALE = 16;
  double edgefactor = 16.; /* nedges / nvertices, i.e., 2*avg. degree */
  if (argc >= 2) SCALE = atoi(argv[1]);
  if (argc >= 3) edgefactor = atof(argv[2]);
  if (argc <= 1 || argc >= 4 || SCALE == 0 || edgefactor == 0) {
    if (rank == 0) {
      fprintf(stderr, "Usage: %s SCALE edgefactor\n  SCALE = log_2(# vertices) [integer, required]\n  edgefactor = (# edges) / (# vertices) = .5 * (average vertex degree) [float, defaults to 16]\n(Random number seed and Kronecker initiator are in main.c)\n", argv[0]);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  uint64_t seed1 = 2, seed2 = 3;
  const double initiator[4] = {.57, .19, .19, .05};


  // 1. CREATE GRAPH
  int64_t nedges;
  int64_t* edges;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  make_graph(SCALE, (int64_t)(edgefactor * pow(2., SCALE)), seed1, seed2, initiator, &nedges, &edges);
  MPI_Barrier(MPI_COMM_WORLD);
  stop = MPI_Wtime();
  if (rank == 0) fprintf(stderr, "Graph size is %" PRIu64 " vertices and %" PRIu64 " edges\n", (uint64_t)pow(2., SCALE), nedges);


  if(1) { // Print graph
    char filename[256];
    sprintf(filename, "out_tup%02d.mat", rank);
    FILE *GraphFile;
    GraphFile = fopen(filename, "w");
    assert(GraphFile != NULL);
    print_graph_tuple(GraphFile, edges, nedges, rank);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  unsigned long my_edges = 0;
  unsigned long global_edges;

  for (i = 0; i < nedges; ++i) {
    if (edges[2 * i] != (uint64_t)(-1)) {
      ++my_edges;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&my_edges, &global_edges, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    fprintf(stderr, "%lu edge%s generated and permuted in %fs (%f Medges/s on %d processor%s)\n", global_edges, (global_edges == 1 ? "" : "s"), (stop - start), global_edges / (stop - start) * 1.e-6, size, (size == 1 ? "" : "s"));
  }

  // 2. CONVERT TO CSR
  csr_graph g;
  double data_struct_start = MPI_Wtime();
  convert_graph_to_csr(nedges, edges, &g);
  double data_struct_stop = MPI_Wtime();
  double data_struct_time = data_struct_stop - data_struct_start;
  if (rank == 0) { fprintf(stdout, "construction_time:              %g s\n", data_struct_time); }

  // 3. RUN FENNEL ON CSR
  partition_graph_data_structure(&g); 

  int64_t perm_local_size = 0;
  int64_t* local_vertex_perm = NULL;

  // 4. PERMUTE TUPLE GRAPH BASED ON COMPUTED PERMUTATION
  int* perms = (int*)malloc(g.nglobalverts * sizeof(int));
  for (i=0; i<g.nglobalverts; ++i) {
    int node = i;
    int node_owner = VERTEX_OWNER(node);
    int node_local_idx = VERTEX_LOCAL(node);
    int parts_idx = node_owner*g.nlocalverts + node_local_idx;
    int v_part = parts[parts_idx];
    perms[i] = v_part;
  }

  for (i=0; i < g.nglobalverts; ++i) {
    if (parts[i] == rank) { ++perm_local_size; }
  }

  // test perm size
  int perm_global_size;
  MPI_Allreduce(&perm_local_size, &perm_global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0) { fprintf(stdout, "Permutation size:              %d\n", perm_global_size); }
  assert(perm_global_size == g.nglobalverts);

  int64_t* perm_sizes = (int64_t*)xmalloc(size * sizeof(int64_t));
  MPI_Allgather((void*)&perm_local_size, 1, INT64_T_MPI_TYPE, perm_sizes, 1, INT64_T_MPI_TYPE, MPI_COMM_WORLD);
  assert(perm_local_size == perm_sizes[rank]);

  int64_t* perm_displs = (int64_t*)xmalloc((size + 1) * sizeof(int64_t));
  perm_displs[0] = 0;
  for (i = 1; i < size + 1; ++i) perm_displs[i] = perm_displs[i - 1] + perm_sizes[i - 1];

  local_vertex_perm = (int64_t*)malloc(perm_local_size * sizeof(int64_t));

  // this way codes for which vertex to receive for.
  // we want to code for which vertex to /send/ to
  /*
  for (i=0; i < g.nglobalverts; ++i) {
    if (perms[i] == rank) { 
      local_vertex_perm[perm_id] = i; 
      perm_id++;
      assert(perm_id <= perm_local_size);
    }
  }
  */
  int perm_id = 0;
  int64_t* perm_idxs = (int64_t*)xmalloc((size + 1) * sizeof(int64_t));
  perm_idxs[0] = 0;
  for (i = 1; i < size + 1; ++i) perm_idxs[i] = perm_displs[i];

  for (i=0; i < g.nglobalverts; ++i) {
    int part = perms[i];
    if (i >= perm_displs[rank] && i < perm_displs[rank+1]) {
      local_vertex_perm[perm_id] = perm_idxs[part];
      perm_id++;
    }
    perm_idxs[part]++;
  }
/*
  fprintf (stdout, "%d :",rank+1);

  for (i=0; i<perm_local_size; ++i) {
    fprintf (stdout, "[%d]%d  ", (int)perm_displs[rank]+i+1, (int)local_vertex_perm[i]+1);
  }
  fprintf (stdout, "\n");

  if(rank==0) {
    for (i=0; i<g.nglobalverts; ++i) {
    fprintf (stdout, "%d ", perms[i]+1);
  }
    fprintf (stdout, "\n");
}
*/

  apply_permutation_mpi(MPI_COMM_WORLD, perm_local_size, local_vertex_perm, g.nglobalverts, nedges, edges);
  if(1) { // Print graph
    char filename[256];
    sprintf(filename, "out_permed%02d.mat", rank);
    FILE *GraphFile;
    GraphFile = fopen(filename, "w");
    assert(GraphFile != NULL);
    print_graph_tuple(GraphFile, edges, nedges, rank);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // 5. CONVERT PERMUTED GRAPH TO CSR
  csr_graph g2;
  data_struct_start = MPI_Wtime();
  convert_graph_to_csr(nedges, edges, &g2);
  data_struct_stop = MPI_Wtime();
  data_struct_time = data_struct_stop - data_struct_start;
  if (rank == 0) { fprintf(stdout, "permuted csr construction_time:              %g s\n", data_struct_time); }
  if(1) { // Print graph
    char filename[256];
    sprintf(filename, "out_pcsr%02d.mat", rank);
    FILE *GraphFile;
    GraphFile = fopen(filename, "w");
    assert(GraphFile != NULL);
    print_graph_csr(GraphFile, g2.rowstarts, g2.column, g2.nlocalverts);
    MPI_Barrier(MPI_COMM_WORLD);
  }


  free_csr_graph(&g);






/* Get roots for BFS runs. */
  int num_bfs_roots = 2;
  int64_t* bfs_roots = (int64_t*)xmalloc(num_bfs_roots * sizeof(int64_t));
  find_bfs_roots(&num_bfs_roots, &g2, seed1, seed2, bfs_roots);

  /* Number of edges visited in each BFS; a double so get_statistics can be
   * used directly. */
  double* edge_counts = (double*)xmalloc(num_bfs_roots * sizeof(double));

  /* Run BFS. */
  int validation_passed = 1;
  double* bfs_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
  double* validate_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
  int64_t* pred = (int64_t*)xMPI_Alloc_mem(g.nlocalverts * sizeof(int64_t));

  int bfs_root_idx;
 
  for (bfs_root_idx = 0; bfs_root_idx < num_bfs_roots; ++bfs_root_idx) {
    int64_t root = bfs_roots[bfs_root_idx];

    if (rank == 0) fprintf(stderr, "Running BFS %d\n", bfs_root_idx);

    /* Clear the pred array. */
    memset(pred, 0, g2.nlocalverts * sizeof(int64_t));

    /* Do the actual BFS. */
    double bfs_start = MPI_Wtime();
    int64_t nvisited = 0;
    run_mpi_bfs(&g2, root, &pred[0], &nvisited);
    double bfs_stop = MPI_Wtime();
    bfs_times[bfs_root_idx] = bfs_stop - bfs_start;

    /* Validate result. */
    if (rank == 0) fprintf(stderr, "Validating BFS %d\n", bfs_root_idx);

    double validate_start = MPI_Wtime();
    int validation_passed_one = validate_bfs_result(&g2, root, pred, nvisited);
    double validate_stop = MPI_Wtime();
    validate_times[bfs_root_idx] = validate_stop - validate_start;

    if (!validation_passed_one) {
      validation_passed = 0;
      if (rank == 0) fprintf(stderr, "Validation failed for this BFS root; skipping rest.\n");
      break;
    }

    /* Calculate number of input edges visited. */
    {
      int64_t edge_visit_count = 0;
      size_t v_local;
      for (v_local = 0; v_local < g2.nlocalverts; ++v_local) {
        if (pred[v_local] != -1) {
          size_t ei, ei_end = g2.rowstarts[v_local + 1];
          for (ei = g2.rowstarts[v_local]; ei < ei_end; ++ei) {
            /* Test here is so that each input edge is counted exactly once, even
             * though non-self-loops are duplicated in the CSR data structure. */
            if (g2.column[ei] <= VERTEX_TO_GLOBAL(v_local)) {
              ++edge_visit_count;
            }
          }
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &edge_visit_count, 1, INT64_T_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
      edge_counts[bfs_root_idx] = (double)edge_visit_count;
    }
  }

  MPI_Free_mem(pred);
  free(bfs_roots);
  free_csr_graph(&g2);

  if (rank == 0) {
    if (!validation_passed) {
      fprintf(stdout, "No results printed for invalid run.\n");
    } else {
      int i;
      fprintf(stdout, "SCALE:                          %d\n", SCALE);
      fprintf(stdout, "edgefactor:                     %.2g\n", edgefactor);
      fprintf(stdout, "NBFS:                           %d\n", num_bfs_roots);
      fprintf(stdout, "num_mpi_processes:              %d\n", size);
      fprintf(stdout, "construction_time:              %g s\n", data_struct_time);
      double stats[s_LAST];
      get_statistics(bfs_times, num_bfs_roots, stats);
      fprintf(stdout, "min_time:                       %g s\n", stats[s_minimum]);
      fprintf(stdout, "firstquartile_time:             %g s\n", stats[s_firstquartile]);
      fprintf(stdout, "median_time:                    %g s\n", stats[s_median]);
      fprintf(stdout, "thirdquartile_time:             %g s\n", stats[s_thirdquartile]);
      fprintf(stdout, "max_time:                       %g s\n", stats[s_maximum]);
      fprintf(stdout, "mean_time:                      %g s\n", stats[s_mean]);
      fprintf(stdout, "stddev_time:                    %g\n", stats[s_std]);
      double* secs_per_edge = (double*)xmalloc(num_bfs_roots * sizeof(double));
      for (i = 0; i < num_bfs_roots; ++i) secs_per_edge[i] = bfs_times[i] / edge_counts[i];
      get_statistics(secs_per_edge, num_bfs_roots, stats);
      fprintf(stdout, "min_TEPS:                       %g TEPS\n", 1. / stats[s_maximum]);
      fprintf(stdout, "firstquartile_TEPS:             %g TEPS\n", 1. / stats[s_thirdquartile]);
      fprintf(stdout, "median_TEPS:                    %g TEPS\n", 1. / stats[s_median]);
      fprintf(stdout, "thirdquartile_TEPS:             %g TEPS\n", 1. / stats[s_firstquartile]);
      fprintf(stdout, "max_TEPS:                       %g TEPS\n", 1. / stats[s_minimum]);
      fprintf(stdout, "harmonic_mean_TEPS:             %g TEPS\n", 1. / stats[s_mean]);
    }
  }







  MPI_Finalize();
  return 0;
}

int print_graph_tuple(FILE* out, int64_t* result, int64_t n_local, int rank) {
  int i;
  for (i=0; i<n_local; ++i) {
      fprintf (out, "%d %d %d\n", (int)result[i*2]+1, (int)result[i*2 + 1]+1, 1);
  }
  return 1;
}

int print_graph_csr(FILE* out, size_t *rowptr, int64_t *colidx, int n_local) {
  int i, k;
  int64_t dst;
  for (i=0; i<n_local; ++i) {
    for (k = rowptr[i]; k < rowptr[i+1]; ++k) {
      dst = colidx[k]; 
      fprintf (out, "%d %d %d\n",(int)VERTEX_TO_GLOBAL_ALT(rank,i)+1,(int)dst+1,1);
    }
  }
  return 1;
}

int print_parts(FILE* out, int* parts, int n, int n_local) {
  int i;
  for (i=0; i<n; ++i) {
    int node = i;
    int node_owner = VERTEX_OWNER(node);
    int node_local_idx = VERTEX_LOCAL(node);
    int parts_idx = node_owner*n_local + node_local_idx;
    int v_part = parts[parts_idx];
    fprintf (out, "%d %d\n",node+1,v_part+1);
  }
  return 1;
}

void partition_graph_data_structure(csr_graph* const g) { 
// distribute the low-degree vertices as in redistribute.h
// distribute high-degree vertices using the owner function on the /destination/ vertex of each edge
// we will need a special CSR structure that has offset columns. 
  int n = g->nglobalverts;
  int n_local = g->nlocalverts;
  int offset = g->nlocalverts * rank; //!//Does this work?
  int nparts = size;
  int tot_nnz = 0;
  parts = (int*)malloc(n * sizeof(int));
  int *partsize_update = (int*)malloc(nparts * sizeof(int));
  int *old_partsize = (int*)malloc(nparts * sizeof(int));
  int *partsize = (int*)malloc(nparts * sizeof(int));
  int *partscore = (int*)malloc(nparts * sizeof(int));
  int *partcost = (int*)malloc(nparts * sizeof(int));
  int *vorder = (int*)malloc(n_local * sizeof(int)); 
  int oldpart;
  int emptyverts = 0;
  int randidx;
  int cutoff = 50;
  size_t *row;
  size_t vert;
  size_t k,  nnz_row, best_part;
  int64_t *colidx = g->column;
  int64_t node = 0;
  size_t *rowptr = g->rowstarts;
  float curr_score, best_score;
  genRandPerm(vorder, n_local);
  float gamma = 1.5;
  int i, s, l; //,j;

  g_perm = (int*)malloc(n * sizeof(int));
 
  if(1) { // Print graph
    char filename[256];
    sprintf(filename, "out_csr%02d.mat", rank);
    FILE *GraphFile;
    GraphFile = fopen(filename, "w");
    assert(GraphFile != NULL);
    print_graph_csr(GraphFile, rowptr, colidx, n_local);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  memset(parts, -1, n * sizeof(int));
  for (l=0; l<nparts; ++l) {
    partsize[l] = 0;
    old_partsize[l] = 0;
    partsize_update[l] = 0;
  }

  //int num_local_hi = 0;
  int localedges = (int)g->nlocaledges;
  MPI_Allreduce(&localedges, &tot_nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  float alpha = sqrt(2) * (tot_nnz/pow(n,gamma));
  if (rank==1) {fprintf(stdout,"n = %d, n_local = %d, local nnz = %d, total nnz = %d\n",n,n_local,localedges,tot_nnz);}
  int repeat_run;
  int run;
  for (repeat_run = 0; repeat_run < 35; repeat_run++) {
    for (run=0; run<nparts; run++) {
      if (rank == run) { //just partition one process after the other...
      //if (1) {
        for (i = 0; i < n_local; ++i) {
          for (l = 0; l < nparts; l++) {
            partscore[l] = 0;
            partcost[l] = 0;
          }
          vert = (size_t)vorder[i];
          //int global_vert_idx = VERTEX_TO_GLOBAL(rank, vert);
          int local_idx = offset + vert; //VERTEX_LOCAL(global_vert_idx);
          row = &rowptr[vert];
          nnz_row = *(row+1) - *row;
          oldpart = -1;
          if (nnz_row >= cutoff) { parts[local_idx] = nparts; }
          else if(nnz_row > 0) {
            for (k = *row; k < ((*row)+nnz_row); ++k) {
              node = colidx[k]; 
              int node_owner = VERTEX_OWNER(node);
              int node_local_idx = VERTEX_LOCAL(node);
              int parts_idx = node_owner*g->nlocalverts + node_local_idx;
              int node_part = parts[parts_idx]; /////
              if (node_part >= 0 && node_part < nparts) { 
                partscore[parts[parts_idx]]++; 
              }
            }
            for (s = 0; s < nparts; ++s) { 
              float dc = calc_dc(alpha,gamma,partsize[s]);
              int normscore = partscore[s] - (int)nnz_row;
              partcost[s] = normscore - dc; 
            }
            best_part = nparts-1;
            best_score = partcost[nparts-1];
            for (s = nparts-2; s >= 0; --s) { 
              curr_score = partcost[s]; 
              if (curr_score > best_score) {
                best_score = curr_score;  best_part = s;
              }
            }
            oldpart = parts[local_idx];
            parts[local_idx] = best_part;
            partsize[best_part]++;
            if (oldpart >= 0 && oldpart < nparts) { partsize[oldpart]--; }
          } else { // empty vertex, assign randomly
            emptyverts++;
            randidx = irand(nparts);
            oldpart = parts[local_idx];
            parts[local_idx] = randidx;
            partsize[randidx]++;
            if (oldpart >= 0 && oldpart < nparts) { partsize[oldpart]--; }
          }
          /*if (i % 32 == 0) { //(isPowerOfTwo(i)) {
            MPI_Allgather(parts+offset, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
            for (l=0; l<nparts; ++l) { 
              partsize_update[l] = partsize[l] - old_partsize[l];
            }            
            MPI_Allreduce(MPI_IN_PLACE, partsize_update, nparts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            for (l=0; l<nparts; ++l) { 
              old_partsize[l] += partsize_update[l]; 
              partsize[l] = old_partsize[l];
            }
          }*/
        }
      }
      MPI_Allgather(parts+offset, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
      for (l=0; l<nparts; ++l) { 
        partsize_update[l] = partsize[l] - old_partsize[l];
        //fprintf(stdout,"partsize[%d] on rank %d is %d. partsizeupdate is %d. oldsize is %d\n", l, rank, partsize[l], partsize_update[l], old_partsize[l]);
      }
      MPI_Allreduce(MPI_IN_PLACE, partsize_update, nparts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      for (l=0; l<nparts; ++l) { 
        old_partsize[l] += partsize_update[l]; 
        partsize[l] = old_partsize[l];
      }
    }
    int cutedges = mpi_compute_cut(rowptr, colidx, parts, nparts, n_local, offset, cutoff, g);
  }



  // Okay this is absurd but for now, randomly assign the large vertices to a permutation
  // To prevent overly large partitions
  for (i=offset; i<offset+n_local; ++i) {
    if (parts[i] == nparts) {
      parts[i] = irand(nparts);
    }
  }
  MPI_Allgather(parts+offset, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {  // Print Parts
    FILE *PartFile;
    PartFile = fopen("parts.mat", "w");
    assert(PartFile != NULL);
    print_parts(PartFile, parts, n, n_local);
  }

  // Sanity Checks
  for (i=0; i < g->nglobalverts; ++i) {
    assert(parts[i]>=0);
    assert(parts[i]<nparts);
  }
}

// Random permutation generator. Move to another file.
int* genRandPerm(int* orderList, int size) {
  assert(orderList);
  srand(time(NULL));
  // Generate 'identity' permutation
  int i;
  for (i = 0; i < size; i++) { orderList[i] = i; }
  shuffle_int(orderList, size);
  return orderList;
}

void shuffle_int(int *list, int len) {
  int j;
  int tmp;
  while(len) {
      j = irand(len);
      if (j != len - 1) {
        tmp = list[j];
        list[j] = list[len - 1];
        list[len - 1] = tmp;
      }
    len--;
  }
}

int irand(int n) {
  int r, rand_max = RAND_MAX - (RAND_MAX % n);
  while ((r = rand()) >= rand_max);
  return r / (rand_max / n);
}

float calc_dc(float alpha, float gamma, int len) {
  return (alpha*pow(len+0.5,gamma)) - (alpha*pow(len,gamma));
}

int mpi_compute_cut(size_t *rowptr, int64_t *colidx, int *parts, int nparts, int n_local, int offset, int cutoff, csr_graph* const g) {
  size_t vert;
  int nnz_row;
  int v_part;
  int cutedges = 0;
  int mytotedges = 0;
  int mytotlodegedges = 0;
  size_t *row;
  int i;
  size_t k;
  int emptyparts = 0;
  mytotedges = rowptr[n_local];
  for (i = 0; i < n_local; i++) {
    vert = i;
    row = &rowptr[vert];
    nnz_row = (int)(*(row+1) - *(row)); //nnz in row
    if (nnz_row < cutoff) { 
      v_part = parts[vert+offset];
      if (v_part == -1) {
        v_part = 0;
        emptyparts++;
      }
      // count edges to other partitions
      for (k = *row; k < ((*row)+nnz_row); ++k) {
        int node = colidx[k];
        int node_owner = VERTEX_OWNER(node);
        int node_local_idx = VERTEX_LOCAL(node);
        int parts_idx = node_owner*g->nlocalverts + node_local_idx;
        if (parts[parts_idx] < nparts) { mytotlodegedges++; } //count low degree edges
        if (parts[parts_idx] != v_part && parts[parts_idx] < nparts) { cutedges++; } //count low degree cut edges
      }
    }
  }
  int tot_cutedges;
  int tot_lodegedges;
  MPI_Allreduce(&cutedges, &tot_cutedges, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mytotlodegedges, &tot_lodegedges, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  //fprintf(stdout,"offset: %d emptyparts = %d cutedges = %d totcutedges = %d tot edges=%d mylodegedges=%d totlodegedges=%d\n",offset, emptyparts,cutedges,tot_cutedges,mytotedges,mytotlodegedges,tot_lodegedges);
  if (rank == 0) {   fprintf(stdout,"total cutedges = %d, percentage of total:%f \n", tot_cutedges, (float)tot_cutedges/tot_lodegedges); }
  return tot_cutedges;
}

static int compare_doubles(const void* a, const void* b) {
  double aa = *(const double*)a;
  double bb = *(const double*)b;
  return (aa < bb) ? -1 : (aa == bb) ? 0 : 1;
}

static void get_statistics(const double x[], int n, double r[s_LAST]) {
  double temp;
  int i;
  /* Compute mean. */
  temp = 0;
  for (i = 0; i < n; ++i) temp += x[i];
  temp /= n;
  r[s_mean] = temp;
  /* Compute std. dev. */
  temp = 0;
  for (i = 0; i < n; ++i) temp += (x[i] - r[s_mean]) * (x[i] - r[s_mean]);
  temp /= n - 1;
  r[s_std] = sqrt(temp);
  /* Sort x. */
  double* xx = (double*)xmalloc(n * sizeof(double));
  memcpy(xx, x, n * sizeof(double));
  qsort(xx, n, sizeof(double), compare_doubles);
  /* Get order statistics. */
  r[s_minimum] = xx[0];
  r[s_firstquartile] = (xx[(n - 1) / 4] + xx[n / 4]) * .5;
  r[s_median] = (xx[(n - 1) / 2] + xx[n / 2]) * .5;
  r[s_thirdquartile] = (xx[n - 1 - (n - 1) / 4] + xx[n - 1 - n / 4]) * .5;
  r[s_maximum] = xx[n - 1];
  /* Clean up. */
  free(xx);
}

