/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*  Adapted from code by:                                                  */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include "oned_csr.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#define isPowerOfTwo(x) ((x != 0) && !(x & (x - 1)))
#define TO_NEW_IDX(x) g_perm[x]
#define NEW_PART_OF_IDX(x) parts[TO_NEW_IDX(x)]

static oned_csr_graph g;

static int* g_perm;
static int* parts;
static int num_total_lo_deg_nodes = 0;
static int num_local_lo_deg_nodes = 0;
static int num_total_hi_deg_nodes = 0;
static int num_local_hi_deg_nodes = 0;
static int local_lo_deg_offset = 0;

static int64_t* g_oldq;
static int64_t* g_newq;
static unsigned long* g_visited;
static const int coalescing_size = 256;
static int64_t* g_outgoing;
static size_t* g_outgoing_counts /* 2x actual count */;
static MPI_Request* g_outgoing_reqs;
static int* g_outgoing_reqs_active;
static int64_t* g_recvbuf;

// Randperm prototypes. Move to another file
int* genRandPerm(int *orderList, int size);
void shuffle_int(int *list, int len);
int irand(int n);
float calc_dc(float alpha, float gamma, int len);
int mpi_compute_cut(size_t *rowptr, int64_t *colidx, int *parts, int nparts, int n_local, int offset, int cutoff);
void quickSort(int64_t a[], int l, int r);
int64_t partition(int64_t a[], int l, int r);
int remove_duplicates(int64_t a[], int n);
int print_parts(FILE* out, int* parts, int n);
int print_graph(FILE* out, size_t *rowptr, int64_t *colidx, int n_local, int offset);

int print_part_graph(FILE* out, size_t *rowptr, int64_t *colidx, int offset, size_t *hrowptr, int64_t *hcolidx);


  /* Predefined entities you can use in your BFS (from common.h and oned_csr.h):
   *   + rank: global variable containing MPI rank
   *   + size: global variable containing MPI size
   *   + DIV_SIZE: single-parameter macro that divides by size (using a shift
   *     when properly set up)
   *   + MOD_SIZE: single-parameter macro that reduces modulo size (using a
   *     mask when properly set up)
   *   + VERTEX_OWNER: single-parameter macro returning the owner of a global
   *     vertex number
   *   + VERTEX_LOCAL: single-parameter macro returning the local offset of a
   *     global vertex number
   *   + VERTEX_TO_GLOBAL: single-parameter macro converting a local vertex
   *     offset to a global number
   *   + g.nlocalverts: number of vertices stored on the local rank
   *   + g.nglobalverts: total number of vertices in the graph
   *   + g.nlocaledges: number of graph edges stored locally
   *   + g.rowstarts, g.column: zero-based compressed sparse row data
   *     structure for the local part of the graph
   *
   * All macros documented above evaluate their arguments exactly once.
   *
   * The graph is stored using a 1-D, cyclic distribution: all edges incident
   * to vertex v are stored on rank (v % size) (aka VERTEX_OWNER(v)).  Edges
   * that are not self-loops are stored twice, once for each endpoint;
   * duplicates edges are kept.  The neighbors of vertex v can be obtained on
   * rank VERTEX_OWNER(v); they are stored in elements
   * {g.rowstarts[VERTEX_LOCAL(v)] ... g.rowstarts[VERTEX_LOCAL(v) + 1] - 1}
   * (inclusive) of g.column.
   *
   * Upon exit, your BFS must have filled in:
   *   + pred (an array of size g.nlocalverts):
   *     - The predecessor of vertex v in the BFS tree should go into
   *       pred[VERTEX_LOCAL(v)] on rank VERTEX_OWNER(v)
   *     - The predecessor of root is root
   *     - The predecessor of any unreachable vertex is -1
   *
   * The validator will check this for correctness. */

void make_graph_data_structure(const tuple_graph* const tg) {
  convert_graph_to_oned_csr(tg, &g);
  const size_t nlocalverts = g.nlocalverts;
  g_oldq = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t));
  g_newq = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t));
  const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;
  int64_t visited_size = (nlocalverts + ulong_bits - 1) / ulong_bits;
  g_visited = (unsigned long*)xmalloc(visited_size * sizeof(unsigned long));
  g_outgoing = (int64_t*)xMPI_Alloc_mem(coalescing_size * size * 2 * sizeof(int64_t));
  g_outgoing_counts = (size_t*)xmalloc(size * sizeof(size_t)) /* 2x actual count */;
  g_outgoing_reqs = (MPI_Request*)xmalloc(size * sizeof(MPI_Request));
  g_outgoing_reqs_active = (int*)xmalloc(size * sizeof(int));
  g_recvbuf = (int64_t*)xMPI_Alloc_mem(coalescing_size * 2 * sizeof(int64_t));
}

void partition_graph_data_structure() { 
// distribute the low-degree vertices as in redistribute.h
// distribute high-degree vertices using the owner function on the /destination/ vertex of each edge
// we will need a special CSR structure that has offset columns. 
  int n = g.nglobalverts;
  int n_local = g.nlocalverts;
  int offset = g.nlocalverts * rank; //!//Does this work?
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
  int64_t *colidx = g.column;
  int64_t node = 0;
  size_t *rowptr = g.rowstarts;
  float curr_score, best_score;
  genRandPerm(vorder, n_local);
  float gamma = 1.5;
  int i, j, s, l;

  g_perm = (int*)malloc(n * sizeof(int));
 
  if(1) { // Print graph
    char filename[256];
    sprintf(filename, "file%02d.mat", rank);
    FILE *GraphFile;
    GraphFile = fopen(filename, "w");
    assert(GraphFile != NULL);
    print_graph(GraphFile, rowptr, colidx, n_local, offset);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  memset(parts, -1, n * sizeof(int));
  for (l=0; l<nparts; ++l) {
    partsize[l] = 0;
    old_partsize[l] = 0;
    partsize_update[l] = 0;
  }

  int num_local_hi = 0;
  int localedge = (int)g.nlocaledges;
  MPI_Allreduce(&localedge, &tot_nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  float alpha = sqrt(2) * (tot_nnz/pow(n,gamma));
  if (rank==1) {fprintf(stdout,"n = %d, n_local = %d, local nnz=%zu, total nnz=%d\n",n,n_local,g.nlocaledges,tot_nnz);}
  int repeat_run;
  int run;
  for (repeat_run = 0; repeat_run < 25; repeat_run++) {
    for (run=0; run<nparts; run++) {
      if (rank == run) { //just partition one process after the other...
      //if (1) {
        for (i = 0; i < n_local; ++i) {
          for (l = 0; l<nparts; l++) {
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
              int parts_idx = node_owner*g.nlocalverts + node_local_idx;
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
    int cutedges = mpi_compute_cut(rowptr, colidx, parts, nparts, n_local, offset, cutoff);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {  // Print Parts
    FILE *PartFile;
    PartFile = fopen("parts.mat", "w");
    assert(PartFile != NULL);
    print_parts(PartFile, parts, n);
  }

  int* edge_counts_per_owner = (int*)xmalloc(size * sizeof(int)); memset(edge_counts_per_owner, 0, size * sizeof(int));
  int* edge_counts_per_sender = (int*)xmalloc(size * sizeof(int)); memset(edge_counts_per_sender, 0, size * sizeof(int));
  int* node_counts_per_owner = (int*)xmalloc(size * sizeof(int)); memset(node_counts_per_owner, 0, size * sizeof(int));
  int* node_counts_per_sender = (int*)xmalloc(size * sizeof(int)); memset(node_counts_per_sender, 0, size * sizeof(int));

  int* perm_counts = (int*)xmalloc((size + 1) * sizeof(int)); memset(perm_counts, 0, (size+1) * sizeof(int));
  int* perm_displs = (int*)xmalloc((size + 2) * sizeof(int));

  int* edge_displs_per_owner = (int*)xmalloc((size + 1) * sizeof(int)); 
  int* edge_displs_per_sender = (int*)xmalloc((size + 1) * sizeof(int)); 
  int* node_displs_per_owner = (int*)xmalloc((size + 1) * sizeof(int)); 
  int* node_displs_per_sender = (int*)xmalloc((size + 1) * sizeof(int)); 

  //Count node counts to send to each process
  for (i = 0; i < n_local; ++i) {
    int part = parts[offset + i];
    if (part>=0 && part<nparts) { //!assert that this is the case
      node_counts_per_sender[part]++;
      edge_counts_per_sender[part] += g.rowstarts[i+1]-g.rowstarts[i]; //!make sure this is valid
    }
  }

/*
  if(rank == 0) {
      fprintf(stdout, "Parts on rank %d: ",rank);
  for (i = 0; i<n; i++) {
    fprintf(stdout, "%d ",parts[i]);
  }
  fprintf(stdout, "\n");
  }*/
/*
  fprintf(stdout, "Parts on rank %d: ",rank);
  for (i = 0; i<n; i++) {
    if (i == offset) {fprintf(stdout,"(");} else {fprintf(stdout," ");}
    if (parts[i] == rank) {fprintf(stdout, "+"); }
    else {fprintf(stdout, "-");}
    if (i == offset + n_local - 1) {fprintf(stdout,")");} else {fprintf(stdout," ");}
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "Nodectspsender on rank %d: ",rank);
  for (i = 0; i<size; ++i) {
    fprintf(stdout, "%d ",node_counts_per_sender[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "Edgectspsender on rank %d: ",rank);
  for (i = 0; i<size; ++i) {
    fprintf(stdout, "%d ",edge_counts_per_sender[i]);
  }
  fprintf(stdout, "\n");
*/


  for (i = 0; i < n; ++i) {
    int part = parts[i];
    if (part>=0 && part <= nparts) { perm_counts[part]++; }
  }

  perm_displs[0] = 0;
  for (i = 0; i <= size; ++i) {
    perm_displs[i+1] = perm_displs[i] + perm_counts[i];
  }

  local_lo_deg_offset = perm_displs[rank];

/*
  if(rank==0) {
  fprintf(stdout, "perm_counts: ",rank);
  for (i = 0; i<=size; ++i) {
    fprintf(stdout, "%d ",perm_counts[i]);
  }
  fprintf(stdout, "\n");
  }

  if(rank==0) {
  fprintf(stdout, "perm_displs: ",rank);
  for (i = 0; i<=size+1; ++i) {
    fprintf(stdout, "%d ",perm_displs[i]);
  }
  fprintf(stdout, "\n");
  }

*/
  size_t* row_offset = (size_t*)xmalloc(size * sizeof(size_t));

  int** perm_writers = (int**)xmalloc((nparts+1) * sizeof(int*));
  int* row_writers_idx = (int*)xmalloc(nparts * sizeof(int));

  size_t** row_writers = (size_t**)xmalloc(nparts * sizeof(size_t*));
  int64_t** col_writers = (int64_t**)xmalloc(nparts * sizeof(int64_t*));

  /* int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm) */

  MPI_Alltoall(edge_counts_per_sender, 1, MPI_INT, edge_counts_per_owner, 1, MPI_INT, MPI_COMM_WORLD);
  edge_displs_per_owner[0] = 0;
  edge_displs_per_sender[0] = 0;
  for (i = 0; i < size; ++i) {
    edge_displs_per_owner[i + 1] = edge_displs_per_owner[i] + edge_counts_per_owner[i];
    edge_displs_per_sender[i + 1] = edge_displs_per_sender[i] + edge_counts_per_sender[i];
  }

  MPI_Alltoall(node_counts_per_sender, 1, MPI_INT, node_counts_per_owner, 1, MPI_INT, MPI_COMM_WORLD);
  node_displs_per_owner[0] = 0;
  node_displs_per_sender[0] = 0;
  for (i = 0; i < size; ++i) {
    node_displs_per_owner[i + 1] = node_displs_per_owner[i] + node_counts_per_owner[i];
    node_displs_per_sender[i + 1] = node_displs_per_sender[i] + node_counts_per_sender[i];
  }
/*
  fprintf(stdout, "\n");
  fprintf(stdout, "Nodectspowner on rank %d: ",rank);
  for (i = 0; i<size; ++i) {
    fprintf(stdout, "%d ",node_counts_per_owner[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "Edgectspowner on rank %d: ",rank);
  for (i = 0; i<size; ++i) {
    fprintf(stdout, "%d ",edge_counts_per_owner[i]);
  }
  fprintf(stdout, "\n");
*/

  int row_send_len = node_displs_per_sender[size];
  int col_send_len = edge_displs_per_sender[size];
  size_t* row_sendbuffer = (size_t*)xmalloc(row_send_len * sizeof(size_t));
  int64_t* col_sendbuffer = (int64_t*)xmalloc(col_send_len * sizeof(int64_t));

  for (l = 0; l<nparts; ++l) { 
    row_writers[l] = &row_sendbuffer[node_displs_per_sender[l]]; //!fill in displs first, make sure strided right
    col_writers[l] = &col_sendbuffer[edge_displs_per_sender[l]]; 
    perm_writers[l] = &g_perm[perm_displs[l]];
    row_writers_idx[l] = node_displs_per_sender[l];
    row_offset[l] = 0;
  }

  perm_writers[nparts] = &g_perm[perm_displs[nparts]];
/*
  fprintf(stdout, "\n");
  fprintf(stdout, "row_writer on rank %d: ",rank);
  for (i = 0; i<size; ++i) {
    fprintf(stdout, "%d/%d ", row_writers_idx[i], row_send_len);
  }
  fprintf(stdout, "\n");*/

  // global permutation order
  for (i=0; i<n; ++i) {
    int part = parts[i];
    if (part>=0 && part<=nparts) { *(perm_writers[part]++) = i; }
    if (part == nparts) { num_total_hi_deg_nodes++; }
  }

  // permute colidxs
  for (i=0; i < n_local; i++) {
    size_t degree = g.rowstarts[i+1] - g.rowstarts[i];
    for (j=0; j<degree; ++j) {
      g.column[g.rowstarts[i] + j] = TO_NEW_IDX(g.column[g.rowstarts[i] + j]);
    }
  }


  /*fprintf(stdout, "\n");
  fprintf(stdout, "gperm on %d: ",rank);
  for (i = 0; i<n; ++i) {
    fprintf(stdout, "%d ",g_perm[i]);
  }
  fprintf(stdout, "\n");*/

  for (i=0; i < n_local; i++) {
    int part = parts[offset + i]; //is this the right mapping?
    if (part == nparts) { num_local_hi++; }
    else if (part>=0 && part<nparts) {
      size_t degree = g.rowstarts[i+1] - g.rowstarts[i];
      //fprintf(stdout, "%lu/%lu ",g.rowstarts[i], degree);
      *(row_writers[part]++) = degree;
      //Advance column buffer, Translate 
      for (j=0; j<degree; ++j) {
        int64_t targ = g.column[g.rowstarts[i] + j];
        //targ = TO_NEW_IDX(targ);
        *(col_writers[part]++) = targ; 
      }
    }
  }

  num_local_hi_deg_nodes = num_local_hi;
  //num_total_hi_deg_nodes already counted

  num_local_lo_deg_nodes = g.nlocalverts - num_local_hi_deg_nodes;
  num_total_lo_deg_nodes = g.nglobalverts - num_total_hi_deg_nodes;

  /*
  fprintf(stdout, "\n");
  fprintf(stdout, "Nodectspsender on rank %d: ",rank);
  for (i = 0; i<size; ++i) {
    fprintf(stdout, "%d ",node_counts_per_sender[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "rank: %d...   ", rank);
  for (i = 0; i<row_send_len; i++) {
    fprintf(stdout, "%lu ", row_sendbuffer[i]);
  }
  fprintf(stdout,"\n");*/
  // Now we have the send buffers ready. Do an Alltoall to get the 'per_owner' values.
//////

  //find sizes of recvbuff and sendbuff
  int row_recv_len = node_displs_per_owner[size] + 1;
  int col_recv_len = edge_displs_per_owner[size];
  /*for (l = 0; l<nparts; ++l) {
    row_recv_len += node_counts_per_owner[l];
    col_recv_len += edge_counts_per_owner[l];
  }*/

  size_t* row_recvbuffer = (size_t*)xmalloc(row_recv_len * sizeof(size_t));
  int64_t* col_recvbuffer = (int64_t*)xmalloc(col_recv_len * sizeof(int64_t));

  // Then do an all_to_all_v
  /*  int MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
                  const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
                  const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
                  MPI_Comm comm)    */


  MPI_Alltoallv(col_sendbuffer, edge_counts_per_sender, edge_displs_per_sender, MPI_UINT64_T, 
                col_recvbuffer, edge_counts_per_owner, edge_displs_per_owner, MPI_UINT64_T,
                MPI_COMM_WORLD); 
  MPI_Free_mem(col_sendbuffer);

  MPI_Alltoallv(row_sendbuffer, node_counts_per_sender, node_displs_per_sender, MPI_INT64_T, 
                row_recvbuffer, node_counts_per_owner, node_displs_per_owner, MPI_INT64_T,
                MPI_COMM_WORLD); 
  MPI_Free_mem(row_sendbuffer);
//#if 0
  // Then do a section-based scan of the new row buffers to update the nonzero values
  size_t idx_so_far = 0;
  size_t degree;
  for (i=0; i<row_recv_len; ++i) {
    degree = row_recvbuffer[i];
    row_recvbuffer[i] = idx_so_far;
    idx_so_far += degree;
  }
/*
    fprintf(stdout, "\n");
  fprintf(stdout, "rank: %d len: %d...   ", rank, row_recv_len);
  for (i = 0; i<row_recv_len; i++) {
    fprintf(stdout, "(%d)%lu ", i,row_recvbuffer[i]);
  }
  fprintf(stdout,"\n");
*/
  int* hedge_counts_per_owner = (int*)xmalloc((size+1) * sizeof(int)); 
  memset(hedge_counts_per_owner, 0, (size+1) * sizeof(int));
  int* hedge_counts_per_sender = (int*)xmalloc((size+1) * sizeof(int)); 
  memset(hedge_counts_per_sender, 0, (size+1) * sizeof(int));
  int* hnode_counts_per_owner = (int*)xmalloc((size+1) * sizeof(int)); 
  memset(hnode_counts_per_owner, 0, (size+1) * sizeof(int));
  int* hnode_counts_per_sender = (int*)xmalloc((size+1) * sizeof(int));
  for (i=0; i<=size; ++i) { hnode_counts_per_sender[i] = num_local_hi; } 

  int* hedge_displs_per_owner = (int*)xmalloc((size + 2) * sizeof(int)); 
  int* hedge_displs_per_sender = (int*)xmalloc((size + 2) * sizeof(int)); 
  int* hnode_displs_per_owner = (int*)xmalloc((size + 2) * sizeof(int)); 
  int* hnode_displs_per_sender = (int*)xmalloc((size + 2) * sizeof(int)); 

  //Compute colidx sizes
  for (i=0; i < n_local; i++) {
    int part = parts[offset + i]; //is this the right mapping?
    if (part==nparts) {
      size_t degree = g.rowstarts[i+1] - g.rowstarts[i];
      //*(row_writers[part]++) = degree;
      for (j=0; j<degree; ++j) {
        int64_t targ = g.column[g.rowstarts[i] + j];
        int targ_part = parts[targ];
        hedge_counts_per_sender[targ_part]++;
      }
    }
  }

  MPI_Alltoall(hedge_counts_per_sender, 1, MPI_INT, hedge_counts_per_owner, 1, MPI_INT, MPI_COMM_WORLD);
  hedge_displs_per_owner[0] = 0;
  hedge_displs_per_sender[0] = 0;
  for (i = 0; i <= size; ++i) {
    hedge_displs_per_owner[i + 1] = hedge_displs_per_owner[i] + hedge_counts_per_owner[i];
    hedge_displs_per_sender[i + 1] = hedge_displs_per_sender[i] + hedge_counts_per_sender[i];
  }

  MPI_Alltoall(hnode_counts_per_sender, 1, MPI_INT, hnode_counts_per_owner, 1, MPI_INT, MPI_COMM_WORLD);
  hnode_displs_per_owner[0] = 0;
  hnode_displs_per_sender[0] = 0;
  for (i = 0; i <= size; ++i) {
    hnode_displs_per_owner[i + 1] = hnode_displs_per_owner[i] + hnode_counts_per_owner[i];
    hnode_displs_per_sender[i + 1] = hnode_displs_per_sender[i] + hnode_counts_per_sender[i];
  }
/*
  fprintf(stdout, "hndctspersender on rank %d: ",rank);
  for (i = 0; i<=size; ++i) {
    fprintf(stdout, "%d ",hnode_counts_per_sender[i]);
  }
  fprintf(stdout, "\n");

    fprintf(stdout, "hedgectspersender on rank %d: ",rank);
  for (i = 0; i<=size; ++i) {
    fprintf(stdout, "%d ",hedge_counts_per_sender[i]);
  }
  fprintf(stdout, "\n");

      fprintf(stdout, "hedgedisplspersender on rank %d: ",rank);
  for (i = 0; i<=size+1; ++i) {
    fprintf(stdout, "%d ",hedge_displs_per_sender[i]);
  }
  fprintf(stdout, "\n");
*/
  // Then we need to partition the high-deg vertices 
  int num_hi = num_local_hi; //this should be proc local not perm
  size_t* hrowstarts_send = (size_t*)malloc((size+1) * num_hi * sizeof(size_t)); 
  int64_t* hcol_send = (int64_t*)malloc((size+1) * hedge_displs_per_sender[size] * sizeof(int64_t)); 

  memset(hrowstarts_send, 0, (size+1) * num_hi * sizeof(size_t));
  size_t** hrow_writers = (size_t**)xmalloc((nparts + 1) * sizeof(size_t*));
  int64_t** hcol_writers = (int64_t**)xmalloc((nparts + 1) * sizeof(int64_t*));
  
  //fprintf(stdout, "rank: %d numhi: %d...   ", rank, num_hi);


  for (l = 0; l<=nparts; ++l) { 
    hrow_writers[l] = &hrowstarts_send[num_hi*l];
    hcol_writers[l] = &hcol_send[hedge_displs_per_sender[l]]; 
  }

  for (i=0; i < n_local; i++) {
    int part = parts[offset + i]; //is this the right mapping?
    if (part==nparts) {
      size_t degree = g.rowstarts[i+1] - g.rowstarts[i];
      //*(row_writers[part]++) = degree;
      for (j=0; j<degree; ++j) {
        int64_t targ = g.column[g.rowstarts[i] + j];
        int targ_part = parts[targ];
        (*(hrow_writers[targ_part]))++;
        *(hcol_writers[targ_part]++) = targ; 
      }
      for (j=0; j<nparts; j++) {
        hrow_writers[j]++;
      }
    }
  }
/*
  fprintf(stdout, "rank: %d...   ", rank);
  for (i=0; i<size * num_hi; ++i) {
    fprintf(stdout, "%lu ", hrowstarts_send[i]);
  }
  fprintf(stdout,"\n");

  fprintf(stdout, "cols... rank: %d...   ", rank);
  for (i=0; i<hedge_displs_per_sender[size]; ++i) {
    fprintf(stdout, "%d ", (int)hcol_send[i]);
  }
  fprintf(stdout,"\n");*/
  #if 0
#endif

  int hrow_recv_len = hnode_displs_per_owner[size] + 1;
  int hcol_recv_len = hedge_displs_per_owner[size];
  /*for (l = 0; l<nparts; ++l) {
    row_recv_len += node_counts_per_owner[l];
    col_recv_len += edge_counts_per_owner[l];
  }*/

  size_t* hrow_recvbuffer = (size_t*)xmalloc(hrow_recv_len * sizeof(size_t));
  int64_t* hcol_recvbuffer = (int64_t*)xmalloc(hcol_recv_len * sizeof(int64_t));

  MPI_Alltoallv(hcol_send, hedge_counts_per_sender, hedge_displs_per_sender, MPI_UINT64_T, 
                hcol_recvbuffer, hedge_counts_per_owner, hedge_displs_per_owner, MPI_UINT64_T,
                MPI_COMM_WORLD); 
  MPI_Free_mem(hcol_send);

  MPI_Alltoallv(hrowstarts_send, hnode_counts_per_sender, hnode_displs_per_sender, MPI_INT64_T, 
                hrow_recvbuffer, hnode_counts_per_owner, hnode_displs_per_owner, MPI_INT64_T,
                MPI_COMM_WORLD); 
  MPI_Free_mem(hrowstarts_send);

  //update row values
  idx_so_far = 0;
  for (i=0; i<hrow_recv_len; ++i) {
    degree = hrow_recvbuffer[i];
    hrow_recvbuffer[i] = idx_so_far;
    idx_so_far += degree;
    //fprintf(stdout, "%d ", (int)hrow_recvbuffer[i]);
  }

  num_local_lo_deg_nodes = row_recv_len - 1;
  num_local_hi_deg_nodes = hrow_recv_len - 1;
  fprintf(stdout, "rank %d locallo: %d localhi: %d\n", rank, num_local_lo_deg_nodes, num_local_hi_deg_nodes);


  if(0) { // Print graph
    char filename[256];
    sprintf(filename, "file%02d.mat", rank);
    FILE *GraphFile;
    GraphFile = fopen(filename, "w");
    assert(GraphFile != NULL);
    print_part_graph(GraphFile, row_recvbuffer, col_recvbuffer, local_lo_deg_offset, hrow_recvbuffer, hcol_recvbuffer);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  free(partscore);
  free(vorder);
  free(parts);
  free(partsize);
  free(old_partsize);
  free(partsize_update);
  free(partcost);

  free(edge_counts_per_owner);
  free(edge_counts_per_sender);
  free(node_counts_per_owner);
  free(node_counts_per_sender);
  free(edge_displs_per_owner);
  free(edge_displs_per_sender);
  free(node_displs_per_owner);
  free(node_displs_per_sender);
//#endif
}

int mpi_compute_cut(size_t *rowptr, int64_t *colidx, int *parts, int nparts, int n_local, int offset, int cutoff) {
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
        int parts_idx = node_owner*g.nlocalverts + node_local_idx;
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

void free_graph_data_structure(void) {
  free(g_oldq);
  free(g_newq);
  free(g_visited);
  MPI_Free_mem(g_outgoing);
  free(g_outgoing_counts);
  free(g_outgoing_reqs);
  free(g_outgoing_reqs_active);
  MPI_Free_mem(g_recvbuf);
  free_oned_csr_graph(&g);
}

int bfs_writes_depth_map(void) {
  return 0;
}

/* This version is the traditional level-synchronized BFS using two queues.  A
 * bitmap is used to indicate which vertices have been visited.  Messages are
 * sent and processed asynchronously throughout the code to hopefully overlap
 * communication with computation. */
void run_bfs(int64_t root, int64_t* pred) {
  const size_t nlocalverts = g.nlocalverts;

  /* Set up the queues. */
  int64_t* restrict oldq = g_oldq;
  int64_t* restrict newq = g_newq;
  size_t oldq_count = 0;
  size_t newq_count = 0;

  /* Set up the visited bitmap. */
  const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;
  int64_t visited_size = (nlocalverts + ulong_bits - 1) / ulong_bits;
  unsigned long* restrict visited = g_visited;
  memset(visited, 0, visited_size * sizeof(unsigned long));
#define SET_VISITED(v) do {visited[VERTEX_LOCAL((v)) / ulong_bits] |= (1UL << (VERTEX_LOCAL((v)) % ulong_bits));} while (0)
#define TEST_VISITED(v) ((visited[VERTEX_LOCAL((v)) / ulong_bits] & (1UL << (VERTEX_LOCAL((v)) % ulong_bits))) != 0)

  /* Set up buffers for message coalescing, MPI requests, etc. for
   * communication. */
  const int coalescing_size = 256;
  int64_t* restrict outgoing = g_outgoing;
  size_t* restrict outgoing_counts = g_outgoing_counts;
  MPI_Request* restrict outgoing_reqs = g_outgoing_reqs;
  int* restrict outgoing_reqs_active = g_outgoing_reqs_active;
  memset(outgoing_reqs_active, 0, size * sizeof(int));
  int64_t* restrict recvbuf = g_recvbuf;
  MPI_Request recvreq;
  int recvreq_active = 0;

  /* Termination counter for each level: this variable counts the number of
   * ranks that have said that they are done sending to me in the current
   * level.  This rank can stop listening for new messages when it reaches
   * size. */
  int num_ranks_done;

  /* Set all vertices to "not visited." */
  {size_t i; for (i = 0; i < nlocalverts; ++i) pred[i] = -1;}

  /* Mark the root and put it into the queue. */
  if (VERTEX_OWNER(root) == rank) {
    SET_VISITED(root);
    //!// no need to check degree because we are simply enqueueing a received vertex -- sender should have broadcasted hi-degs
    //!// but we need to make sure we can use the same queue when receiving both p2p and broadcasted data
    pred[VERTEX_LOCAL(root)] = root;
    oldq[oldq_count++] = root;
  }

#define CHECK_MPI_REQS \
  /* Check all MPI requests and handle any that have completed. */ \
  do { \
    /* Test for incoming vertices to put onto the queue. */ \
    while (recvreq_active) { \
      int flag; \
      MPI_Status st; \
      MPI_Test(&recvreq, &flag, &st); \
      if (flag) { \
        recvreq_active = 0; \
        int count; \
        MPI_Get_count(&st, MPI_INT64_T, &count); \
        /* count == 0 is a signal from a rank that it is done sending to me
         * (using MPI's non-overtaking rules to keep that signal after all
         * "real" messages. */ \
        if (count == 0) { \
          ++num_ranks_done; \
        } else { \
          int j; \
          for (j = 0; j < count; j += 2) { \
            int64_t tgt = recvbuf[j]; \
            int64_t src = recvbuf[j + 1]; \
            /* Process one incoming edge. */ \
            assert (VERTEX_OWNER(tgt) == rank); \
            if (!TEST_VISITED(tgt)) { \
              SET_VISITED(tgt); \
              pred[VERTEX_LOCAL(tgt)] = src; \
              newq[newq_count++] = tgt; \
            } \
          } \
        } \
        /* Restart the receive if more messages will be coming. */ \
        if (num_ranks_done < size) { \
          MPI_Irecv(recvbuf, coalescing_size * 2, MPI_INT64_T, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvreq); \
          recvreq_active = 1; \
        } \
      } else break; \
    } \
    /* Mark any sends that completed as inactive so their buffers can be
     * reused. */ \
    int c; \
    for (c = 0; c < size; ++c) { \
      if (outgoing_reqs_active[c]) { \
        int flag; \
        MPI_Test(&outgoing_reqs[c], &flag, MPI_STATUS_IGNORE); \
        if (flag) outgoing_reqs_active[c] = 0; \
      } \
    } \
  } while (0)

  while (1) {
    memset(outgoing_counts, 0, size * sizeof(size_t));
    num_ranks_done = 1; /* I never send to myself, so I'm always done */
    
    /* Start the initial receive. */
    if (num_ranks_done < size) {
      MPI_Irecv(recvbuf, coalescing_size * 2, MPI_INT64_T, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvreq);
      recvreq_active = 1;
    }

    /* Step through the current level's queue. */
    size_t i;
    for (i = 0; i < oldq_count; ++i) {
      CHECK_MPI_REQS;
      assert (VERTEX_OWNER(oldq[i]) == rank);
      assert (pred[VERTEX_LOCAL(oldq[i])] >= 0 && pred[VERTEX_LOCAL(oldq[i])] < g.nglobalverts);
      int64_t src = oldq[i];
      /* Iterate through its incident edges. */
      size_t j, j_end = g.rowstarts[VERTEX_LOCAL(oldq[i]) + 1];
      for (j = g.rowstarts[VERTEX_LOCAL(oldq[i])]; j < j_end; ++j) {
        int64_t tgt = g.column[j];
        int owner = VERTEX_OWNER(tgt);
        /* If the other endpoint is mine, update the visited map, predecessor
         * map, and next-level queue locally; otherwise, send the target and
         * the current vertex (its possible predecessor) to the target's owner.
         * */
        if (owner == rank) {
          if (!TEST_VISITED(tgt)) {
            SET_VISITED(tgt);
            pred[VERTEX_LOCAL(tgt)] = src;
            newq[newq_count++] = tgt;
          }
        } else {
          while (outgoing_reqs_active[owner]) CHECK_MPI_REQS; /* Wait for buffer to be available */
          size_t c = outgoing_counts[owner];
          outgoing[owner * coalescing_size * 2 + c] = tgt;
          outgoing[owner * coalescing_size * 2 + c + 1] = src;
          outgoing_counts[owner] += 2;
          if (outgoing_counts[owner] == coalescing_size * 2) {
            MPI_Isend(&outgoing[owner * coalescing_size * 2], coalescing_size * 2, MPI_INT64_T, owner, 0, MPI_COMM_WORLD, &outgoing_reqs[owner]);
            outgoing_reqs_active[owner] = 1;
            outgoing_counts[owner] = 0;
          }
        }
      }
    }
    /* Flush any coalescing buffers that still have messages. */
    int offset;
    for (offset = 1; offset < size; ++offset) {
      int dest = MOD_SIZE(rank + offset);
      if (outgoing_counts[dest] != 0) {
        while (outgoing_reqs_active[dest]) CHECK_MPI_REQS;
        MPI_Isend(&outgoing[dest * coalescing_size * 2], outgoing_counts[dest], MPI_INT64_T, dest, 0, MPI_COMM_WORLD, &outgoing_reqs[dest]);
        outgoing_reqs_active[dest] = 1;
        outgoing_counts[dest] = 0;
      }
      /* Wait until all sends to this destination are done. */
      while (outgoing_reqs_active[dest]) CHECK_MPI_REQS;
      /* Tell the destination that we are done sending to them. */
      MPI_Isend(&outgoing[dest * coalescing_size * 2], 0, MPI_INT64_T, dest, 0, MPI_COMM_WORLD, &outgoing_reqs[dest]); /* Signal no more sends */
      outgoing_reqs_active[dest] = 1;
      while (outgoing_reqs_active[dest]) CHECK_MPI_REQS;
    }
    /* Wait until everyone else is done (and thus couldn't send us any more
     * messages). */
    while (num_ranks_done < size) CHECK_MPI_REQS;

    /* Test globally if all queues are empty. */
    int64_t global_newq_count;
    MPI_Allreduce(&newq_count, &global_newq_count, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    /* Quit if they all are empty. */
    if (global_newq_count == 0) break;

    /* Swap old and new queues; clear new queue for next level. */
    {int64_t* temp = oldq; oldq = newq; newq = temp;}
    oldq_count = newq_count;
    newq_count = 0;
  }
#undef CHECK_MPI_REQS
}

void get_vertex_distribution_for_pred(size_t count, const int64_t* vertex_p, int* owner_p, size_t* local_p) {
  const int64_t* restrict vertex = vertex_p;
  int* restrict owner = owner_p;
  size_t* restrict local = local_p;
  ptrdiff_t i;
#pragma omp parallel for
  for (i = 0; i < (ptrdiff_t)count; ++i) {
    owner[i] = VERTEX_OWNER(vertex[i]);
    local[i] = VERTEX_LOCAL(vertex[i]);
  }
}

int64_t vertex_to_global_for_pred(int v_rank, size_t v_local) { return VERTEX_TO_GLOBAL(v_rank, v_local); }
size_t get_nlocalverts_for_pred(void) { return g.nlocalverts; }

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

void quickSort(int64_t a[], int l, int r)
{
  int64_t j;
  if(l < r) {
    j = partition(a, l, r);
    quickSort(a, l, j-1);
    quickSort(a, j+1, r);
  }
}

int64_t partition(int64_t a[], int l, int r) {
  int64_t pivot, i, j, t;
  pivot = a[l];
  i = l; j = r + 1;
  while(1) {
    do ++i; while(a[i] <= pivot && i <= r );
    do --j; while(a[j] > pivot );
    if( i >= j ) break;
    t = a[i]; a[i] = a[j]; a[j] = t;
  }
  t = a[l]; a[l] = a[j]; a[j] = t;
  return j;
}

int remove_duplicates(int64_t a[], int n) {
  int i=0; int j;
  if (n <= 1) return n;
  for (j = 1; j < n; j++) {            
    if (a[j] != a[i]) {                
      a[++i] = a[j];
    }
  }
  return i+1;
}

int print_parts(FILE* out, int* parts, int n) {
  int i;
  for (i=0; i<n; ++i) {
    int node = i;
    int node_owner = VERTEX_OWNER(node);
    int node_local_idx = VERTEX_LOCAL(node);
    int parts_idx = node_owner*g.nlocalverts + node_local_idx;
    int v_part = parts[parts_idx];
    fprintf (out, "%d %d\n",node+1,v_part+1);
  }
  return 1;
}

int print_graph(FILE* out, size_t *rowptr, int64_t *colidx, int n_local, int offset) {
  int i, k;
  int64_t dst;
  for (i=0; i<n_local; ++i) {
    for (k = rowptr[i]; k < rowptr[i+1]; ++k) {
      dst = colidx[k]; 
      fprintf (out, "%d %d %d\n",(int)VERTEX_TO_GLOBAL(rank,i)+1,(int)dst+1,1);
    }
  }
  return 1;
}

int print_part_graph(FILE* out, size_t *rowptr, int64_t *colidx, int offset, size_t *hrowptr, int64_t *hcolidx) {
  int i, k;
  int64_t dst;

  //print lo-deg
  for (i=0; i < num_local_lo_deg_nodes; ++i) {
    for (k = rowptr[i]; k < rowptr[i+1]; ++k) {
      dst = colidx[k]; 
      fprintf (out, "%d %d %d\n",offset+i+1,(int)dst+1,1);
    }
  }
  
  //print hi-deg
  for (i=0; i < num_local_hi_deg_nodes; ++i) {
    for (k = hrowptr[i]; k < hrowptr[i+1]; ++k) {
      dst = hcolidx[k]; 
      fprintf (out, "%d %d %d\n",num_total_lo_deg_nodes+i+1,(int)dst+1,1);
    }
  }

  return 1;
}
