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

static oned_csr_graph g;
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
  int *parts = (int*)malloc(n * sizeof(int));
  int *partsize_update = (int*)malloc(nparts * sizeof(int));
  int *old_partsize = (int*)malloc(nparts * sizeof(int));
  int *partsize = (int*)malloc(nparts * sizeof(int));
  int *partscore = (int*)malloc(nparts * sizeof(int));
  int *partcost = (int*)malloc(nparts * sizeof(int));
  int *vorder = (int*)malloc(n_local * sizeof(int)); 
  int oldpart;
  int emptyverts = 0;
  int randidx;
  int cutoff = 100;
  size_t *row;
  size_t vert;
  size_t k,  nnz_row, best_part;
  int64_t *colidx = g.column;
  int64_t node = 0;
  size_t *rowptr = g.rowstarts;
  float curr_score, best_score;
  genRandPerm(vorder, n_local);
  float gamma = 1.5;
  int i, s, l;
 
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

  int localedge = (int)g.nlocaledges;
  MPI_Allreduce(&localedge, &tot_nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  float alpha = sqrt(2) * (tot_nnz/pow(n,gamma));
  if (rank==1) {fprintf(stdout,"n = %d, n_local = %d, local nnz=%zu, total nnz=%d\n",n,n_local,g.nlocaledges,tot_nnz);}
  int repeat_run;
  int run;
  for (repeat_run = 0; repeat_run < 16; repeat_run++) {
    for (run=0; run<nparts; run++) {
      if (rank == run) { //for now just partition one after the other...
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
          /*if (i % 16 == 0) {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allgather(parts+offset, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
            MPI_Allreduce(local_partsize, partsize, nparts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
  free(partscore);
  free(vorder);
  free(parts);
  free(partsize);
  free(old_partsize);
  free(partsize_update);
  free(partcost);
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
        if (parts[parts_idx] < nparts) { mytotlodegedges++; }
        if (parts[parts_idx] != v_part && parts[parts_idx] < nparts) { cutedges++; }
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
