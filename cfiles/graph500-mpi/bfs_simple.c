/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
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

int print_graph(FILE* out, size_t *rowptr, int64_t *colidx, int n_local, int offset);


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
  /*
    char name[14] = "mygraph1.mat";
    char name2[14] = "mygraph2.mat";

    char* targ;
    if (rank==1) { targ = name2;} else { targ = name; }
        FILE *GraphFile;
        GraphFile = fopen(targ, "w");
        assert(GraphFile != NULL);
        print_graph(GraphFile, g.rowstarts, g.column, g.nlocalverts, rank*g.nlocalverts);
    //  }
      MPI_Barrier(MPI_COMM_WORLD);  */
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

  double COMP_START_TIMES[1024];
  double COMP_END_TIMES[1024];
  double COMM_START_TIMES[1024];
  double COMM_END_TIMES[1024];
  int iter = 0;

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
    

    //!start time of computation here
    COMP_START_TIMES[iter] = MPI_Wtime();

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

    //! end computation measure here
    COMP_END_TIMES[iter] = MPI_Wtime();

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

    //!time this
    while (num_ranks_done < size) CHECK_MPI_REQS;

    //! end communication measure here
    COMM_END_TIMES[iter] = MPI_Wtime();

    /* Test globally if all queues are empty. */
    int64_t global_newq_count;
    MPI_Allreduce(&newq_count, &global_newq_count, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    /* Quit if they all are empty. */
    if (global_newq_count == 0) break;

    /* Swap old and new queues; clear new queue for next level. */
    {int64_t* temp = oldq; oldq = newq; newq = temp;}
    oldq_count = newq_count;
    newq_count = 0;
    iter++;
  }
#undef CHECK_MPI_REQS



  //ostringstream pertimes;
  //pertimes << LOC_SPMV_TIMES[iterations] + LOC_MERGE_TIMES[iterations] + LOC_TRANS_TIMES[iterations] << ",";
  if (rank == 0) { 
    fprintf(stdout, "$\n");  
    fprintf(stdout, "%d\n", size); 
  }
  size_t  proc, i;
  double zerotime = COMP_START_TIMES[0];
  for (i = 0; i < iter; i++ ) {
    //COMP_START_TIMES[i] = COMP_START_TIMES[i] - zerotime;
    //COMP_END_TIMES[i] = COMP_END_TIMES[i] - zerotime;
    //COMM_END_TIMES[i] = COMM_END_TIMES[i] - zerotime;
  }
  for (proc = 0; proc < size; proc++ ) {
    if (rank == proc) { 
      for (i = 0; i<iter; i++) {
        fprintf (stdout, "%f, %f, %f, ", COMP_START_TIMES[i], COMP_END_TIMES[i], COMM_END_TIMES[i]);
      }
      fprintf(stdout, "\n"); 
    }
    usleep(200);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (rank == 0) {  
    //MPI_Barrier(MPI_COMM_WORLD);
    usleep(200);
    fprintf(stdout, "$\n"); 
  }
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

int64_t vertex_to_global_for_pred(int v_rank, size_t v_local) {
  return VERTEX_TO_GLOBAL(v_rank, v_local);
}

size_t get_nlocalverts_for_pred(void) {
  return g.nlocalverts;
}

int print_graph(FILE* out, size_t *rowptr, int64_t *colidx, int n_local, int offset) {
  int i, k;
  int64_t dst;

  for (i=0; i<n_local; ++i) {
    for (k = rowptr[i]; k < rowptr[i+1]; ++k) {
      dst = colidx[k]; 
      fprintf (out, "%d %d %d\n",offset+i+1,(int)dst+1,1);
    }
  }
  return 1;
}
