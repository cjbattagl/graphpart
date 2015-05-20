
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

#include <time.h>
#include <math.h>

///////////////////
//#define TO_NEW_IDX(x) g_perm[x]
//#define NEW_PART_OF_IDX(x) parts[TO_NEW_IDX(x)]
///////////////////
static int* g_perm;
static int* parts;
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

static int num_hi_deg_verts;
static int64_t* hi_column;
static size_t* hi_rowstarts;
static int* hi_deg_ids;
static int local_hi_nnzs;

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

void make_graph_data_structure(const tuple_graph* const tg) {
  convert_graph_to_oned_csr(tg, &g);
  const size_t nlocalverts = g.nlocalverts;
}

void remake_graph_data_structure(const tuple_graph* const tg) {
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


int64_t get_permed_vertex(int64_t id) { 
  return ((g_perm[id] % g.nlocalverts) * size + floor(g_perm[id]/(g.nglobalverts/size))); 
}

void print_graph() {
    char filename[256];
    sprintf(filename, "out_pcsr%02d.mat", rank);
    FILE *GraphFile;
    GraphFile = fopen(filename, "w");
    assert(GraphFile != NULL);
    print_graph_csr(GraphFile, g.rowstarts, g.column, g.nlocalverts);
    MPI_Barrier(MPI_COMM_WORLD);
    fclose(GraphFile);
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

  //////Set up Hi-deg visited bitmap
  int64_t hi_visited_size = (num_hi_deg_verts + ulong_bits - 1) / ulong_bits; //number of longs needed to fit all bits
  unsigned long* hi_deg_visited;
  unsigned long* hi_deg_visited_new;
  unsigned long* hi_deg_visited_next;

  hi_deg_visited = (unsigned long*)xMPI_Alloc_mem(hi_visited_size * sizeof(unsigned long));
  hi_deg_visited_new = (unsigned long*)xMPI_Alloc_mem(hi_visited_size * sizeof(unsigned long));
  hi_deg_visited_next = (unsigned long*)xMPI_Alloc_mem(hi_visited_size * sizeof(unsigned long));

  memset(hi_deg_visited, 0, hi_visited_size * sizeof(unsigned long));
  memset(hi_deg_visited_new, 0, hi_visited_size * sizeof(unsigned long));
  memset(hi_deg_visited_next, 0, hi_visited_size * sizeof(unsigned long));

#define SET_HI_VISITED(v) do {hi_deg_visited_next[hi_deg_ids[v] / ulong_bits] |= (1UL << (hi_deg_ids[v] % ulong_bits));} while (0)
#define TEST_HI_VISITED(v) ((hi_deg_visited_new[hi_deg_ids[v] / ulong_bits] & (1UL << (hi_deg_ids[v] % ulong_bits))) != 0)
#define TEST_HI_NEXT_VISITED(v) ((hi_deg_visited_next[hi_deg_ids[v] / ulong_bits] & (1UL << (hi_deg_ids[v] % ulong_bits))) != 0)
#define TEST_HI_EVER_VISITED(v) ((hi_deg_visited[hi_deg_ids[v] / ulong_bits] & (1UL << (hi_deg_ids[v] % ulong_bits))) != 0)
#define IS_HI_DEG_VERTEX(v) (hi_deg_ids[tgt] >= 0)

  //MPI_Allreduce(MPI_IN_PLACE, hi_deg_visited, hi_visited_size, MPI_UNSIGNED_LONG, MPI_BOR, MPI_COMM_WORLD);


  /* Set up buffers for message coalescing, MPI requests, etc. for
   * communication. */
  const int coalescing_size = 256;
  int64_t* restrict outgoing = g_outgoing;
  size_t* restrict outgoing_counts = g_outgoing_counts;
  MPI_Request* restrict outgoing_reqs = g_outgoing_reqs;
  int* restrict outgoing_reqs_active = g_outgoing_reqs_active;
  assert(outgoing_reqs_active);
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

    size_t i,j;

    /* Step through incoming hi-degree vertices */
    for (i = 0; i < num_hi_deg_verts; ++i) {
      if (TEST_HI_VISITED(i)) {
        //we never reach here
        int64_t src = hi_deg_ids[i];

        for (j = hi_rowstarts[hi_deg_ids[i]]; j<hi_rowstarts[hi_deg_ids[i]+1]; ++j) {
        //assert(local_hi_nnzs >= hi_rowstarts[hi_deg_ids[i]+1]);
        int64_t tgt = hi_column[j];
        int owner = VERTEX_OWNER(tgt);
          if (IS_HI_DEG_VERTEX(tgt)) {
            if (!TEST_HI_EVER_VISITED(tgt)) {
              SET_HI_VISITED(tgt); //Hi degree, unvisited
              if (owner == rank) {
                if (!TEST_VISITED(tgt)) {
                  SET_VISITED(tgt);
                  pred[VERTEX_LOCAL(tgt)] = src;
                }
              }
            }
          }
          else { //Lo degree
            assert(owner==rank);
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
      }
    }
    /* Step through the current level's queue. */
    for (i = 0; i < oldq_count; ++i) {
      CHECK_MPI_REQS;
      assert (VERTEX_OWNER(oldq[i]) == rank);
      assert (pred[VERTEX_LOCAL(oldq[i])] >= 0 && pred[VERTEX_LOCAL(oldq[i])] < g.nglobalverts);
      int64_t src = oldq[i];
      /* Iterate through its incident edges. */
      size_t j, j_end = g.rowstarts[VERTEX_LOCAL(oldq[i]) + 1];


      int this_nnz = g.rowstarts[VERTEX_LOCAL(oldq[i]) + 1] - g.rowstarts[VERTEX_LOCAL(oldq[i])];
      if (this_nnz >= F_CUTOFF) { num_sends--; num_bcasts++; }

      // BROADCAST
      /*
      if (this_nnz >= F_CUTOFF) {
        num_processed+=this_nnz; 
        if (!TEST_HI_VISITED(tgt)) {
          SET_HI_VISITED(tgt);
          pred[VERTEX_LOCAL(tgt)] = src;
        }
      }*/
      for (j = g.rowstarts[VERTEX_LOCAL(oldq[i])]; j < j_end; ++j) {
        int64_t tgt = g.column[j];
        int owner = VERTEX_OWNER(tgt);

        if (IS_HI_DEG_VERTEX(tgt)) {
          //if (TEST_HI_EVER_VISITED(tgt)) { fprintf(stdout, " %d ", (int)tgt); }
          if (!TEST_HI_EVER_VISITED(tgt)) {
            SET_HI_VISITED(tgt);
            if (owner == rank) { pred[VERTEX_LOCAL(tgt)] = src; }
            num_processed++; 
          }
        }
        else {
          //UPDATE COMMUNICATION COUNTS
          //if (owner == rank) { }
          //else { num_sends++;
          //  if (this_nnz >= F_CUTOFF) { num_sends--; }
          //}
          num_processed++; 
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

    //HI BMAP |= NEW_HI BMAP
    for (i=0; i<hi_visited_size; ++i) { hi_deg_visited[i] |= hi_deg_visited_new[i]; }
    //MEMSET NEW_HI -> 0
    memset(hi_deg_visited_new, 0, hi_visited_size * sizeof(unsigned long));
    //BROADCAST NEXT_HI -> NEW_HI BMAPS
    MPI_Allreduce(hi_deg_visited_next, hi_deg_visited_new, hi_visited_size, MPI_UNSIGNED_LONG, MPI_BOR, MPI_COMM_WORLD);
    //MEMSET NEXT_HI -> 0
    memset(hi_deg_visited_next, 0, hi_visited_size * sizeof(unsigned long));

    //for (i=0; i<hi_visited_size; ++i) { hi_deg_visited_old[i] |= hi_deg_visited_new[i]; }

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

void distribute_hi_degrees() {
  int n = g.nglobalverts;
  int n_local = g.nlocalverts;
  int offset = g.nlocalverts * rank; //!//Does this work?
  int nparts = size;
  int cutoff = F_CUTOFF;
  size_t *row;
  //size_t vert = 0;
  size_t k,  nnz_row, best_part;
  int64_t *colidx = g.column;
  int64_t node = 0;
  size_t *rowptr = g.rowstarts;
  int i,j,l;

  int* hi_deg_total_sends = (int*)malloc(nparts * sizeof(int));
  int* hi_deg_total_recvs = (int*)malloc(nparts * sizeof(int));

  for (l=0; l<nparts; ++l) { hi_deg_total_sends[l] = 0; }

  MPI_Allreduce(MPI_IN_PLACE, &num_hi_deg_verts, nparts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  int num_local_hi_verts = 0;

  if (1 && rank==0) { fprintf(stdout,"n = %d, n_local = %d, num hi = %d\n",n,n_local,num_hi_deg_verts); }

//////////////////
#if 0
  const int ulong_bits = sizeof(unsigned long) * CHAR_BIT; //number of bits in ulong
  int64_t visited_size = (num_hi_deg_verts + ulong_bits - 1) / ulong_bits; //number of longs needed to fit all bits
  unsigned long* hi_deg_visited;
  hi_deg_visited = (unsigned long*)xMPI_Alloc_mem(visited_size * sizeof(unsigned long));
  memset(hi_deg_visited, 0, visited_size * sizeof(unsigned long));
#define SET_VISITED(v) do {hi_deg_visited[hi_deg_ids[v] / ulong_bits] |= (1UL << (hi_deg_ids[v] % ulong_bits));} while (0)
#define TEST_VISITED(v) ((hi_deg_visited[hi_deg_ids[v] / ulong_bits] & (1UL << (hi_deg_ids[v] % ulong_bits))) != 0)
  MPI_Allreduce(MPI_IN_PLACE, hi_deg_visited, visited_size, MPI_UNSIGNED_LONG, MPI_BOR, MPI_COMM_WORLD);
#endif
//////////////////

  hi_deg_ids = (int*)malloc(n * sizeof(int));
  memset(hi_deg_ids, -1, n * sizeof(int));

  // Count total sends so we can allocate basic data structures
  for (i = 0; i < n_local; ++i) {
    row = &rowptr[i];
    nnz_row = *(row+1) - *row;
    if (nnz_row >= cutoff) {
      num_local_hi_verts++;

      size_t local_idx = offset + i; 

      //fprintf(stdout, "%d ",local_idx);
      hi_deg_ids[local_idx] = num_local_hi_verts-1;
      // count total edges we are sending to each process
      for (k = *row; k < ((*row)+nnz_row); ++k) {
        node = colidx[k]; 
        int node_owner = VERTEX_OWNER(node);
        hi_deg_total_sends[node_owner]++;
      }
    }
  }

  MPI_Allreduce(hi_deg_total_sends, hi_deg_total_recvs, nparts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allgather(MPI_IN_PLACE, n_local, MPI_INT, hi_deg_ids, n_local, MPI_INT, MPI_COMM_WORLD);

  for (l=0; l<nparts; ++l) {
    //fprintf(stdout, "SEND %d: %d\n", rank, hi_deg_total_sends[l]);
    //fprintf(stdout, "RECV %d: %d\n", rank, hi_deg_total_recvs[l]);
  }

  for (l=0; l<n; ++l) {
    //fprintf(stdout, "%d ", hi_deg_ids[l]);
    //fprintf(stdout, "RECV %d: %d\n", rank, hi_deg_total_recvs[l]);
  }

  //size_t total_recvs = 0;
  size_t total_sends = 0;
  for (l=0; l<nparts; ++l) {
    //total_recvs += hi_deg_total_recvs[l];
    total_sends += hi_deg_total_sends[l];
  }

  if (VERBY) { fprintf(stdout,"total_sends = %d,  num_hi_deg_verts = %d, num hi = %d\n",(int)total_sends, num_hi_deg_verts, num_local_hi_verts); }

  hi_rowstarts = (size_t*)malloc((num_hi_deg_verts+1)*sizeof(size_t));
  int* sendcounts = (int *)malloc( size * sizeof(int) );
  int* recvcounts = (int *)malloc( size * sizeof(int) );
  int* rdispls = (int *)malloc( (size+1) * sizeof(int) );
  int* sdispls = (int *)malloc( (size+1) * sizeof(int) );

  int64_t* send_buffer = (int64_t*)malloc(total_sends * sizeof(int64_t));
  int64_t* send_buffer_writers = (int64_t*)malloc(nparts * sizeof(int64_t));

  int so_far = 0;
  for (l=0; l<nparts; ++l) {
    send_buffer_writers[l] = so_far;
    sdispls[l] = so_far;
    sendcounts[l] = hi_deg_total_sends[l];
    so_far += hi_deg_total_sends[l];
  }

  MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
  rdispls[0] = 0; 
  sdispls[0] = 0; 
  for (i = 0; i < size; ++i) { 
    rdispls[i + 1] = rdispls[i] + recvcounts[i]; 
    sdispls[i + 1] = sdispls[i] + sendcounts[i]; 
  }

  size_t* row_ptr_sends = ( size_t* )malloc(size * num_local_hi_verts * sizeof(size_t*));
  int* row_ptr_sends_map = (int*)malloc(size*sizeof(int*));
  for (i=0; i<size; i++) { 
    row_ptr_sends_map[i] = (num_local_hi_verts)*i; 
    //row_ptr_sends[row_ptr_sends_map[i]] = 0;
  }

  int* targ_proc_counts = (int*)malloc(size * sizeof(int*));
  memset(targ_proc_counts, 0, size * sizeof(int));
  size_t new_vert_id = 0;

  // Count number of incident edges to send to each processor, to compute rowstarts
  // Bin together edges to send to each processor, to compute col_idx
  for (i = 0; i < n_local; ++i) {
    row = &rowptr[i];
    nnz_row = *(row+1) - *row;
    if (nnz_row >= cutoff) {
      //int local_idx = offset + vert; 
      // count total edges we are sending to each process
      for (k = *row; k < ((*row)+nnz_row); ++k) {
        node = colidx[k]; 
        int node_owner = VERTEX_OWNER(node);
        send_buffer[send_buffer_writers[node_owner]] = node;
        send_buffer_writers[node_owner]++;
        targ_proc_counts[node_owner]++;
      }

      for (k = 0; k<size; ++k) {
        int vert_idx = row_ptr_sends_map[k] + new_vert_id;
        row_ptr_sends[vert_idx] = targ_proc_counts[k];
      }

      new_vert_id++;
      memset(targ_proc_counts, 0, size * sizeof(int));
    }
  }

  if(VERBY){
    for(j=0; j<size; j++){
      for (i=0; i<num_local_hi_verts; i++) {
        fprintf(stdout, " %d ", (int)row_ptr_sends[row_ptr_sends_map[j] + i]);
      }
      fprintf(stdout, "\n");
    }
  }

  if(VERBY){
    fprintf(stdout, "SEND  %d: ", rank);
    for (l=0; l<nparts; ++l) {
      fprintf(stdout, " %d ", sendcounts[l]);
    }
    fprintf(stdout, "\nRECV %d: ", rank);
    for (l=0; l<nparts; ++l) {
      fprintf(stdout, " %d ", recvcounts[l]);
    }
    fprintf(stdout, "\nSEND_DSPL  %d: ", rank);
    for (l=0; l<=nparts; ++l) {
      fprintf(stdout, " %d ", sdispls[l]);
    }
    fprintf(stdout, "\nRECV_DSPL %d: ", rank);
    for (l=0; l<=nparts; ++l) {
      fprintf(stdout, " %d ", rdispls[l]);
    }
    fprintf(stdout, "\n");
  }

  /////////Send columns
  hi_column = (int64_t*)malloc(rdispls[size] * sizeof(int64_t));
  local_hi_nnzs = rdispls[size];
  MPI_Alltoallv(send_buffer, sendcounts, sdispls, MPI_INT64_T, hi_column, recvcounts, rdispls, MPI_INT64_T, MPI_COMM_WORLD);
  
  /////////Send rowstarts
  int* r_sendcounts = (int *)malloc( size * sizeof(int) );
  int* r_recvcounts = (int *)malloc( size * sizeof(int) );
  int* r_rdispls = (int *)malloc( (size+1) * sizeof(int) );
  int* r_sdispls = (int *)malloc( (size+1) * sizeof(int) );
  for (i = 0; i < size; ++i) { r_sendcounts[i] = num_local_hi_verts; }
  MPI_Alltoall(r_sendcounts, 1, MPI_INT, r_recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
  r_rdispls[0] = 0; 
  r_sdispls[0] = 0; 
  for (i = 0; i < size; ++i) { 
    r_rdispls[i + 1] = r_rdispls[i] + r_recvcounts[i]; 
    r_sdispls[i + 1] = r_sdispls[i] + r_sendcounts[i]; 
  }
  MPI_Alltoallv(row_ptr_sends, r_sendcounts, r_sdispls, MPI_INT64_T, hi_rowstarts, r_recvcounts, r_rdispls, MPI_INT64_T, MPI_COMM_WORLD);
  //for (i=0; i<num_hi_deg_verts; ++i) { fprintf(stdout, "%d ", hi_rowstarts[i]); }
  //fprintf(stdout, "\n");
  size_t temp = hi_rowstarts[0];
  size_t sum_so_far = 0;
  hi_rowstarts[0] = 0;
  for (i=0; i<num_hi_deg_verts; i++) {
    sum_so_far += temp;
    temp = hi_rowstarts[i+1];
    hi_rowstarts[i+1] = sum_so_far;
  }

  for (i=1; i<nparts; ++i) {
    size_t hi_offset = r_rdispls[i];
    for (j=0; j<n_local; ++j) {
      if (hi_deg_ids[j+(n_local*i)] >= 0) {
        hi_deg_ids[j+(n_local*i)] += hi_offset;
      }
    }
  }

  for (l=0; l<n; ++l) {
    //if (hi_deg_ids[l] >= 0) {fprintf(stdout, "%d ", hi_deg_ids[l]);}
    //fprintf(stdout, "RECV %d: %d\n", rank, hi_deg_total_recvs[l]);
  }
  //for (i=0; i<=num_hi_deg_verts; ++i) { fprintf(stdout, "%d ", hi_rowstarts[i]); }
  //fprintf(stdout, "\n");
  if(rank==0) { printf("hi deg verts: %d, local hi nnzs: %d\n", (int)num_hi_deg_verts, local_hi_nnzs); }

}

void permute_tuple_graph(tuple_graph* tg) {
  if (tg->edgememory_size > 0) {
    int64_t i;
    int perm_local_size = g.nglobalverts;
    int* perms = (int*)malloc(g.nglobalverts * sizeof(int));
    for (i=0; i<g.nglobalverts; ++i) {
      int node = i;
      int node_owner = VERTEX_OWNER(node);
      int node_local_idx = VERTEX_LOCAL(node);
      int parts_idx = node_owner*g.nlocalverts + node_local_idx;
      int v_part = parts[parts_idx];
      perms[i] = v_part;
    }

    int64_t* perm_sizes = (int64_t*)xmalloc((size+1) * sizeof(int64_t));
    for (i=0; i<size; ++i) { perm_sizes[i] = 0; }
    for (i=0; i < g.nglobalverts; ++i) { perm_sizes[parts[i]]++; }
    int64_t* perm_displs = (int64_t*)xmalloc((size + 1) * sizeof(int64_t));
    perm_displs[0] = 0;
    for (i = 1; i < size + 1; ++i) perm_displs[i] = perm_displs[i - 1] + perm_sizes[i - 1];

    int64_t* local_vertex_perm = NULL;
    local_vertex_perm = (int64_t*)malloc(perm_local_size * sizeof(int64_t));

    int perm_id = 0;
    int64_t* perm_idxs = (int64_t*)xmalloc((size + 1) * sizeof(int64_t));
    perm_idxs[0] = 0;
    for (i = 1; i < size + 1; ++i) perm_idxs[i] = perm_displs[i];

    for (i=0; i < g.nglobalverts; ++i) {
      int part = perms[i];
      g_perm[perm_id] = perm_idxs[part];
      perm_id++;
      perm_idxs[part]++;
    }

    packed_edge* result = tg->edgememory;
    int64_t v0, v1;
    int64_t v0_new, v1_new;
    int src_proc, tgt_proc;
    int64_t partsize = (int)(g.nglobalverts / size);
    packed_edge* edge;
    for (i = 0; i < (int64_t)tg->edgememory_size; ++i) {
      //fprintf(stderr, "%d ", (int)tg->edgememory_size);
      edge = &result[i];
      v0 = get_v0_from_edge(edge);
      v1 = get_v1_from_edge(edge);
      //we have to do this because it mods by source vertex....
      src_proc = floor(g_perm[v0]/(g.nglobalverts/size));
      tgt_proc = floor(g_perm[v1]/(g.nglobalverts/size));
      v0_new = g_perm[v0];
      v1_new = g_perm[v1];
      v0_new = (v0_new % partsize) * size + src_proc;
      v1_new = (v1_new % partsize) * size + tgt_proc;
      //fprintf(stderr, "%d %d %d %d %d %d %d %d\n", (int)v0, (int)v1, (int)g_perm[v0], (int)g_perm[v1], (int)src_proc, (int)targ_proc, (int)v0_new, (int)v1_new);
      write_edge(edge, v0_new, v1_new);
    }
    free(perms);
    free(parts);
  }
  //free_oned_csr_graph(&g);
  free_graph_data_structure();
  remake_graph_data_structure(tg);
}

void partition_graph_data_structure() { 
  int n = g.nglobalverts;
  int n_local = g.nlocalverts;
  int offset = g.nlocalverts * rank; //!//Does this work?
  int nparts = size;
  int tot_nnz = 0;
  parts = (int*)malloc(n * sizeof(int));
  int *partsize_update = (int*)malloc(nparts * sizeof(int));
  int *old_partsize = (int*)malloc(nparts * sizeof(int));
  int *partsize = (int*)malloc(nparts * sizeof(int));

  int *partnnz_update = (int*)malloc(nparts * sizeof(int));
  int *old_partnnz = (int*)malloc(nparts * sizeof(int));
  int *partnnz = (int*)malloc(nparts * sizeof(int));

  int *partscore = (int*)malloc(nparts * sizeof(int));
  int *partcost = (int*)malloc(nparts * sizeof(int));
  int *vorder = (int*)malloc(n_local * sizeof(int)); 

  int oldpart;
  int emptyverts = 0;
  num_hi_deg_verts = 0;
  int randidx;
  int cutoff = F_CUTOFF;
  size_t *row;
  size_t vert;
  size_t k,  nnz_row, best_part;
  int64_t *colidx = g.column;
  int64_t node = 0;
  size_t *rowptr = g.rowstarts;
  float curr_score, best_score;
  float gamma = F_GAMMA;
  int i, s, l; //,j;

  g_perm = (int*)malloc(n * sizeof(int));
 
  if(MAT_OUT) { // Print graph
    char filename[256];
    sprintf(filename, "out_csr%02d.mat", rank);
    FILE *GraphFile;
    GraphFile = fopen(filename, "w");
    assert(GraphFile != NULL);
    print_graph_csr(GraphFile, rowptr, colidx, n_local);
    MPI_Barrier(MPI_COMM_WORLD);
    fclose(GraphFile);
  }

  memset(parts, -1, n * sizeof(int));
  for (l=0; l<nparts; ++l) {
    partsize[l] = 0;
    old_partsize[l] = 0;
    partsize_update[l] = 0;
    partnnz[l] = 0;
    old_partnnz[l] = 0;
    partnnz_update[l] = 0;
  }

  int localedges = (int)g.nlocaledges;
  MPI_Allreduce(&localedges, &tot_nnz, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  float alpha = sqrt(2) * (tot_nnz/pow(n,gamma));
  
  //fprintf(stdout,"n = %d, n_local = %d, local nnz = %d, total nnz = %d\n",n,n_local,localedges,tot_nnz);
  int repeat_run;
  //int run;
  genRandPerm(vorder, n_local);
  for (repeat_run = 0; repeat_run < NUM_STREAMS; repeat_run++) {
    //for (run=0; run<nparts; run++) {
      //if (rank == run) { //just partition one process after the other...
        //First run: initialize with hashed partition
        if (repeat_run == 0) {
          for (i = 0; i < n_local; ++i) {
            vert = (size_t)vorder[i];
            row = &rowptr[vert];
            nnz_row = *(row+1) - *row;
            int local_idx = offset + vert; 

            if (nnz_row >= cutoff) { parts[local_idx] = nparts; num_hi_deg_verts++; }
            else {
              //fprintf(stdout," %d ",local_idx);
              randidx = irand(nparts);
              //oldpart = parts[local_idx];
              parts[local_idx] = randidx;
              partsize[randidx]++;
              partnnz[randidx] += nnz_row;
            }
          }
        }
        else {    //Additional runs: run FENNEL algorithm. todo: tempering
          alpha *= ALPHA_EXP_RATE;
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
                int parts_idx = node_owner*g.nlocalverts + node_local_idx;
                int node_part = parts[parts_idx]; /////
                if (node_part >= 0 && node_part < nparts) { partscore[node_part]++; }
              }
              for (s = 0; s < nparts; ++s) { 
                //float dc = calc_dc(alpha,gamma,partsize[s]); 
                //int normscore = partscore[s] - (int)nnz_row;
                //partcost[s] = normscore - dc; 

                partcost[s] = partscore[s] - alpha*(gamma/2)*pow(partsize[s],gamma-1);
              }
              best_part = 0;
              best_score = partcost[0];
              for (s = 1; s < nparts; ++s) { 
                curr_score = partcost[s]; 
                if (curr_score > best_score) {
                  best_score = curr_score;  best_part = s;
                }
              }
              oldpart = parts[local_idx];
              parts[local_idx] = best_part;
              partsize[best_part]++; partnnz[best_part]+=nnz_row;
              if (oldpart >= 0 && oldpart < nparts) { partsize[oldpart]--; partnnz[oldpart]-=nnz_row; }
            } else { // empty vertex, assign randomly
              if (parts[local_idx]==-1) {
                emptyverts++;
                randidx = irand(nparts);
                oldpart = parts[local_idx];
                parts[local_idx] = randidx;
                partsize[randidx]++;
                if (oldpart >= 0 && oldpart < nparts) { partsize[oldpart]--; partnnz[oldpart]-=nnz_row; }
              }
            }
            if (i % 512 == 0) { //(isPowerOfTwo(i)) {
              //MPI_Allgather(parts+offset, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
              MPI_Allgather(MPI_IN_PLACE, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
              for (l=0; l<nparts; ++l) { 
                partsize_update[l] = partsize[l] - old_partsize[l];
              }            
              MPI_Allreduce(MPI_IN_PLACE, partsize_update, nparts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
              for (l=0; l<nparts; ++l) { 
                old_partsize[l] += partsize_update[l]; 
                partsize[l] = old_partsize[l];
              }
            }
          }
        }
      //}
      //MPI_Allgather(parts+offset, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
      //double allgstart= MPI_Wtime();
      //MPI_Allgather(MPI_IN_PLACE, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
      //double allgstop = MPI_Wtime();
      //if (rank == 0) { fprintf(stderr, "allgather time:               %f s\n", allgstop - allgstart); }
      for (l=0; l<nparts; ++l) { 
        partsize_update[l] = partsize[l] - old_partsize[l];
        partnnz_update[l] = partnnz[l] - old_partnnz[l];
        //fprintf(stdout,"partsize[%d] on rank %d is %d. partsizeupdate is %d. oldsize is %d\n", l, rank, partsize[l], partsize_update[l], old_partsize[l]);
      }
      MPI_Allreduce(MPI_IN_PLACE, partsize_update, nparts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, partnnz_update, nparts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      for (l=0; l<nparts; ++l) { 
        old_partsize[l] += partsize_update[l]; 
        old_partnnz[l] += partnnz_update[l]; 
        partsize[l] = old_partsize[l];
        partnnz[l] = old_partnnz[l];
      }

      //double allgstart= MPI_Wtime();
      MPI_Allgather(MPI_IN_PLACE, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
      //double allgstop = MPI_Wtime();
      //if (rank == 0) { fprintf(stderr, "allgather time:               %f s\n", allgstop - allgstart); }

    //sanity check: manually compute partsizes
    if (SANITY) {
      int max_partsize = 0;
      int min_partsize = n;
      int *check_partsize = (int*)malloc((nparts+1)*sizeof(int));
      int max_partnnz = 0;
      int min_partnnz = tot_nnz;
      for (l=0; l<=nparts; ++l) { 
        check_partsize[l] = 0; 
      }
      for (i=0; i<n; ++i) {
        int mypart = parts[i];
        if (mypart == -1) { fprintf(stderr, "-1 :("); }
        assert(mypart>=0 && mypart<=nparts);
        check_partsize[mypart]++;
      }
      for (l=0; l<nparts; ++l) { 
        if (rank==l && VERBY) { fprintf(stdout,"partsize[%d] on rank %d is %d. check_partsize is %d\n", l, rank, partsize[l], check_partsize[l]); }
        assert(check_partsize[l] == partsize[l]); 
        if (check_partsize[l] > max_partsize) { max_partsize = check_partsize[l]; }
        if (check_partsize[l] < min_partsize) { min_partsize = check_partsize[l]; }
        if (partnnz[l] > max_partnnz) { max_partnnz = partnnz[l]; }
        if (partnnz[l] < min_partnnz) { min_partnnz = partnnz[l]; }
        //if (rank==0) { fprintf(stdout,"%d / %d\n",max_partsize, min_partsize); }
        //if (rank==0) { fprintf(stdout,"%d / %d\n",max_partnnz, min_partnnz); }
      }
      //if (rank==0) { fprintf(stdout,"max partsize = %d, min partsize = %d max/min partnnz = %d, %d ", max_partsize, min_partsize, max_partnnz, min_partnnz); }
      if (rank==0) { fprintf(stdout,"n balance: %f, nnz balance: %f\t", (float)max_partsize / min_partsize, (float)max_partnnz / min_partnnz); }
    }

    if (1 || SANITY) {
      mpi_compute_cut(rowptr, colidx, parts, nparts, n_local, offset, cutoff);
    }
  }

  // Okay this is absurd but for now, randomly assign the large vertices to a permutation
  // To prevent overly large partitions
  for (i=offset; i<offset+n_local; ++i) {
    if (parts[i] == nparts) {
      if (HI_RAND) { parts[i] = irand(nparts); }
      else { parts[i] = nparts; }
    }
  }
  MPI_Allgather(MPI_IN_PLACE, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
  //MPI_Allgather(parts+offset, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (MAT_OUT &&(rank == 0)) {  // Print Parts
    FILE *PartFile;
    PartFile = fopen("parts.mat", "w");
    assert(PartFile != NULL);
    print_parts(PartFile, parts, n, n_local);
    fclose(PartFile);
  }

  if (SANITY) {
    // Sanity Checks
    for (i=0; i < g.nglobalverts; ++i) {
      assert(parts[i]>=0);
      assert(parts[i]<nparts);
    }
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
  return (alpha*pow(len + F_DELTA ,gamma)) - (alpha*pow(len,gamma));
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
  if (rank == 0) {   fprintf(stdout,"total cutedges = %d, pct of total:%f pct of worstcase:%f \n", tot_cutedges, (float)tot_cutedges/tot_lodegedges, ((float)tot_cutedges/tot_lodegedges)/((float)(nparts-1)/nparts)); }
  return tot_cutedges;
}

int print_graph_tuple(FILE* out, tuple_graph* tg, int rank) {
  int64_t v0, v1;
  packed_edge* result = tg->edgememory;
  packed_edge* edge;
  int i;
  for (i=0; i < tg->edgememory_size; ++i) {
    edge = &result[i];
    v0 = get_v0_from_edge(edge);
    v1 = get_v1_from_edge(edge);
    fprintf (out, "%d %d %d\n", (int)v0+1, (int)v1+1, 1);
  }
  return 1;
}

int print_graph_csr(FILE* out, size_t *rowptr, int64_t *colidx, int n_local) {
  int i, k;
  int64_t src, dst;
  for (i=0; i<n_local; ++i) {
    for (k = rowptr[i]; k < rowptr[i+1]; ++k) {
      src = n_local*rank + i;
      //src = VERTEX_TO_GLOBAL(rank,i);
      //colidxs are correct, but to correctly vis
      //communication, we need to map them
      dst = colidx[k]; 
      dst = ((dst-VERTEX_OWNER(dst))/size) + n_local*VERTEX_OWNER(dst);
      //src = (int)VERTEX_TO_GLOBAL(rank,i);
      //fprintf (out, "%d %d %d\n",(int)VERTEX_TO_GLOBAL(rank,i)+1,(int)dst+1,1);
      fprintf (out, "%d %d %d\n",(int)src+1,(int)dst+1,1);
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
