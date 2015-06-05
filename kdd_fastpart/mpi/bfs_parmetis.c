
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

#include <parmetis.h>

///////////////////
//#define TO_NEW_IDX(x) g_perm[x]
//#define NEW_PART_OF_IDX(x) parts[TO_NEW_IDX(x)]
///////////////////

static int64_t* g_perm;
static PART_TYPE* parts;
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

static int64_t num_hi_deg_verts;

void free_graph_data_structure(void) {
  /*free(g_oldq);
  free(g_newq);
  free(g_visited);
  MPI_Free_mem(g_outgoing);
  free(g_outgoing_counts);
  free(g_outgoing_reqs);
  free(g_outgoing_reqs_active);
  MPI_Free_mem(g_recvbuf);*/
  free_oned_csr_graph(&g);
}

void make_graph_data_structure(const tuple_graph* const tg) {
  convert_graph_to_oned_csr(tg, &g);
//  const size_t nlocalverts = g.nlocalverts;
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
void permute_tuple_graph(tuple_graph* tg) { }

void partition_graph_data_structure() { 
  fprintf(stdout, " %d ", (int)sizeof(idx_t));
  size_t n = g.nglobalverts;
  size_t n_local = g.nlocalverts;
  int64_t localedges = (int64_t)g.nlocaledges;

  size_t offset = g.nlocalverts * rank; //!//Does this work?
  int64_t tot_nnz = 0;
  //parts = (PART_TYPE*)malloc(n * sizeof(PART_TYPE));

  size_t k,  nnz_row, best_part;
  int64_t *colidx = g.column;
  size_t *rowptr = g.rowstarts;
  size_t i, s, l; //,j;

  // need to convert colidx, rowptr to 32-bit ints.
  idx_t* colidx_32 = (idx_t*)malloc(localedges*sizeof(idx_t));
  idx_t* rowptr_32 = (idx_t*)malloc((n_local+1)*sizeof(idx_t));
  for (i=0; i<localedges; ++i) { colidx_32[i] = (idx_t)colidx[i]; }
  for (i=0; i<=n_local; ++i) { rowptr_32[i] = (idx_t)rowptr[i]; }

  int result;
// Needed by parmetis
 
  //idx_t *vtxdist=NULL;
  idx_t *xadj=NULL;
  idx_t *adjncy=NULL;
  idx_t *vwgt=NULL, *adjwgt=NULL;
  idx_t wgtflag=0;
  idx_t numflag=0;
  idx_t ncon=1;
  idx_t nparts=size;
  //real_t *tpwgts=NULL, 

  real_t ubvec;
  idx_t options[4], edgecut;
  idx_t part[n_local];
  MPI_Comm comm;

  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  MPI_Allreduce(&localedges, &tot_nnz, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

  //fprintf(stdout, " %d ", (int)sizeof(idx_t));

  double streamstart= MPI_Wtime();
  idx_t vtxdist[size+1];// = new idx_t[size];
// For AdaptiveRepart
  float itr = 1000.0;
  idx_t *vsize=NULL;
  idx_t vert_so_far = 0;
  for (i=0; i<=size; ++i) { vtxdist[i] = vert_so_far; vert_so_far+=n_local; }

  for (i=0; i<size; ++i) { fprintf(stdout, " %d ", (int)vtxdist[i]);}
    fprintf(stdout,"\n");
  //MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT64_T, vtxdist, 1, MPI_INT64_T, MPI_COMM_WORLD);

  /*
  vtxdist[0] = 0;
  vtxdist[1] = 5;
  vtxdist[2] = 10;
  vtxdist[3] = 15;
*/
  ubvec = 1.05;

  options[0] = 0;
  options[1] = 0;
  options[2] = 0;
  options[3] = 0;

  real_t tpwgts[size];

  for (i=0; i<n_local; ++i) { part[i] = rank; }
  for (i=0; i<size; ++i) { tpwgts[i] = 1.0/(float)size; }
  xadj = rowptr_32;
  adjncy = colidx_32;
  //xadj = new idx_t[6];
  //adjncy = new idx_t[13];

  if ( rank == 0 ){
    fprintf(stdout,"parmetis initialized.\n");
  }
 
result = ParMETIS_V3_PartKway( vtxdist, xadj, adjncy, vwgt, adjwgt, 
                                 &wgtflag, &numflag, &ncon, 
                                 &nparts, tpwgts, &ubvec, options, 
                                 &edgecut, part, &comm );  
/*result = ParMETIS_V3_AdaptiveRepart( vtxdist, xadj, adjncy, vwgt, vsize, 
                                 adjwgt, &wgtflag, &numflag, &ncon, 
                                 &nparts, tpwgts, &ubvec, &itr, options, 
                                 &edgecut, part, &comm );  */


  double streamstop = MPI_Wtime();
  double streamtime = streamstop - streamstart;
  if (rank == 0) { fprintf(stderr, "stream time: %f,  per-stream time: %f \n", streamtime, streamtime/NUM_STREAMS); }
  for (i=0;i<size;++i) {
  MPI_Barrier(MPI_COMM_WORLD);
    if ( rank == i ){
      fprintf(stderr,"%d edgecut: %d\n", rank, edgecut);// MPI_PROC_ID << " edgecut " << edgecut << '\n';
    }
}
  //MPI_Allgather(MPI_IN_PLACE, n_local, MPI_PART_TYPE, parts, n_local, MPI_PART_TYPE, MPI_COMM_WORLD);
}

// Random permutation generator. Move to another file.
int64_t* genRandPerm(int64_t* orderList, int64_t size) {
  assert(orderList);
  srand(time(NULL));
  // Generate 'identity' permutation
  int i;
  for (i = 0; i < size; i++) { orderList[i] = i; }
  shuffle_int(orderList, size);
  return orderList;
}

void shuffle_int(int64_t *list, int len) {
  int j;
  int64_t tmp;
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

float calc_dc(float alpha, float gamma, int64_t len) {
  return (alpha*pow(len + F_DELTA ,gamma)) - (alpha*pow(len,gamma));
}

int64_t mpi_compute_cut(size_t *rowptr, int64_t *colidx, PART_TYPE* parts, int nparts, int64_t n_local, int64_t offset, int cutoff) {
  size_t vert;
  int64_t nnz_row;
  int64_t v_part;
  int64_t cutedges = 0;
  int64_t mytotedges = 0;
  int64_t mytotlodegedges = 0;
  size_t *row;
  int64_t i;
  size_t k;
  int64_t emptyparts = 0;
  mytotedges = rowptr[n_local];
  for (i = 0; i < n_local; i++) {
    vert = i;
    row = &rowptr[vert];
    nnz_row = (int64_t)(*(row+1) - *(row)); //nnz in row
    if (nnz_row < cutoff) { 
      v_part = parts[vert+offset];
      if (v_part == -1) {
        v_part = 0;
        emptyparts++;
      }
      // count edges to other partitions
      for (k = *row; k < ((*row)+nnz_row); ++k) {
        int64_t node = colidx[k];
        int64_t node_owner = VERTEX_OWNER(node);
        int64_t node_local_idx = VERTEX_LOCAL(node);
        int64_t parts_idx = node_owner*g.nlocalverts + node_local_idx;
        if (parts[parts_idx] < nparts) { mytotlodegedges++; } //count low degree edges
        if (parts[parts_idx] != v_part && parts[parts_idx] < nparts) { cutedges++; } //count low degree cut edges
      }
    }
  }
  int64_t tot_cutedges;
  int64_t tot_lodegedges;
  MPI_Allreduce(&cutedges, &tot_cutedges, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mytotlodegedges, &tot_lodegedges, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  //fprintf(stdout,"offset: %d emptyparts = %d cutedges = %d totcutedges = %d tot edges=%d mylodegedges=%d totlodegedges=%d\n",offset, emptyparts,cutedges,tot_cutedges,mytotedges,mytotlodegedges,tot_lodegedges);
  if (rank == 0) {   fprintf(stdout,"total cutedges = %" PRId64 ", pct of total:%f pct of worstcase:%f \n", tot_cutedges, (float)tot_cutedges/tot_lodegedges, ((float)tot_cutedges/tot_lodegedges)/((float)(nparts-1)/nparts)); }
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

int print_parts(FILE* out, PART_TYPE* parts, int n, int n_local) {
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

