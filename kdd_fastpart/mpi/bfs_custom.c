
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
      PART_TYPE v_part = parts[parts_idx];
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
  size_t n = g.nglobalverts;
  size_t n_local = g.nlocalverts;
  size_t offset = g.nlocalverts * rank; //!//Does this work?
  int nparts = size;
  int64_t tot_nnz = 0;
  parts = (PART_TYPE*)malloc(n * sizeof(PART_TYPE));
  int64_t *partsize_update = (int64_t*)malloc(nparts * sizeof(int64_t));
  int64_t *old_partsize = (int64_t*)malloc(nparts * sizeof(int64_t));
  int64_t *partsize = (int64_t*)malloc(nparts * sizeof(int64_t));

  int64_t *partnnz_update = (int64_t*)malloc(nparts * sizeof(int64_t));
  int64_t *old_partnnz = (int64_t*)malloc(nparts * sizeof(int64_t));
  int64_t *partnnz = (int64_t*)malloc(nparts * sizeof(int64_t));

  int64_t *partscore = (int64_t*)malloc(nparts * sizeof(int64_t));
  int64_t *partcost = (int64_t*)malloc(nparts * sizeof(int64_t));
  int64_t *vorder = (int64_t*)malloc(n_local * sizeof(int64_t)); 

  PART_TYPE oldpart;
  int64_t emptyverts = 0;
  num_hi_deg_verts = 0;
  PART_TYPE randidx;
  int cutoff = F_CUTOFF;
  size_t *row;
  size_t vert;
  size_t k,  nnz_row, best_part;
  int64_t *colidx = g.column;
  size_t *rowptr = g.rowstarts;
  float curr_score, best_score;
  float gamma = F_GAMMA;
  size_t i, s, l; //,j;

  //g_perm = (int64_t*)malloc(n * sizeof(int64_t));
 
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

  //memset(parts, -1, n * sizeof(PART_TYPE));

  for (l=0; l<nparts; ++l) {
    partsize[l] = 0;
    old_partsize[l] = 0;
    partsize_update[l] = 0;
    partnnz[l] = 0;
    old_partnnz[l] = 0;
    partnnz_update[l] = 0;
  }

  int64_t localedges = (int64_t)g.nlocaledges;
  MPI_Allreduce(&localedges, &tot_nnz, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  float alpha = sqrt(2) * (tot_nnz/pow(n,gamma));
  
  //fprintf(stdout,"n = %d, n_local = %d, local nnz = %d, total nnz = %d\n",n,n_local,localedges,tot_nnz);
  int64_t repeat_run;
  int64_t parts_idx;
  int64_t local_idx;
  double allgpartstime = 0;
  double streamtime = 0;
  genRandPerm(vorder, n_local);

  double streamstart= MPI_Wtime();
  for (repeat_run = 0; repeat_run < NUM_STREAMS; repeat_run++) {
    //First run: initialize with hashed partition
    if (repeat_run == 0) {
      for (i = 0; i < n_local; ++i) {
        vert = (size_t)vorder[i];
        row = &rowptr[vert];
        nnz_row = *(row+1) - *row;
        local_idx = offset + vert; 
        randidx = (PART_TYPE)irand(nparts);
        parts[local_idx] = randidx;
        partsize[randidx]++;
        partnnz[randidx] += nnz_row;
      }
    }
    else {
      alpha *= ALPHA_EXP_RATE;
      for (i = 0; i < n_local; ++i) {
        memset(partscore, 0, nparts * sizeof(int64_t));
        memset(partcost, 0, nparts * sizeof(int64_t));
        vert = (size_t)vorder[i];
        local_idx = offset + vert; //VERTEX_LOCAL(global_vert_idx);
        row = &rowptr[vert];
        nnz_row = *(row+1) - *row;
        oldpart = -1;
        if(nnz_row > 0) {
          for (k = *row; k < ((*row)+nnz_row); ++k) {
            parts_idx = VERTEX_OWNER(colidx[k])*g.nlocalverts + VERTEX_LOCAL(colidx[k]);
            PART_TYPE node_part = parts[parts_idx]; /////
            if (parts[parts_idx] >= 0) { partscore[node_part]++; }
          }
          for (s = 0; s < nparts; ++s) { partcost[s] = partscore[s] - alpha*(gamma/2)*pow(partsize[s],gamma-1); }
          best_part = 0;
          best_score = partcost[0];
          for (s = 1; s < nparts; ++s) { 
            curr_score = partcost[s]; 
            if (curr_score > best_score) { best_score = curr_score;  best_part = s; }
          }
          oldpart = parts[local_idx];
          parts[local_idx] = best_part;
          partsize[best_part]++; partnnz[best_part]+=nnz_row;
          if (oldpart >= 0 && oldpart < nparts) { partsize[oldpart]--; partnnz[oldpart]-=nnz_row; }
        } 
#if 0
        else { // empty vertex, assign randomly
          if (parts[local_idx]==-1) {
            emptyverts++;
            randidx = (PART_TYPE)irand(nparts);
            oldpart = parts[local_idx];
            parts[local_idx] = randidx; partsize[randidx]++;
            if (oldpart >= 0 && oldpart < nparts) { partsize[oldpart]--; partnnz[oldpart]-=nnz_row; }
          }
        }
#endif
        if (i % (n_local / 2000) == 0) {
          for (l=0; l<nparts; ++l) { partsize_update[l] = partsize[l] - old_partsize[l]; }            
          MPI_Allreduce(MPI_IN_PLACE, partsize_update, nparts, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
          for (l=0; l<nparts; ++l) { old_partsize[l] += partsize_update[l]; partsize[l] = old_partsize[l]; }
        }
      }
    }
    for (l=0; l<nparts; ++l) { 
      partsize_update[l] = partsize[l] - old_partsize[l];
      partnnz_update[l] = partnnz[l] - old_partnnz[l];
      //fprintf(stdout,"partsize[%d] on rank %d is %d. partsizeupdate is %d. oldsize is %d\n", l, rank, partsize[l], partsize_update[l], old_partsize[l]);
    }
    MPI_Allreduce(MPI_IN_PLACE, partsize_update, nparts, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, partnnz_update, nparts, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    for (l=0; l<nparts; ++l) { 
      old_partsize[l] += partsize_update[l]; 
      old_partnnz[l] += partnnz_update[l]; 
      partsize[l] = old_partsize[l];
      partnnz[l] = old_partnnz[l];
    }

    double allgstart= MPI_Wtime();
    MPI_Allgather(MPI_IN_PLACE, n_local, MPI_PART_TYPE, parts, n_local, MPI_PART_TYPE, MPI_COMM_WORLD);
    double allgstop = MPI_Wtime();
    allgpartstime += allgstop - allgstart;
    //if (rank == 0) { fprintf(stderr, "allgather time:               %f s\n", allgstop - allgstart); }

    //sanity check: manually compute partsizes
    if (1){//(0 || SANITY) {
      int64_t max_partsize = 0;
      int64_t min_partsize = n;
      int64_t *check_partsize = (int64_t*)malloc((nparts+1)*sizeof(int64_t));
      int64_t max_partnnz = 0;
      int64_t min_partnnz = tot_nnz;
      for (l=0; l<=nparts; ++l) { 
        check_partsize[l] = 0; 
      }
      for (i=0; i<n; ++i) {
        PART_TYPE mypart = parts[i];
        if (mypart == -1) { fprintf(stderr, "-1 :("); }
        assert(mypart>=0 && mypart<=nparts);
        check_partsize[mypart]++;
      }
      for (l=0; l<nparts; ++l) { 
        if (rank==l && VERBY) { fprintf(stdout,"partsize[%d] on rank %d is %" PRId64 ". check_partsize is %" PRId64 "\n", (int)l, (int)rank, partsize[l], check_partsize[l]); }
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

    if (1) {//(repeat_run % 2 == 0 && repeat_run > 0 || SANITY) {
      mpi_compute_cut(rowptr, colidx, parts, nparts, n_local, offset, cutoff);
    }
  }
  double streamstop= MPI_Wtime();
  streamtime += streamstop - streamstart;
  if (rank == 0) { fprintf(stderr, "stream time: %f,  per-stream time: %f \n", streamtime, streamtime/NUM_STREAMS); }

#if 0
  for (i=offset; i<offset+n_local; ++i) {
    if (parts[i] == nparts) {
      if (HI_RAND) { parts[i] = irand(nparts); }
      else { parts[i] = nparts; }
    }
  }
#endif
  MPI_Allgather(MPI_IN_PLACE, n_local, MPI_PART_TYPE, parts, n_local, MPI_PART_TYPE, MPI_COMM_WORLD);
  //MPI_Allgather(parts+offset, n_local, MPI_INT, parts, n_local, MPI_INT, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
#if 0
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
#endif
  if (rank == 0) { fprintf(stderr, "allgather parts time:               %f s\n", allgpartstime); }
  //if (rank == 0) { fprintf(stderr, "comp  time:               %f s\n", comptime); }

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
