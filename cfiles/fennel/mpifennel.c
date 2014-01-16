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

#include <mpi.h>

#include "mpifennel.h"
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
mpi_usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args) {
  fprintf (out, "Usage:\n");
  fprintf (out, "mpifennel [options] <in-filename> <in-format>\n");
  fprintf (out, "<in-filename>: name of file containing the sparse matrix to read in\n");
  fprintf (out, "<in-format>: format of the input file (\"HB\", \"ML\", \"MM\")\n");
  fprintf (out, "<out-filename>: name of file to which to output results\n");
  fprintf (out, "[options]: -e -- expand symmetric into unsymmetric storage\n");
  fprintf (out, " -v -- verbose mode\n");
  fprintf (out, "EX: mpirun -np 4 ./mpifennel -v -e 'test.mtx' 'MM'\n");
}

int main (int argc, char *argv[]) {
  MPI_Init(NULL, NULL);

  // MPI STATS ///////
  int numprocs, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  
  int m,n,nnz;
  void *values;
  int *colidx, *rowptr;
  enum symmetry_type_t symmetry_type;
  enum symmetric_storage_location_t symmetric_storage_location;
  enum value_type_t value_type;
  
  // For now we are only simulating a parallel file system.
  // The root process reads in the data and then scatters the graph
  if (rank == 0) {
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
        MPI_Abort(MPI_COMM_WORLD,-1);
        bebop_exit (EXIT_FAILURE);
    }

    const char* matrix_filename = argv[1];
    register_usage_function (mpi_usage);
    arglist = register_arginfo (arglist, 'v', NULLARG, NULL, "If specified, ac"
                           "tivate verbose mode");
                         

    arglist = register_arginfo (arglist, 'e', NULLARG, NULL, "If specified, ex"
                           "pand the input matrix from symmetric storage "
                           "into unsymmetric storage");
    arglist = register_arginfo (arglist, 'a', NULLARG, NULL, "If specified, va"
                           "lidate the input matrix, without outputting a"
                           "nything");

    get_options (argc, argv, arglist, NULL);

    if (argc - optind != 2) {
        fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
           "e specified, but there should be %d ***\n", argc - optind, 2);
        dump_usage (stderr, argv[0], arglist, NULL);
        MPI_Abort(MPI_COMM_WORLD,-1);
        bebop_exit (EXIT_FAILURE); /* stops logging */
    }

    opts.input_filename = argv[optind];
    opts.input_file_format = argv[optind+1];
    int parts = rank; //atoi(argv[optind+2]);
  
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
        MPI_Abort(MPI_COMM_WORLD,-1);
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
        MPI_Abort(MPI_COMM_WORLD,-1);
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
        MPI_Abort(MPI_COMM_WORLD,-1);
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

    unpack_csr_matrix (A->repr, &m, &n, &nnz, &values, &colidx, &rowptr, &symmetry_type,
                 &symmetric_storage_location, &value_type);
                 
  }
  
  // Only process 0 has n, so broadcast ...
  MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nnz,1,MPI_INT,0,MPI_COMM_WORLD); 
  MPI_Barrier(MPI_COMM_WORLD);

  // Allocate local rowptr so that we can scatter global rowptr
  int p = numprocs;
  int bound = ceil((float)n/p);
  int nnz_local;
  int* ir_local = (int*)malloc(sizeof(int)*bound);
  
  int* startidxs;
  int* endidxs;
  int* nnz_per_proc;
  int* ir2;
  
  // Now scatter repr to all processes from root
  if (rank == 0) {
    // bound = ceil(N/p)
    // Resize ir to p*bound
    // Send ir[rank*bound .. (rank+1)*bound - 1] to p_rank  (Scatter of width bound)
    // Send cidx[ir[rank*bound] .. ir[(rank+1)*bound - 1]] to p_rank  (Using MPI_Send ... )
    
    // Create ir2, a resized rowptr
    ir2 = (int*)malloc(sizeof(int) * p*bound);
    
    // Copy rowptr to ir2
    // Todo: switch to memcpy
    int i;
    for (i = 0; i<n; i++) { ir2[i] = rowptr[i]; }
    for (i = n; i<bound*p; i++) { ir2[i] = rowptr[n]; } //extend ir2 to have empty nodes
    
    // Scatter row ptr
    //MPI_Scatter(ir2, bound, MPI_INT, ir_local, bound, MPI_INT, 0, MPI_COMM_WORLD);
  }
  
  //fprintf(stdout,"bound = %d\n",bound);
  
  if (rank == 0) { fprintf(stdout,"=== Sending vertices to processes ===\n"); }

  MPI_Scatter(ir2, bound, MPI_INT, ir_local, bound, MPI_INT, 0, MPI_COMM_WORLD);
  assert(n);
  assert(p);

  if (rank == 0) {
    // OK.. now to scatter the colidx we need to first compute the subindices that
    // belong to each process... we'll store this in an array that we can then scatter
    // so that processes can allocate their slices
    
    startidxs = (int*)malloc(sizeof(int)*p);
    endidxs = (int*)malloc(sizeof(int)*p);
    nnz_per_proc = (int*)malloc(sizeof(int)*p);

    for (int i = 0; i<p; i++) { 
      int startvert = i*bound; 
      int endvert = (i+1)*bound - 1;
      startidxs[i] = (i > 0) ? ir2[startvert] : 0;
      endidxs[i] = ir2[endvert];
      nnz_per_proc[i] = endidxs[i] - startidxs[i];
      //fprintf(stdout,"Bound for p%d from %d to %d \n",i,startvert,endvert);
      //fprintf(stdout,"startidxs[%d] = %d\n",i,startidxs[i]);
    }
  }
  
  MPI_Scatter(nnz_per_proc, 1, MPI_INT, &nnz_local, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // Every processor should have nnz_local at this point, so allocate local colidx
  // so that we can scatter colidx
  int* colidx_local = (int*)malloc(sizeof(int)*nnz_local);

  // Now scatter colidx to all processes from root
  if (rank == 0) {
    fprintf(stdout,"=== Sending nnzs to processes ===\n");
    for (int i=0; i<p; i++) {
      //fprintf(stdout,"Sending from %d to %d, offset %d, length %d\n",rank,i,startidxs[i],nnz_per_proc[i]);
      MPI_Send(colidx + startidxs[i], nnz_per_proc[i], MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
  MPI_Recv(colidx_local, nnz_local, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  

  // ir_local points to colidx idxs, so we need to normalize
  // so that it points to colidx_local
  int ir_offset = ir_local[0];
  for (int i=0; i<bound; i++) { ir_local[i] -= ir_offset; }
  // Reduce nnz_local to assert that it is the same
  
  //////////////// SANITY CHECKS ///////////////////////////////////////
  //MPI_Barrier(MPI_COMM_WORLD);
  //fprintf(stdout,"First nz on proc %d: %d. global node id=%d\n",rank,ir_offset,rank*bound);
  int tot_nnz_for_assert;
  MPI_Reduce(&nnz_local, &tot_nnz_for_assert, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  //if (rank == 0) { fprintf(stdout,"Total nnz: %d\n",tot_nnz_for_assert); }
  int sum = 0;
  int tot_counted_nnz_for_assert;
  for (int i=0; i<bound-1; i++) { sum += ir_local[i+1] - ir_local[i]; }
  MPI_Reduce(&sum, &tot_counted_nnz_for_assert, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  //fprintf(stdout,"IR counted nnz on proc %d: %d, difference of %d\n",rank,sum,nnz_local-sum);
  if (nnz_local-sum != 0) { 
    fprintf(stdout,"===== ERROR: not all nnzs were communicated ===\n");
    MPI_Abort(MPI_COMM_WORLD,-1);
    bebop_exit (EXIT_FAILURE);
  }
  for (int i=0; i<nnz_local; i++) { 
    if (colidx_local[i] < 0 || colidx_local[i] > n) { fprintf(stdout,"%d ",colidx_local[i]); }
  }
  //////////////////////////////////////////////////////////////////////

  int local_vertex_offset = rank*bound; // add to all local vertex id's to get global id
  // ********** Run MPI FENNEL ***************************************
  if (rank==0) { fprintf (stdout, "===== Running mpi fennel =====\n"); }
  mpi_run_fennel(ir_local, colidx_local, n, bound, nnz, nnz_local, local_vertex_offset, 1.5); //todo: nparts, gamma as inputs
  // *************************************************************
  //destroy_sparse_matrix (A);
  
  /* END */
  MPI_Finalize();
  return 0;
}

int mpi_fennel_kernel(int n, int n_local, int offset, int nparts, int *partsize, 
    int *rowptr, int *colidx, int **parts, float alpha, float gamma, int *emptyverts) {
      
  int *partscore = (int*)malloc(nparts * sizeof(int));
  int *row;
  int vert, k, s, nnz_row, best_part, randidx, nededges = 0, node = 0;
  float curr_score, best_score;
  int *vorder = genRandPerm(n_local);
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
        for (s = 0; s < nparts; s++) { if (parts[s][node] == 1) { partscore[s]++; }}
      }
        
      // choose optimal partition (initializing first)
      best_score = (partscore[0]-nnz_row) - calc_dc(alpha,gamma,partsize[0]);
      best_part = 0;
      for (s = 1; s < nparts; s++) {
        curr_score = (partscore[s]-nnz_row) - calc_dc(alpha,gamma,partsize[s]);
        if (curr_score > best_score) { best_score = curr_score; best_part = s; }
      }
      for (s = 0; s < nparts; s++) { 
        if (parts[s][offset + vert] == 1) {
          oldpart = s;
        }
        parts[s][offset + vert] = 0; 
      }
      parts[best_part][offset + vert] = 1;
      //int sum=0;
      //for (s = 0; s < nparts; s++) { sum += parts[s][offset + vert]; }
      //assert(sum==1);
      partsize[best_part]++;
      if (oldpart >= 0) {
        partsize[oldpart]--;
      }
      
    } else { // empty vertex for some reason... assign it to random permutation
      emptyverts++;
      randidx = irand(nparts);
      for (s = 1; s < nparts; s++) {
        if (parts[s][offset + vert] == 1) {
          oldpart = s;
        }
        parts[s][offset + vert] = 0; 
      }
      parts[randidx][offset + vert] = 1;
      partsize[randidx]++;
      if (oldpart >= 0) {
        partsize[oldpart]--;
      }
    }
  }
  
  free(partscore);
  free(vorder);
}

int mpi_run_fennel(int* rowptr, int* colidx, int n, int n_local, 
    int nnz, int nnz_local, int v_offset, float gamma) {

  int nparts, rank;
  float seconds = 0.0f;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nparts);

  float alpha = sqrt(2) * (nnz/pow(n,gamma));
  //float alpha = nnz * pow(2,(gamma-1))/pow(m,gamma);
  if (rank==0) { fprintf (stdout, "p0:----> gamma = %f, alpha = %f\n",gamma,alpha); }
          
  // Allocate partition vectors
  if (rank==0) { fprintf (stdout, "p0:----> Gen %d partition vectors\n",nparts); }
  int** parts = (int**)malloc(nparts * sizeof(int*));
  for(int i = 0; i < nparts; i++) {
    parts[i] = (int*)malloc(n * sizeof(int));
    for( int j=0; j<n; j++ ) { parts[i][j] = 0; } //fill with zeros
    assert(parts[i]);
  }
  assert(parts);

  // Iterate over vorder
  if (rank==0) { fprintf (stdout, "p0:----> Begin outer loop... \n"); }
  int vert, k, s, node = 0;
  int *row;
  int nnz_row = 0;
  int emptyverts;
  
  //initialize part sizes
  int *partsize = (int*)malloc(nparts * sizeof(int));
  for (s = 0; s < nparts; s++) { partsize[s] = 0; }

  seconds = get_seconds();
  mpi_fennel_kernel(n, n_local, v_offset, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts);
  seconds = get_seconds() - seconds;
  MPI_Barrier(MPI_COMM_WORLD);

  // Compute load balance
  int max_load = partsize[0];
  int min_load = partsize[0];
  for (s = 1; s < nparts; s++) {
    if (partsize[s] > max_load) {max_load = partsize[s];}
    if (partsize[s] < min_load) {min_load = partsize[s];}
  }
 
  //fprintf (stdout, "\n===== Fennel Complete in %g seconds =====\n", seconds);
  fprintf (stdout, "----> Partition sizes: ");
  for (s = 0; s < nparts; s++) { fprintf (stdout, "| %d |", partsize[s]); }
  MPI_Barrier(MPI_COMM_WORLD);
  fprintf (stdout, "\n----> Load balance: %d / %d = %1.3f\n",max_load,min_load,(float)max_load/min_load);

  // Compute cut quality
  int cutedges = 0;
  int emptyparts = 0; //empty assignments
  int redparts = 0; //redundant assignments

  MPI_Barrier(MPI_COMM_WORLD);
  //for (s=0; s<nparts; s++) {
     //MPI_Allgather(parts[0], n, MPI_INT, parts+v_offset, n/nparts, MPI_INT, MPI_COMM_WORLD);
  //}
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank==-1) {
  cutedges = mpi_compute_cut(&emptyparts, &redparts, rowptr, colidx, parts, nparts, n, NULL);

  fprintf (stdout, "----> Percent edges cut = %d / %d = %1.3f\n",cutedges,nnz,(float)cutedges/nnz);
  fprintf (stdout, "----> Percent of random: %1.3f\n\n",((float)cutedges/nnz)/((float)(nparts-1)/nparts));
  fprintf (stdout, "===== Sanity Check =====\n");
  fprintf (stdout, "----> Unassigned vertices (error): %d\n",emptyparts);
  fprintf (stdout, "----> Overassigned vertices (error): %d\n", redparts);
  fprintf (stdout, "----> Empty vertices: %d\n", emptyverts);
  }
/*  
  //FILE *PartFile;
  //PartFile = fopen("parts.mat", "w");
  //assert(PartFile != NULL);
  
  //FILE *LambdaFile;
  //LambdaFile = fopen("lambda.txt", "a");
  //assert(LambdaFile != NULL);

  
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
    mpi_fennel_kernel(n, n_local, v_offset, nparts, partsize, rowptr, colidx, parts, alpha, gamma, &emptyverts);
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
  */
}
  
static int mpi_compute_cut(int *emptyparts, int *redparts, int *rowptr, int *colidx, int **parts, int nparts, int n, FILE *out) {
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