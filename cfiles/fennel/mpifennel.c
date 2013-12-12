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
  fprintf (out, "EX: (MPIRUN) ./mpifennel -v -e 'test.mtx' 'MM' 4\n");
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
  
  // For now we are simulating a parallel file system, so root process reads
  // in the data and then scatters the graph
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

    if (argc - optind != 3) {
        fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
           "e specified, but there should be %d ***\n", argc - optind, 3);
        dump_usage (stderr, argv[0], arglist, NULL);
        bebop_exit (EXIT_FAILURE); /* stops logging */
    }

    opts.input_filename = argv[optind];
    opts.input_file_format = argv[optind+1];
    int parts = atoi(argv[optind+2]);
    //opts.output_filename = argv[optind+2];
  
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
                 
    // Only process 0 has n, so broadcast ...
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD); 
  }
  
  MPI_Barrier(MPI_COMM_WORLD);


  // Every processor should have n at this point, so allocate local ir
  // so that we can scatter ir
  int p = numprocs;
  int bound = ceil(n/p);
  int nnz_local;
  int* ir_local = (int*)malloc(sizeof(int)*bound);
  
  int* startidxs;
  int* endidxs;
  int* nnz_per_proc;

  // Now scatter repr to all processes from root
  if (rank == 0) {
    // bound = ceil(N/p)
    // Resize ir to p*bound
    // Send ir[rank*bound .. (rank+1)*bound - 1] to p_rank  (Scatter of width bound)
    // Send cidx[ir[rank*bound] .. ir[(rank+1)*bound - 1]] to p_rank  (Using MPI_Send ... )
    
    // Create ir2, a resized rowptr
    int* ir2 = (int*)malloc(sizeof(int) * p*bound);
    
    // Copy rowptr to ir2
    // Todo: switch to memcpy
    int i;
    for (i = 0; i<n; i++) { ir2[i] = rowptr[i]; }
    for (i = n; i<bound*n; i++) { ir2[i] = rowptr[n-1]; } //extend ir2 to have empty nodes
    
    // Scatter row ptr
    MPI_Scatter(ir2, bound, MPI_INT, ir_local, bound, MPI_INT, 0, MPI_COMM_WORLD);
}
/*
    // OK.. now to scatter the colidx we need to first compute the subindices that
    // belong to each process... we'll store this in an array that we can then scatter
    // so that processes can allocate their slices
    
    startidxs = (int*)malloc(sizeof(int)*p);
    endidxs = (int*)malloc(sizeof(int)*p);
    nnz_per_proc = (int*)malloc(sizeof(int)*p);

    for (i = 0; i<p; i++) { 
      int startvert = p*bound; 
      int endvert = (p+1)*bound - 1;
      startidxs[i] = ir2[startvert];
      endidxs[i] = ir2[endvert];
      nnz_per_proc[i] = endidxs[i] - startidxs[i];
    }
    
    MPI_Scatter(nnz_per_proc, 1, MPI_INT, &nnz_local, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }


  // Every processor should have nnz_local at this point, so allocate local colidx
  // so that we can scatter colidx
  int* colidx_local = (int*)malloc(sizeof(int)*nnz_local);
  
  // Now scatter colidx to all processes from root
  if (rank == 0) {
    int i;
    for (i=0; i<p; i++) {
      MPI_Send(colidx + startidxs[i], nnz_per_proc[i], MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
  
  MPI_Recv(colidx_local, nnz_local, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
  
  // ir_local points to colidx idxs, so we need to normalize
  // so that it points to colidx_local
  int ir_offset = ir_local[0];
  for (int i=0; i<bound; i++) { ir_local[i] -= ir_offset; }
  
  MPI_Barrier(MPI_COMM_WORLD);

  // ********** Run FENNEL ***************************************
  //fprintf (stdout, "\n===== Running fennel =====\n");
  //run_mpi_fennel(repr, parts, 1.5, rank); //todo: nparts, gamma as inputs
  // *************************************************************
  //destroy_sparse_matrix (A);
  
  /* END */
  MPI_Finalize();
  return 0;
}