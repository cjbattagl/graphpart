#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/read_mm.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <bebop/util/config.h>
#include <bebop/util/get_options.h>
#include <bebop/util/init.h>
#include <bebop/util/log.h>
#include <bebop/util/malloc.h>
#include <bebop/util/timer.h>
#include <bebop/util/util.h>

#include "fennel.h"

struct 
cmdlineopts_t
{
  char* input_filename;
  char* input_file_format;
  char* output_filename;
  char* output_file_format;
} opts;

static void
usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args)
{
  fprintf (out, "Usage:\n");
  fprintf (out, "fennel <in-filename> <in-format> <out-filename> <out-format> [options]\n");
  fprintf (out, "<in-filename>:   name of file containing the sparse matrix to read in\n");
  fprintf (out, "<in-format>:     format of the input file (\"HB\" for Harwell-Boeing, \"ML\" for\n");
  fprintf (out, "                 Matlab or \"MM\" for MatrixMarket)\n");
  fprintf (out, "<out-filename>:  name of file to which to output the sparse matrix\n");
  fprintf (out, "<out-format>:    format of the output file (\"HB\" for Harwell-Boeing, \"ML\" for\n");
  fprintf (out, "                 Matlab or \"MM\" for MatrixMarket)\n");
  fprintf (out, "[options]: -a -- validate the input matrix only, without outputting anything\n");
  fprintf (out, "           -e -- expand symmetric into unsymmetric storage (this option is\n");
  fprintf (out, "                 unnecessary for output to Matlab format, as Matlab format is already\n");
  fprintf (out, "                 expanded)\n");
  fprintf (out, "           -v  -- verbose mode\n");
  fprintf (out, "           -h  -- print this usage message and exit\n\n");
  fprintf (out, "EX: ./fennel -a 'test.mtx' 'MM'\n");
  fprintf (out, "EX: ./fennel -v 'test.mtx' 'MM' 'out.mtx' 'MM'\n");
}

int main (int argc, char *argv[])
{
  extern int optind;

  enum sparse_matrix_file_format_t informat = 0;
  enum sparse_matrix_file_format_t outformat = 0;
  struct sparse_matrix_t* A = NULL;
  struct arginfo *arglist = NULL;
  double seconds = 0.0;
  int errcode = 0;

  bebop_default_initialize (argc, argv, &errcode);
  if (errcode != 0)
    {
      fprintf (stderr, "*** Failed to initialize BeBOP Utility Library "
	       "(error code %d) ***\n", errcode);
      bebop_exit (EXIT_FAILURE);
    }

  const char* matrix_filename = argv[1];
  /* Set the get_options usage function to "usage", instead of using the default
   * usage function.  This is necessary because the command-line arguments include
   * things that are not "options" in the strict sense, because they do not follow
   * a "-[a-z]".
   */
  register_usage_function (usage);

  arglist = register_arginfo (arglist, 'v', NULLARG, NULL, "If specified, ac"
			      "tivate verbose mode");
  arglist = register_arginfo (arglist, 'e', NULLARG, NULL, "If specified, ex"
			      "pand the input matrix from symmetric storage "
			      "into unsymmetric storage");
  arglist = register_arginfo (arglist, 'a', NULLARG, NULL, "If specified, va"
			      "lidate the input matrix, without outputting a"
			      "nything");
  get_options (argc, argv, arglist, NULL);

  if (got_arg_p (arglist, 'a'))
    {
      int errcode = do_validate_matrix (argc, argv, arglist);
      destroy_arginfo_list (arglist);
      deinit_timer();
      bebop_exit (errcode); /* stops logging */
    }

  if (argc - optind != 4)
    {
      fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
	       "e specified, but there should be %d ***\n", argc - optind, 4);
      dump_usage (stderr, argv[0], arglist, NULL);
      bebop_exit (EXIT_FAILURE); /* stops logging */
    }

  opts.input_filename = argv[optind];
  opts.input_file_format = argv[optind+1];
  opts.output_filename = argv[optind+2];
  opts.output_file_format = argv[optind+3];

  if (strcmp (opts.input_file_format, "HB") == 0 || 
      strcmp (opts.input_file_format, "hb") == 0)
    informat = HARWELL_BOEING;
  else if (strcmp (opts.input_file_format, "MM") == 0 ||
	   strcmp (opts.input_file_format, "mm") == 0)
    informat = MATRIX_MARKET;
  else if (strcmp (opts.input_file_format, "ML") == 0 ||
           strcmp (opts.input_file_format, "ml") == 0)
    informat = MATLAB;
  else
    {
      fprintf (stderr, "*** Unsupported input file format \"%s\" ***\n", 
	       opts.input_file_format);
      dump_usage (stderr, argv[0], arglist, NULL);
      destroy_arginfo_list (arglist);
      bebop_exit (EXIT_FAILURE); /* stops logging */
    }

  if (strcmp (opts.output_file_format, "HB") == 0 || 
      strcmp (opts.output_file_format, "hb") == 0)
    outformat = HARWELL_BOEING;
  else if (strcmp (opts.output_file_format, "MM") == 0 ||
	   strcmp (opts.output_file_format, "mm") == 0)
    outformat = MATRIX_MARKET;
  else if (strcmp (opts.output_file_format, "ML") == 0 ||
           strcmp (opts.output_file_format, "ml") == 0)
    outformat = MATLAB;
  else
    {
      fprintf (stderr, "*** Unsupported output file format \"%s\" ***\n", 
	       opts.output_file_format);
      dump_usage (stderr, argv[0], arglist, NULL);
      destroy_arginfo_list (arglist);
      bebop_exit (EXIT_FAILURE); /* stops logging */
    }

  if (got_arg_p (arglist, 'v'))
    {
      printf ("Loading sparse matrix...");
      fflush (stdout); /* Otherwise the message may not be printed until 
			  after loading is complete */
    }
  seconds = get_seconds();
  A = load_sparse_matrix (informat, opts.input_filename);
  seconds = get_seconds() - seconds;
  if (A == NULL)
    {
      fprintf (stderr, "*** Failed to load input matrix file \"%s\" ***\n", 
	       opts.input_filename);
      destroy_arginfo_list (arglist);
      bebop_exit (1); /* stops logging */
    }
  if (got_arg_p (arglist, 'v'))
    printf ("done, in %g seconds\n", seconds);

  fprintf (stdout, "\n===== Running fennel =====\n");
  //run_driver (&(A->repr), "rhist__sequential.out");

  destroy_sparse_matrix (A);
  return 0;
}

static void setup_mv (int n, double** p_b, double** p_x)
{
  if (n <= 0) return;

  if (p_b) {
    if (!(*p_b)) { // not yet allocated
      *p_b = (double *)malloc (n * sizeof (double));
      assert (*p_b);
    }
    for (int i = 0; i < n; ++i)
      (*p_b)[i] = 1.0;
  }
  if (p_x) {
    if (!(*p_x)) { // not yet allocated
      *p_x = (double *)malloc (n * sizeof (double));
      assert (*p_x);
    }
  }
}

/**
 * Perform the matrix validation operation specified by the "-a" command-line option.
 */
static int
do_validate_matrix (int argc, char *argv[], struct arginfo* arglist)
{
  extern int optind;

  enum sparse_matrix_file_format_t informat = 0;
  struct sparse_matrix_t* A = NULL;
  double seconds = 0.0;
  int errcode = 0;

  if (argc - optind != 2)
    {
      fprintf (stderr, "*** Incorrect number of command-line arguments: %d ar"
	       "e specified, but there should be %d ***\n", argc - optind, 2);
      dump_usage (stderr, argv[0], arglist, NULL);
      return -1;
    }

  opts.input_filename = argv[optind];
  opts.input_file_format = argv[optind+1];

  if (strcmp (opts.input_file_format, "HB") == 0 || 
      strcmp (opts.input_file_format, "hb") == 0)
    informat = HARWELL_BOEING;
  else if (strcmp (opts.input_file_format, "MM") == 0 ||
	   strcmp (opts.input_file_format, "mm") == 0)
    informat = MATRIX_MARKET;
  else if (strcmp (opts.input_file_format, "ML") == 0 ||
           strcmp (opts.input_file_format, "ml") == 0)
    informat = MATLAB;
  else
    {
      fprintf (stderr, "*** Unsupported input file format \"%s\" ***\n", 
	       opts.input_file_format);
      dump_usage (stderr, argv[0], arglist, NULL);
      return -1;
    }

  if (got_arg_p (arglist, 'v'))
    {
      printf ("Loading sparse matrix...");
      fflush (stdout); /* Otherwise the message may not be printed until 
			  after loading is complete */
    }
  seconds = get_seconds();
  A = load_sparse_matrix (informat, opts.input_filename);
  seconds = get_seconds() - seconds;
  if (A == NULL)
    {
      fprintf (stderr, "*** Failed to load input matrix file \"%s\" ***\n", 
	       opts.input_filename);
      destroy_arginfo_list (arglist);
      return 1;
    }
  if (got_arg_p (arglist, 'v'))
    printf ("done, in %g seconds\n", seconds);


  if (got_arg_p (arglist, 'v'))
    {
      printf ("Validating sparse matrix...");
      fflush (stdout); 
    }
  seconds = get_seconds();
  errcode = valid_sparse_matrix (A);
  seconds = get_seconds() - seconds;
  if (got_arg_p (arglist, 'v'))
    printf ("done, in %g seconds\n", seconds);

  if (valid_sparse_matrix (A))
    printf ("\n### Sparse matrix is valid ###\n\n");
  else 
    printf ("\n### Invalid sparse matrix! ###\n\n");

  destroy_sparse_matrix (A);
  return 0;
}