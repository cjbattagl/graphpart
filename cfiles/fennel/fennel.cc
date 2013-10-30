#include "fennel.h"
#include <bebop/util/config.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/bcoo_matrix.h>
#include <bebop/smc/bcsr_matrix.h>
#include <bebop/smc/coo_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/jad_matrix.h>

int
main (int argc, char *argv[])
{
  /* Process arguments */
  if (argc != 2) {
    fprintf (stderr, "usage: %s <input-matrix>\n", argv[0]);
    return -1;
  }

}
