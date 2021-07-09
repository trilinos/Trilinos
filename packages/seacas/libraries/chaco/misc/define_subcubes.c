/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"    // for FALSE, TRUE
#include "params.h"  // for MAXSETS
#include "structs.h" // for set_info

int define_subcubes(int              nsets_real, /* actual number of sets being created */
                    int              ndims_tot,  /* total hypercube dimensions */
                    int              ndims,      /* # dimension in this cut */
                    struct set_info *set,        /* data for set being divided */
                    struct set_info *set_info,   /* data for all sets */
                    int *            subsets,    /* subsets to be created */
                    int              inert,      /* using inertial method? */
                    int *            pstriping,  /* cut in single direction? */
                    int              hop_mtx_special[MAXSETS][MAXSETS] /* nonstandard hop values */
)
{
  extern int KL_METRIC; /* 2 => using hops so generate hop matrix */
  int        hop_flag;  /* use special hop matrix? */
  int        nsets;     /* number of sets being created */
  int        setnum;    /* global number of subset */
  int        bits;      /* number of bits in which two sets differ */
  int        i, j, k;   /* loop counters */
  int        gray();

  nsets    = 1 << ndims;
  hop_flag = FALSE;

  for (k = nsets - 1; k >= 0; k--) { /* Backwards to not overwrite current set. */

    setnum                 = set->setnum | (k << (ndims_tot - set->ndims));
    set_info[setnum].ndims = set->ndims - ndims;
    subsets[k]             = setnum;
  }

  *pstriping = (inert && nsets_real > 2);

  if (*pstriping) { /* Gray code for better mapping. */
    for (k = 0; k < nsets; k++) {
      subsets[k] = gray(subsets[k]);
    }

    if (KL_METRIC == 2) {
      hop_flag = TRUE;
      for (i = 0; i < nsets; i++) {
        hop_mtx_special[i][i] = 0;
        for (j = 0; j < i; j++) {
          hop_mtx_special[i][j] = 0;
          bits                  = (subsets[i]) ^ (subsets[j]);
          while (bits) {
            if (bits & 1) {
              ++hop_mtx_special[i][j];
            }
            bits >>= 1;
          }
          hop_mtx_special[j][i] = hop_mtx_special[i][j];
        }
      }
    }
  }

  return (hop_flag);
}
