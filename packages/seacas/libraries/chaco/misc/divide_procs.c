/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"    // for min
#include "params.h"  // for MAXSETS
#include "structs.h" // for set_info

int divide_procs(int              architecture, /* 0 => hypercube, d => d-dimensional mesh */
                 int              ndims,        /* normal dimension of each cut */
                 int              ndims_tot,    /* total number of hypercube dimensions */
                 struct set_info *info_set,     /* data for all sets */
                 struct set_info *divide_set,   /* data for set being divided */
                 int *            subsets,      /* subsets to be created */
                 int              inert,        /* using inertial method? */
                 int *            pndims_real,  /* actual ndims for this cut */
                 int *            pnsets_real,  /* # sets created by this cut */
                 int *            pstriping,    /* cut in single direction? */
                 int *            cut_dirs,     /* direction of each cut if mesh */
                 int *            mesh_dims,    /* size of full mesh */
                 int              hops_special[][MAXSETS] /* hop matrix for nonstandard cases */
)
{
  int nsets_real = -1; /* number of sets to divide into */
  int ndims_real = -1; /* number of eigenvectors to use */
  int striping   = -1; /* cut in single direction? */
  int flag       = -1; /* unusual partition => use special hops */
  int ndim_poss;       /* largest dimensionality possible */
  int idims;           /* true dimensionality of subgrid */
  int i;               /* loop counter */
  int define_submeshes(), define_subcubes();

  if (architecture > 0) { /* Mesh, complicated case. */
    nsets_real = divide_set->span[0] * divide_set->span[1] * divide_set->span[2];
    nsets_real = min(1 << ndims, nsets_real);
    ndims_real = ndims;
    while (1 << ndims_real > nsets_real) {
      --ndims_real;
    }

    ndim_poss = 0;
    idims     = 0;
    for (i = 0; i < 3; i++) {
      if (divide_set->span[i] >= 2) {
        ndim_poss++;
        idims++;
      }
      if (divide_set->span[i] >= 4) {
        ndim_poss++;
      }
      if (divide_set->span[i] >= 8) {
        ndim_poss++;
      }
    }
    ndims_real = min(ndim_poss, ndims_real);

    if (idims > 1) {
      nsets_real = 1 << ndims_real;
    }

    flag = define_submeshes(nsets_real, architecture, mesh_dims, divide_set, info_set, subsets,
                            inert, &striping, cut_dirs, hops_special);
    if (striping) {
      ndims_real = 1;
    }
  }

  else if (architecture == 0) { /* Hypercube, easy case. */
    ndims_real = min(ndims, divide_set->ndims);
    nsets_real = 1 << ndims_real;

    flag = define_subcubes(nsets_real, ndims_tot, ndims_real, divide_set, info_set, subsets, inert,
                           &striping, hops_special);

    if (striping) {
      ndims_real = 1;
    }
  }

  *pndims_real = ndims_real;
  *pnsets_real = nsets_real;
  *pstriping   = striping;

  return (flag);
}
