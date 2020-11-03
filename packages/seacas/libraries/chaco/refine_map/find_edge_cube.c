/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "refine_map.h" // for refine_edata

struct refine_edata *
find_edge_cube(int                  node,     /* processor node */
               int                  dim,      /* direction of edge from node */
               struct refine_edata *edata,    /* data structure for edge preferences */
               int                  nsets_tot /* total number of processors */
)
{
  struct refine_edata *eguy;  /* returned pointer to edge info */
  int                  index; /* computed index into edata */

  /* Squeeze out bit dim from node number. */
  index = node ^ ((node >> dim) << dim);
  index ^= ((node >> (dim + 1)) << dim);
  index += dim * nsets_tot / 2;

  eguy = &(edata[index]);

  return (eguy);
}
