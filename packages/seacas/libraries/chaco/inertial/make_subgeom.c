/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

void make_subgeom(int     igeom,     /* 1, 2 or 3 dimensional geometry? */
                  float **coords,    /* x, y and z coordinates of vertices */
                  float **subcoords, /* x, y and z coordinates in subgraph */
                  int     subnvtxs,  /* number of vertices in subgraph */
                  int    *loc2glob   /* maps from subgraph to graph numbering */
)
{
  int i; /* loop counter */

  if (igeom == 1) {
    for (i = 1; i <= subnvtxs; i++) {
      subcoords[0][i] = coords[0][loc2glob[i]];
    }
  }
  else if (igeom == 2) {
    for (i = 1; i <= subnvtxs; i++) {
      subcoords[0][i] = coords[0][loc2glob[i]];
      subcoords[1][i] = coords[1][loc2glob[i]];
    }
  }
  else if (igeom > 2) {
    for (i = 1; i <= subnvtxs; i++) {
      subcoords[0][i] = coords[0][loc2glob[i]];
      subcoords[1][i] = coords[1][loc2glob[i]];
      subcoords[2][i] = coords[2][loc2glob[i]];
    }
  }
}
