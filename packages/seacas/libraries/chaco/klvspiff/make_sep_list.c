/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

int make_sep_list(int *bspace,      /* list of vtxs to be moved */
                  int  list_length, /* current length of bspace */
                  int *sets         /* processor each vertex is assigned to */
)
{
  int vtx;  /* vertex in list */
  int i, k; /* loop counters */

  /* Compress out the actual boundary vertices. */
  k = 0;
  for (i = 0; i < list_length; i++) {
    vtx = bspace[i];
    if (vtx < 0) {
      vtx = -vtx;
    }
    if (sets[vtx] == 2) {
      bspace[k++] = vtx;
    }
  }

  bspace[k] = 0;
  return (k);
}
