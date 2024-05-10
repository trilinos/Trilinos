/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for smalloc
#include "structs.h" // for vtx_data

/* Make coarse graph vertex coordinates be center-of-mass of their */
/* fine graph constituents. */

void makeccoords(struct vtx_data **graph,     /* array of vtx data for graph */
                 int               cnvtxs,    /* number of vertices in coarse graph */
                 int              *cv2v_ptrs, /* vtxs corresponding to each cvtx */
                 int              *cv2v_vals, /* indices into cv2v_vals */
                 int               igeom,     /* dimensions of geometric data */
                 float           **coords,    /* coordinates for vertices */
                 float           **ccoords    /* coordinates for coarsened vertices */
)
{
  double mass; /* total mass of merged vertices */
  float *cptr; /* loops through ccoords */
  int    cvtx; /* coarse graph vertex */
  int    vtx;  /* vertex being merged */
  int    i, j; /* loop counters */

  for (i = 0; i < igeom; i++) {
    ccoords[i] = cptr = smalloc((cnvtxs + 1) * sizeof(float));
    for (cvtx = cnvtxs; cvtx; cvtx--) {
      *(++cptr) = 0;
    }
  }
  if (igeom == 1) {
    for (cvtx = 1; cvtx <= cnvtxs; cvtx++) {
      mass = 0;
      for (j = cv2v_ptrs[cvtx]; j < cv2v_ptrs[cvtx + 1]; j++) {
        vtx = cv2v_vals[j];
        mass += graph[vtx]->vwgt;
        ccoords[0][cvtx] += graph[vtx]->vwgt * coords[0][vtx];
      }
      ccoords[0][cvtx] /= mass;
    }
  }
  else if (igeom == 2) {
    for (cvtx = 1; cvtx <= cnvtxs; cvtx++) {
      mass = 0;
      for (j = cv2v_ptrs[cvtx]; j < cv2v_ptrs[cvtx + 1]; j++) {
        vtx = cv2v_vals[j];
        mass += graph[vtx]->vwgt;
        ccoords[0][cvtx] += graph[vtx]->vwgt * coords[0][vtx];
        ccoords[1][cvtx] += graph[vtx]->vwgt * coords[1][vtx];
      }
      ccoords[0][cvtx] /= mass;
      ccoords[1][cvtx] /= mass;
    }
  }
  else if (igeom > 2) {
    for (cvtx = 1; cvtx <= cnvtxs; cvtx++) {
      mass = 0;
      for (j = cv2v_ptrs[cvtx]; j < cv2v_ptrs[cvtx + 1]; j++) {
        vtx = cv2v_vals[j];
        mass += graph[vtx]->vwgt;
        ccoords[0][cvtx] += graph[vtx]->vwgt * coords[0][vtx];
        ccoords[1][cvtx] += graph[vtx]->vwgt * coords[1][vtx];
        ccoords[2][cvtx] += graph[vtx]->vwgt * coords[2][vtx];
      }
      ccoords[0][cvtx] /= mass;
      ccoords[1][cvtx] /= mass;
      ccoords[2][cvtx] /= mass;
    }
  }
}
