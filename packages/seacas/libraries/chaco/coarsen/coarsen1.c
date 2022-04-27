/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for smalloc, sfree
#include "structs.h"

void coarsen1(struct vtx_data  **graph,      /* array of vtx data for graph */
              int                nvtxs,      /* number of vertices in graph */
              int                nedges,     /* number of edges in graph */
              struct vtx_data ***pcgraph,    /* coarsened version of graph */
              int               *pcnvtxs,    /* number of vtxs in coarsened graph */
              int               *pcnedges,   /* number of edges in coarsened graph */
              int              **pv2cv,      /* pointer to v2cv */
              int                igeom,      /* dimension for geometric information */
              float            **coords,     /* coordinates for vertices */
              float            **ccoords,    /* coordinates for coarsened vertices */
              int                using_ewgts /* are edge weights being used? */
)
{
  extern double coarsen_time;
  extern double match_time;
  double        time;    /* time routine is entered */
  int          *v2cv;    /* maps from vtxs to cvtxs */
  int          *mflag;   /* flag indicating vtx matched or not */
  int           cnvtxs;  /* number of vtxs in coarse graph */
  int           nmerged; /* number of edges contracted */
  double        seconds(void);
  int           maxmatch();
  void          makev2cv(), makefgraph();

  time = seconds();

  /* Allocate and initialize space. */
  v2cv  = smalloc((nvtxs + 1) * sizeof(int));
  mflag = smalloc((nvtxs + 1) * sizeof(int));

  /* Find a maximal matching in the graph. */
  nmerged = maxmatch(graph, nvtxs, nedges, mflag, using_ewgts, igeom, coords);
  match_time += seconds() - time;

  /* Now construct coarser graph by contracting along matching edges. */
  /* Pairs of values in mflag array indicate matched vertices. */
  /* A zero value indicates that vertex is unmatched. */

  /*
      makecgraph(graph, nvtxs, pcgraph, pcnvtxs, pcnedges, mflag,
                    *pv2cv, nmerged, using_ewgts, igeom, coords, ccoords);
      makecgraph2(graph, nvtxs, nedges, pcgraph, pcnvtxs, pcnedges, mflag,
                    *pv2cv, nmerged, using_ewgts, igeom, coords, ccoords);
  */

  makev2cv(mflag, nvtxs, v2cv);

  sfree(mflag);

  cnvtxs = nvtxs - nmerged;
  makefgraph(graph, nvtxs, nedges, pcgraph, cnvtxs, pcnedges, v2cv, using_ewgts, igeom, coords,
             ccoords);

  *pcnvtxs = cnvtxs;
  *pv2cv   = v2cv;
  coarsen_time += seconds() - time;
}
