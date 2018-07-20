/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "smalloc.h" // for smalloc, sfree
#include "structs.h"

void coarsen1(struct vtx_data ** graph,      /* array of vtx data for graph */
              int                nvtxs,      /* number of vertices in graph */
              int                nedges,     /* number of edges in graph */
              struct vtx_data ***pcgraph,    /* coarsened version of graph */
              int *              pcnvtxs,    /* number of vtxs in coarsened graph */
              int *              pcnedges,   /* number of edges in coarsened graph */
              int **             pv2cv,      /* pointer to v2cv */
              int                igeom,      /* dimension for geometric information */
              float **           coords,     /* coordinates for vertices */
              float **           ccoords,    /* coordinates for coarsened vertices */
              int                using_ewgts /* are edge weights being used? */
              )
{
  extern double coarsen_time;
  extern double match_time;
  double        time;    /* time routine is entered */
  int *         v2cv;    /* maps from vtxs to cvtxs */
  int *         mflag;   /* flag indicating vtx matched or not */
  int           cnvtxs;  /* number of vtxs in coarse graph */
  int           nmerged; /* number of edges contracted */
  double        seconds();
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
