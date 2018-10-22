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

#include "structs.h" // for vtx_data

/* Interpolate the eigenvector from the coarsened to the original graph.
   This may require regenerating mflag and v2cv arrays from merged array.

   I also need to undo the merged edges in the reverse order from that in
   which they were collapsed.
*/

void ch_interpolate(double **         vecs,       /* approximate eigenvectors for graph */
                    double **         cvecs,      /* exact eigenvectors for coarse graph */
                    int               ndims,      /* number of vectors to interpolate */
                    struct vtx_data **graph,      /* array of vtx data for graph */
                    int               nvtxs,      /* number of vertices in graph */
                    int *             v2cv,       /* mapping from vtxs to cvtxs */
                    int               using_ewgts /* are edge weights being used in fine graph? */
)
{
  double *vec, *cvec;  /* pointers into vecs and vecs */
  int *   eptr;        /* loops through edge lists */
  float * ewptr;       /* loops through edge weights */
  float   ewgt;        /* value for edge weight */
  double  ewsum;       /* sum of incident edge weights */
  double  sum;         /* sum of values of neighbors */
  int     neighbor;    /* neighboring vertex */
  int     degree;      /* number of neighbors */
  int     npasses = 1; /* number of Gauss-Seidel iterations */
  int     pass;        /* loops through Gauss-Seidel iterations */
  int     i, j, k;     /* loop counters */

  /* Uncompress the coarse eigenvector by replicating matched values. */
  for (j = 1; j <= ndims; j++) {
    vec  = vecs[j];
    cvec = cvecs[j];
    for (i = 1; i <= nvtxs; i++) {
      vec[i] = cvec[v2cv[i]];
    }
  }

  /* Now do a single pass of Gauss-Seidel to smooth eigenvectors. */

  for (pass = 1; pass <= npasses; pass++) {
    if (using_ewgts) {
      for (j = 1; j <= ndims; j++) {
        vec = vecs[j];
        for (i = 1; i <= nvtxs; i++) {
          eptr   = graph[i]->edges;
          ewptr  = graph[i]->ewgts;
          sum    = 0;
          ewsum  = 0;
          degree = graph[i]->nedges - 1;
          for (k = degree; k; k--) {
            neighbor = *(++eptr);
            ewgt     = *(++ewptr);
            sum += ewgt * vec[neighbor];
            ewsum += ewgt;
          }
          vec[i] = sum / ewsum;
        }
      }
    }

    else {
      for (j = 1; j <= ndims; j++) {
        vec = vecs[j];
        for (i = 1; i <= nvtxs; i++) {
          eptr   = graph[i]->edges;
          sum    = 0;
          degree = graph[i]->nedges - 1;
          for (k = degree; k; k--) {
            neighbor = *(++eptr);
            sum += vec[neighbor];
          }
          vec[i] = sum / degree;
        }
      }
    }
  }
}
