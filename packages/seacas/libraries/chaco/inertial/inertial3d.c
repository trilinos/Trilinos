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

#include "defs.h"    // for TRUE
#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for vtx_data
#include <stdio.h>   // for printf

void inertial3d(struct vtx_data **graph,        /* graph data structure */
                int               nvtxs,        /* number of vtxs in graph */
                int               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
                int               nsets,        /* number of sets to divide into */
                float *x, float *y, float *z,   /* x, y and z coordinates of vertices */
                int *   sets,                   /* set each vertex gets assigned to */
                double *goal,                   /* desired set sizes */
                int     using_vwgts             /* are vertex weights being used? */
)
{
  extern int    DEBUG_INERTIAL;     /* debug flag for inertial method */
  extern double inertial_axis_time; /* time spent computing inertial axis */
  extern double median_time;        /* time spent finding medians */
  double        tensor[3][3];       /* inertia tensor */
  double        evec[3];            /* eigenvector */
  double *      value;              /* values along selected direction to sort */
  double        xcm, ycm, zcm;      /* center of mass in each direction */
  double        xx, yy, zz;         /* elements of inertial tensor */
  double        xy, xz, yz;         /* elements of inertial tensor */
  double        xdif, ydif;         /* deviation from center of mass */
  double        zdif;               /* deviation from center of mass */
  double        eval, res;          /* eigenvalue and error in eval calculation */
  double        vwgt_sum;           /* sum of all the vertex weights */
  double        time;               /* timing parameter */
  int *         space;              /* space required by median routine */
  int           i;                  /* loop counter */
  double        seconds();

  void ch_eigenvec3(), ch_evals3(), rec_median_1();

  /* Compute center of mass and total mass. */
  time = seconds();
  xcm = ycm = zcm = 0.0;
  if (using_vwgts) {
    vwgt_sum = 0;
    for (i = 1; i <= nvtxs; i++) {
      vwgt_sum += graph[i]->vwgt;
      xcm += graph[i]->vwgt * x[i];
      ycm += graph[i]->vwgt * y[i];
      zcm += graph[i]->vwgt * z[i];
    }
  }
  else {
    vwgt_sum = nvtxs;
    for (i = 1; i <= nvtxs; i++) {
      xcm += x[i];
      ycm += y[i];
      zcm += z[i];
    }
  }

  xcm /= vwgt_sum;
  ycm /= vwgt_sum;
  zcm /= vwgt_sum;

  /* Generate 6 elements of Inertial tensor. */
  xx = yy = zz = xy = xz = yz = 0.0;
  if (using_vwgts) {
    for (i = 1; i <= nvtxs; i++) {
      xdif = x[i] - xcm;
      ydif = y[i] - ycm;
      zdif = z[i] - zcm;
      xx += graph[i]->vwgt * xdif * xdif;
      yy += graph[i]->vwgt * ydif * ydif;
      zz += graph[i]->vwgt * zdif * zdif;
      xy += graph[i]->vwgt * xdif * ydif;
      xz += graph[i]->vwgt * xdif * zdif;
      yz += graph[i]->vwgt * ydif * zdif;
    }
  }
  else {
    for (i = 1; i <= nvtxs; i++) {
      xdif = x[i] - xcm;
      ydif = y[i] - ycm;
      zdif = z[i] - zcm;
      xx += xdif * xdif;
      yy += ydif * ydif;
      zz += zdif * zdif;
      xy += xdif * ydif;
      xz += xdif * zdif;
      yz += ydif * zdif;
    }
  }

  /* Compute eigenvector with maximum eigenvalue. */

  tensor[0][0] = xx;
  tensor[1][1] = yy;
  tensor[2][2] = zz;
  tensor[0][1] = tensor[1][0] = xy;
  tensor[0][2] = tensor[2][0] = xz;
  tensor[1][2] = tensor[2][1] = yz;
  ch_evals3(tensor, &res, &res, &eval);
  ch_eigenvec3(tensor, eval, evec, &res);

  inertial_axis_time += seconds() - time;

  if (DEBUG_INERTIAL > 0) {
    printf("Principle Axis = (%g, %g, %g)\n  Eval=%g, Residual=%e\n", evec[0], evec[1], evec[2],
           eval, res);
  }

  /* Allocate space for value array. */

  value = smalloc((nvtxs + 1) * sizeof(double));

  /* Calculate value to sort/split on for each cell. */
  /* This is inner product with eigenvector. */
  for (i = 1; i <= nvtxs; i++) {
    value[i] = (x[i] - xcm) * evec[0] + (y[i] - ycm) * evec[1] + (z[i] - zcm) * evec[2];
  }

  /* Now find the median value and partition based upon it. */
  space = smalloc(nvtxs * sizeof(int));
  time  = seconds();
  rec_median_1(graph, value, nvtxs, space, cube_or_mesh, nsets, goal, using_vwgts, sets, TRUE);
  median_time += seconds() - time;

  sfree(space);
  sfree(value);
}
