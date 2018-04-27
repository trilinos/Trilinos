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

#include "params.h"  // for MAXSETS
#include "smalloc.h" // for smalloc, sfree

void sorts3d(
    /* Sort the lists needed to find the splitter. */
    double *vals[8][MAXSETS],    /* lists of values to sort */
    int *   indices[8][MAXSETS], /* indices of sorted lists */
    int     nvtxs                /* number of vertices */
)
{
  int *space;       /* space for mergesort routine */
  int  nsets  = 8;  /* number of sets */
  int  nlists = 13; /* number of directions to sort */
  int *temp[13];    /* place holders for indices */
  int  i, j;        /* loop counter */

  void mergesort();

  space = smalloc(nvtxs * sizeof(int));

  for (i = 0; i < nlists; i++) {
    temp[i] = smalloc(nvtxs * sizeof(int));
  }

  mergesort(vals[0][1], nvtxs, temp[0], space);
  mergesort(vals[0][2], nvtxs, temp[1], space);
  mergesort(vals[0][4], nvtxs, temp[2], space);
  mergesort(vals[0][3], nvtxs, temp[3], space);
  mergesort(vals[1][2], nvtxs, temp[4], space);
  mergesort(vals[0][5], nvtxs, temp[5], space);
  mergesort(vals[1][4], nvtxs, temp[6], space);
  mergesort(vals[0][6], nvtxs, temp[7], space);
  mergesort(vals[2][4], nvtxs, temp[8], space);
  mergesort(vals[0][7], nvtxs, temp[9], space);
  mergesort(vals[1][6], nvtxs, temp[10], space);
  mergesort(vals[2][5], nvtxs, temp[11], space);
  mergesort(vals[3][4], nvtxs, temp[12], space);

  sfree(space);

  indices[0][1] = indices[2][3] = indices[4][5] = indices[6][7] = temp[0];
  indices[0][2] = indices[1][3] = indices[4][6] = indices[5][7] = temp[1];
  indices[0][4] = indices[1][5] = indices[2][6] = indices[3][7] = temp[2];
  indices[0][3] = indices[4][7] = temp[3];
  indices[1][2] = indices[5][6] = temp[4];
  indices[0][5] = indices[2][7] = temp[5];
  indices[1][4] = indices[3][6] = temp[6];
  indices[0][6] = indices[1][7] = temp[7];
  indices[2][4] = indices[3][5] = temp[8];
  indices[0][7]                 = temp[9];
  indices[1][6]                 = temp[10];
  indices[2][5]                 = temp[11];
  indices[3][4]                 = temp[12];

  for (i = 0; i < nsets; i++) {
    for (j = i + 1; j < nsets; j++) {
      indices[j][i] = indices[i][j];
    }
  }
}
