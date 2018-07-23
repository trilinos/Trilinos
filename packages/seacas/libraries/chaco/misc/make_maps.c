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

#include <stdio.h>

int make_maps(int *setlists,  /* linked list of set vertices */
              int *list_ptrs, /* head of each linked list */
              int  set,       /* set value denoting subgraph */
              int *glob2loc,  /* graph -> subgraph numbering map */
              int *loc2glob   /* subgraph -> graph numbering map */
)
{
  int i, j; /* loop counter */

  j = 0;
  i = list_ptrs[set];

  if (glob2loc != NULL) {
    while (i != 0) {
      loc2glob[++j] = i;
      glob2loc[i]   = j;
      i             = setlists[i];
    }
  }

  else {
    while (i != 0) {
      loc2glob[++j] = i;
      i             = setlists[i];
    }
  }

  return (j);
}

void make_maps2(int *assignment, /* set assignments for graph */
                int  nvtxs,      /* length of assignment */
                int  set,        /* set value denoting subgraph */
                int *glob2loc,   /* graph -> subgraph numbering map */
                int *loc2glob    /* subgraph -> graph numbering map */
)
{
  int i, j; /* loop counter */

  j = 0;
  if (glob2loc != NULL) {
    for (i = 1; i <= nvtxs; i++) {
      if (assignment[i] == set) {
        j++;
        glob2loc[i] = j;
        loc2glob[j] = i;
      }
    }
  }
  else {
    for (i = 1; i <= nvtxs; i++) {
      if (assignment[i] == set) {
        j++;
        loc2glob[j] = i;
      }
    }
  }
}
