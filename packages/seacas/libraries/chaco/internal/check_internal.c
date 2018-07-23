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

#include "defs.h"
#include "internal.h"
#include "structs.h"
#include <stdio.h>

void check_internal(struct vtx_data **graph,      /* graph data structure */
                    int               nvtxs,      /* number of vertices in graph */
                    struct bidint *   int_list,   /* sorted list of internal weights per set */
                    struct bidint *   set_list,   /* head of vtx_elems lists */
                    struct bidint *   vtx_elems,  /* start of vertex data */
                    int *             total_vwgt, /* total weight in each set */
                    int *             assign,     /* current assignment */
                    int               nsets_tot   /* total number of sets */
)
{
  struct bidint *ptr, *ptr2;         /* elements in int_list */
  struct bidint *old_ptr, *old_ptr2; /* elements in set_list */
  int            vwgt_sum;           /* sum of vertex weights */
  int            set, set2;          /* sets two vertices are in */
  int            sum;                /* sum of internal weights */
  int            nseen;              /* number of vertices found in set_lists */
  int            old_val, val;       /* consecutive values in int_list */
  int            vtx;                /* vertex in set list */
  int            internal;           /* is a vertex internal or not? */
  int            size;               /* array spacing */
  int            j, k;               /* loop counters */

  k       = 0;
  size    = (int)(&(int_list[1]) - &(int_list[0]));
  nseen   = 0;
  old_val = -1;
  old_ptr = &(int_list[nsets_tot]);
  for (ptr = int_list[nsets_tot].next; ptr != NULL; ptr = ptr->next) {
    set = ((int)(ptr - int_list)) / size;
    val = ptr->val;
    if (val < old_val) {
      printf("int_list out of order, k=%d, set = %d, old_val=%d, val = %d\n", k, set, old_val, val);
    }
    if (ptr->prev != old_ptr) {
      printf(" int_list back link screwed up, set=%d, k=%d, old_ptr=%ld, ptr->prev = %ld\n", set, k,
             (long)old_ptr, (long)ptr->prev);
    }
    old_ptr = ptr;
    old_val = val;

    vwgt_sum = 0;
    sum      = 0;
    old_ptr2 = &(set_list[set]);
    for (ptr2 = set_list[set].next; ptr2 != NULL; ptr2 = ptr2->next) {
      vtx = ((int)(ptr2 - vtx_elems)) / size;
      vwgt_sum += graph[vtx]->vwgt;
      if (ptr2->prev != old_ptr2) {
        printf(" set_list back link screwed up, set=%d, k=%d, old_ptr2=%ld, ptr2->prev = %ld\n",
               set, k, (long)old_ptr2, (long)ptr2->prev);
      }
      old_ptr2 = ptr2;

      ++nseen;
      if (assign[vtx] != set) {
        printf("assign[%d] = %d, but in set_list[%d]\n", vtx, assign[vtx], set);
      }
      internal = TRUE;
      for (j = 1; j < graph[vtx]->nedges && internal; j++) {
        set2     = assign[graph[vtx]->edges[j]];
        internal = (set2 == set);
      }
      if (internal) {
        sum += graph[vtx]->vwgt;
      }
    }
    if (sum != val) {
      printf("set = %d, val = %d, but I compute internal = %d\n", set, val, sum);
    }
    if (vwgt_sum != total_vwgt[set]) {
      printf(" vwgt_sum = %d, but total_vwgt[%d] = %d\n", vwgt_sum, set, total_vwgt[set]);
    }
    k++;
  }
  if (k != nsets_tot) {
    printf(" Only %d sets in int_sets list, but nsets_tot = %d\n", k, nsets_tot);
  }
  if (nseen != nvtxs) {
    printf(" Only %d vertices found in int_sets lists, but nvtxs = %d\n", nseen, nvtxs);
  }
}
