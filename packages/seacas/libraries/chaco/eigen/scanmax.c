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

#include "structs.h" // for scanlink
#include <stdio.h>   // for NULL

/* Return maximum entries of vector over range */
void scanmax(double *vec,               /* vector to scan */
             int beg, int end,          /* index range */
             struct scanlink **scanlist /* pntr to list holding results of scan */
             )
{
  extern double    DOUBLE_MAX;
  struct scanlink *top;
  struct scanlink *curlnk;
  struct scanlink *prevlnk;
  double           val;
  int              i;

  curlnk = *scanlist;
  while (curlnk != NULL) {
    curlnk->indx = 0;
    curlnk->val  = -DOUBLE_MAX;
    curlnk       = curlnk->pntr;
  }

  /* Note: Uses current top link (which would need to be deleted anyway) each time
           an insertion to the list is required. */

  for (i = beg; i <= end; i++) {
    /* consider each element for insertion */
    top = *scanlist;
    val = vec[i];
    if (val > top->val) {
      if (top->pntr == NULL) {
        /* the list is only one long, so just replace */
        top->val  = val;
        top->indx = i;
      }
      else {
        /* beats top element; scan for insertion point */
        if (val > (top->pntr)->val) {
          /* 2nd link becomes list pntr; otherwise stays same */
          *scanlist = top->pntr;
        }
        prevlnk = curlnk = top;
        while ((val > curlnk->val) && (curlnk->pntr != NULL)) {
          prevlnk = curlnk;
          curlnk  = curlnk->pntr;
        }
        if (val > curlnk->val) {
          /* got to end of list; add top to bottom */
          curlnk->pntr = top;
          top->val     = val;
          top->indx    = i;
          top->pntr    = NULL;
        }
        else {
          /* stopped within list; insert top here */
          prevlnk->pntr = top;
          top->val      = val;
          top->indx     = i;
          top->pntr     = curlnk;
        }
      }
    }
  }
}
