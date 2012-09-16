/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*  NOTE: See file, README, for associated documentation. (RTH) */






/*******************  Zoltan_DD_Stats()  ***************************/

void Zoltan_DD_Stats (
 Zoltan_DD_Directory *dd)   /* directory state information */
{
   int node_count = 0;     /* counts Nodes in local directory      */
   int maxlength  = 0;     /* length of longest linked list        */
   int list_count = 0;     /* number of linked lints in hash table */

   int      length;
   int      i;
   DD_NodeIdx nodeidx;
   DD_Node *ptr;
   char     str[100];      /* used to build message string */
   char    *yo = "Zoltan_DD_Stats";


   /* Input sanity check */
   if (dd == NULL) {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument.");
      return;
   }

   if (dd->debug_level > 4)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL);

   /* walk down each list in hash table to find every Node */
   for (i = 0; i < dd->table_length; i++) {
      length = 0;                    /* reset length for next count */
      if (dd->table[i] != -1)
         list_count++;               /* count of distict linked lists */

      for (nodeidx = dd->table[i]; nodeidx != -1;
           nodeidx = dd->nodelist[nodeidx].next) {
         ptr = dd->nodelist + nodeidx;
         if (dd->debug_level > 6) {
            sprintf(str, "GID " ZOLTAN_ID_SPEC ", Owner %d, Table Index %d.",
                    *ptr->gid, ptr->owner, i);
            ZOLTAN_PRINT_INFO (dd->my_proc, yo, str);
         }
         length++;                  /* linked list length */
         node_count++;              /* count of Nodes */
      }
      if (length > maxlength)
         maxlength = length;        /* save length of longest linked list */
   }

   sprintf(str, "Hash table size %d, %d nodes on %d lists, max list length %d.",
           dd->table_length, node_count, list_count, maxlength);
   ZOLTAN_PRINT_INFO (dd->my_proc, yo, str);

   if (dd->debug_level > 4)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
