// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include <stdio.h>
#include <stdlib.h>

#include "zoltan_dd_const.h"


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
