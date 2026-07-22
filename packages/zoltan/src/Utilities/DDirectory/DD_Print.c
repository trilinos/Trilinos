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






/********************   Zoltan_DD_Print()  ***********************/

int Zoltan_DD_Print (
 Zoltan_DD_Directory *dd)        /* contains directory state information */
{
   int      i,j;
   DD_NodeIdx nodeidx;
   DD_Node *ptr;
   char    *yo = "Zoltan_DD_Print";


   /* input sanity checks */
   if (dd == NULL) {
      ZOLTAN_PRINT_ERROR (ZOLTAN_DD_NO_PROC, yo, "NULL dd input argument");
      return ZOLTAN_FATAL;
   }
   if (dd->debug_level > 4)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL);

   /* walk linked list printing each node */
   for (i = 0; i < dd->table_length; i++)
      for (nodeidx = dd->table[i]; nodeidx != -1;
           nodeidx = dd->nodelist[nodeidx].next) {
         ptr = dd->nodelist + nodeidx;
         printf ("ZOLTAN DD Print(%d): \tList %3d, \tGID ", dd->my_proc, i);
         printf("(");
         for (j = 0 ; j < dd->gid_length; j++)
            printf(ZOLTAN_ID_SPEC " ", ptr->gid[j]);
         printf(") ");
         if (dd->lid_length > 0)  {
            printf("\tLID (");
            for (j = 0; j < dd->lid_length; j++)
               printf( ZOLTAN_ID_SPEC " ", ptr->gid[j+dd->gid_length]);
            printf(") ");
         }
         printf ("\tPart %d\n", ptr->partition);
         printf ("\tOwner %d\n", ptr->owner);
      }

   if (dd->debug_level > 4)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;
   return ZOLTAN_OK;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
