/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "DD.h"



/*  NOTE: See file, README, for associated documentation. (RTH) */






/*******************  Zoltan_DD_Stats()  ***************************/

void Zoltan_DD_Stats (
 Zoltan_DD_Directory *dd)   /* directory state information */
   {
   int node_count = 0 ;     /* counts Nodes in local directory      */
   int maxlength  = 0;      /* length of longest linked list        */
   int list_count = 0 ;     /* number of linked lints in hash table */

   int      length ;
   int      i ;
   DD_Node *ptr ;
   char     str[100] ;      /* used to build message string */
   char    *yo = "Zoltan_DD_Stats" ;


   /* Input sanity check */
   if (dd == NULL)
      {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument.") ;
      return ;
      }

   if (dd->debug_level > 1)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL) ;

   /* walk down each list in hash table to find every Node */
   for (i = 0 ; i < dd->table_length ; i++)
      {
      length = 0 ;                    /* reset length for next count */
      if (dd->table[i] != NULL)
         list_count++ ;               /* count of distict linked lists */

      for (ptr = dd->table[i] ; ptr != NULL ; ptr = ptr->next)
         {
         if (dd->debug_level > 1)
            {
            sprintf (str, "GID %4u, Owner %d, Table Index %d.", *ptr->gid,
             ptr->owner, i) ;
            ZOLTAN_PRINT_INFO (dd->my_proc, yo, str) ;
            }
         length++ ;                  /* linked list length */
         node_count++ ;              /* count of Nodes */
         }
      if (length > maxlength)
         maxlength = length ;        /* save length of longest linked list */
      }

   sprintf (str, "Hash table size %d, %d nodes on %d lists, max list length %d.",
    dd->table_length, node_count, list_count, maxlength) ;
   ZOLTAN_PRINT_INFO (dd->my_proc, yo, str) ;

   if (dd->debug_level > 1)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;
   }
