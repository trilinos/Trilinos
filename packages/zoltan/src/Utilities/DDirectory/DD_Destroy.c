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






/**********************  Zoltan_DD_Destroy()  *********************/

void Zoltan_DD_Destroy (
 Zoltan_DD_Directory **dd)     /* contains directory state information */
   {
   int i ;
   DD_Node *ptr ;
   DD_Node *next ;

   int my_proc ;
   int debug_level ;
   char *yo = "ZOLTAN_DD_Destroy" ;

   /* input sanity check */
   if (dd == NULL || *dd == NULL)
      {
      ZOLTAN_PRINT_ERROR (0, yo, "Input argument dd is NULL.") ;
      return ;
      }

   if ((*dd)->debug_level > 1)
      ZOLTAN_TRACE_IN ((*dd)->my_proc, yo, NULL) ;

   /* for each linked list head, walk its list freeing memory */
   for (i = 0 ; i < (*dd)->table_length ; i++)
      for (ptr = (*dd)->table[i] ; ptr != NULL ; ptr = next)
         {
         next = ptr->next ;            /* save before deletion         */
         ZOLTAN_FREE (&ptr) ;              /* destroy node                 */
         }

   /* execute user registered cleanup function, if needed */
   if ((*dd)->cleanup != NULL)
       (*dd)->cleanup() ;

   MPI_Comm_free (&((*dd)->comm)) ;    /* free MPI Comm, ignore errors */

   debug_level = (*dd)->debug_level ;  /* save for final debug print   */
   my_proc     = (*dd)->my_proc ;      /* save for final debug print   */

   ZOLTAN_FREE (dd) ;                      /* free directory structure     */

   if (debug_level > 1)
      ZOLTAN_TRACE_OUT (my_proc, yo, NULL) ;
   }
