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

#include "DD_Const.h"



/*  NOTE: See file, README, for associated documentation. (RTH) */






/**********************  Zoltan_DD_Destroy()  *********************/

void Zoltan_DD_Destroy (Zoltan_DD_Directory *dd)
   {
   int i ;
   DD_Node *ptr ;
   DD_Node *next ;

   int my_proc ;
   int debug_level ;

   /* input sanity check */
   if (dd == NULL)
      return ;

   /* for each linked list head, walk its list freeing memory */
   for (i = 0 ; i < dd->table_length ; i++)
      for (ptr = dd->table[i] ; ptr != NULL ; ptr = next)
         {
         next = ptr->next ;         /* save before deletion         */
         LB_FREE (&ptr) ;           /* destroy node                 */
         }
   debug_level = dd->debug_level ;  /* to control final debug print */
   my_proc     = dd->my_proc ;

   if (dd->cleanup != NULL)       
      dd->cleanup() ;    /* execute user registered cleanup function */

   MPI_Comm_free (&(dd->comm)) ;    /* free MPI Comm, ignore errors */
   LB_FREE (&dd) ;                  /* free directory structure     */

   if (debug_level> 0)
      printf ("ZOLTAN_DD_DESTROY(%d): Successful\n", my_proc) ;
   }
