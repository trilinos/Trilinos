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






/*******************  Zoltan_DD_Stats()  ***************************/

void Zoltan_DD_Stats (Zoltan_DD_Directory *dd)
   {
   int total  = 0 ;
   int maxlen = 0;
   int chains = 0 ;

   int len ;
   int i ;
   DD_Node *ptr ;
   char str[100] ;  /* used to build message string */
   char *yo = "Zoltan_DD_Stats" ;


   /* Input sanity check */
   if (dd == NULL)
      {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument") ;
      return ;
      }

   for (i = 0 ; i < dd->table_length ; i++)
      {
      len = 0 ;
      if (dd->table[i] != NULL)
         chains++ ;

      for (ptr = dd->table[i] ; ptr != NULL ; ptr = ptr->next)
         {
         sprintf (str, "GID %4u, Owner %d, Table Index %d\n",
          *ptr->gid, ptr->owner, i) ;
         ZOLTAN_PRINT_INFO (dd->my_proc, yo, str) ;

         len++ ;
         total++ ;
         }
      if (len > maxlen)
         maxlen = len ;
      }

   if (dd->debug_level >= 0)
      {
      sprintf (str, "Hash table size %d, %d nodes with %d chains, max chain %d",
       dd->table_length, total, chains, maxlen) ;
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, str) ;
      }
   }
