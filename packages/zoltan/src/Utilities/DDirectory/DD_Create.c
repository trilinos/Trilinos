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
#include "zoltan_align.h"



/*  NOTE: See file, README, for associated documentation. (RTH) */




/*******************  Zoltan_DD_Create()  ***************************/

int Zoltan_DD_Create (
 Zoltan_DD_Directory **dd,    /* contains directory state and pointers */
 MPI_Comm comm,               /* Dup'ed and saved for future use       */
 int num_gid,                 /* Number of entries in a global ID.     */
 int num_lid,                 /* Number of entries in a local ID.      
                                 If zero, ignore LIDs                  */
 int user_length,             /* Optional user data length, 0 ignore   */
 int table_length,            /* sizeof hash table, use default if 0   */
 int debug_level)             /* control actions to errors, normally 0 */
   {
   int size ;
   int my_proc;
   char *yo = "Zoltan_DD_Create" ;

   if (MPI_Comm_rank(comm, &my_proc) != MPI_SUCCESS) 
      {
      ZOLTAN_PRINT_ERROR (-1, yo, "MPI_Comm_rank failed.");
      return ZOLTAN_DD_MPI_ERROR;
      }

   if (debug_level > 1)
      ZOLTAN_TRACE_IN (my_proc, yo, NULL);

   /* input sanity check */
   if (dd == NULL || num_gid < 1 || table_length < 0 || num_lid < 0)
      {
      ZOLTAN_PRINT_ERROR (my_proc, yo, "Invalid input argument.") ;
      if (debug_level > 1)
         ZOLTAN_TRACE_OUT (my_proc, yo, NULL);
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   /* malloc memory for the directory structure + hash table */
   size = (table_length == 0) ? ZOLTAN_DD_HASH_TABLE_COUNT : table_length ;
   *dd  = (Zoltan_DD_Directory *) ZOLTAN_MALLOC (sizeof (Zoltan_DD_Directory)
        + size * sizeof (DD_Node*)) ;
   if (*dd == NULL)
      {
      ZOLTAN_PRINT_ERROR (my_proc, yo, "Can not malloc hash table.") ;
      if (debug_level > 1)
        ZOLTAN_TRACE_OUT(my_proc, yo, NULL);
      return ZOLTAN_DD_MEMORY_ERROR ;
      }

   /* NULL heads of link list in hash table */
   memset ((char *) (*dd)->table, '\0', size * sizeof (DD_Node*)) ;

   /* save useful constants into directory for convenience */
   (*dd)->debug_level      = debug_level ;  /* [0,3], default 0          */
   (*dd)->gid_length       = num_gid ;      /* saved input Num_GID       */
   (*dd)->lid_length       = num_lid ;      /* saved input Num_LIB       */
   (*dd)->table_length     = size ;         /* # of linked list heads    */
   (*dd)->user_data_length = user_length ;  /* optional user data length */
   (*dd)->hash             = Zoltan_DD_Hash2;/* default hash algorithm   */
   (*dd)->cleanup          = NULL ;         /* user registered cleanup   */
   (*dd)->max_id_length    = (num_gid > num_lid) ? num_gid : num_lid ;

   /* frequently used dynamic allocation computed sizes */
   size = (num_gid + num_lid + user_length) * sizeof (ZOLTAN_ID_PTR) ;
   (*dd)->node_size       = size + sizeof(DD_Node) ;
   (*dd)->update_msg_size = size + sizeof(DD_Update_Msg) ;

   size = num_gid * sizeof (ZOLTAN_ID_PTR) ;
   (*dd)->remove_msg_size = size + sizeof(DD_Remove_Msg) ;

   size = (user_length + (*dd)->max_id_length) * sizeof (ZOLTAN_ID_PTR) ;
   (*dd)->find_msg_size   = size + sizeof (DD_Find_Msg) ;

   /* force alignment */
   (*dd)->update_msg_size = Zoltan_Align((*dd)->update_msg_size);
   (*dd)->remove_msg_size = Zoltan_Align((*dd)->remove_msg_size);
   (*dd)->find_msg_size   = Zoltan_Align((*dd)->find_msg_size);

   /* duplicate MPI comm to prevent future comm changes from disrupting  */
   /* directory communications & save the associated comm size & rank    */
   if (MPI_Comm_dup  (comm,  &((*dd)->comm))    != MPI_SUCCESS
    || MPI_Comm_size (comm,  &((*dd)->nproc))   != MPI_SUCCESS
    || MPI_Comm_rank (comm,  &((*dd)->my_proc)) != MPI_SUCCESS)
         {
         ZOLTAN_PRINT_ERROR (my_proc, yo, "MPI Problem, unable to continue.") ;
         if (debug_level > 1)
           ZOLTAN_TRACE_OUT(my_proc, yo, NULL);
         return ZOLTAN_DD_MPI_ERROR ;
         }

   if (debug_level > 1)
      ZOLTAN_TRACE_OUT (my_proc, yo, NULL);

   return ZOLTAN_DD_NORMAL_RETURN ;
   }
