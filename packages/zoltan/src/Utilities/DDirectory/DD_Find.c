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



static int DD_Find_Local (Zoltan_DD_Directory *dd, LB_ID_PTR gid,
 LB_ID_PTR lid, LB_ID_PTR user, int *partition, int *owner) ;



/********************  Zoltan_DD_Find()  **************************/

int Zoltan_DD_Find (Zoltan_DD_Directory *dd,
 LB_ID_PTR gid,       /* Incoming list of GIDs to get owners proc      */
 LB_ID_PTR lid,       /* Outgoing corresponding list of LIDs           */
 LB_ID_PTR data,      /* Outgoing optional corresponding user data     */
 int *partition,
 int  count,          /* Count of GIDs in above list                   */
 int *owner)          /* Outgoing corresponding list of data locations */
   {
   struct Comm_Obj *plan  = NULL ;  /* efficient MPI communication     */
   char            *rbuff = NULL ;  /* receive buffer                  */
   char            *sbuff = NULL ;  /* send buffer                     */
   int             *procs = NULL ;  /* list of processors to contact   */
   DD_Find_Msg     *ptr   = NULL ;
   int              i ;
   int              nrec ;          /* number of messages to receive   */
   int              err ;           /* return error condition          */
   int              errcount ;
   char             str[100] ;      /* diagnostic message string */
   char            *yo = "Zoltan_DD_Find" ;

   if (dd != NULL && dd->debug_level > 1)
      ZOLTAN_TRACE_ENTER(dd->my_proc, yo, NULL);

   /* input sanity check */
   if (dd == NULL || count < 0 || ((owner == NULL || gid == NULL) && count > 0))
      {
      ZOLTAN_PRINT_ERROR ((dd == NULL ? -1 : dd->my_proc), yo, 
                          "invalid input argument") ;
      if (dd != NULL && dd->debug_level > 1)
         ZOLTAN_TRACE_EXIT((dd == NULL ? -1 : dd->my_proc), yo, NULL);
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   /* allocate memory for processors to contact for directory info */
   procs = (int *) LB_MALLOC (sizeof (int) * ((count == 0) ? 1 : count)) ;
   if (procs == NULL)
      {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "unable to malloc proc list") ;
      if (dd->debug_level > 1)
         ZOLTAN_TRACE_EXIT(dd->my_proc, yo, NULL);
      return ZOLTAN_DD_MEMORY_ERROR ;
      }

   /* allocate total memory for DD_Find_Msg send buffer */
   sbuff = (char *) LB_MALLOC (dd->find_msg_size * ((count == 0) ? 1 : count)) ;
   if (sbuff == NULL)
      {
      LB_FREE (&procs) ;
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "unable to malloc send buffer") ;
      if (dd->debug_level > 1)
         ZOLTAN_TRACE_EXIT(dd->my_proc, yo, NULL);
      return ZOLTAN_DD_MEMORY_ERROR ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After mallocs.");

   /* for each GID, fill DD_Find_Msg buffer and contact list */
   for (i = 0 ; i < count ; i++)
      {
      procs[i] = dd->hash (gid + i*dd->gid_length, dd->gid_length, dd->nproc) ;

      ptr = (DD_Find_Msg *) (sbuff + i*dd->find_msg_size) ;
      ptr->index = i ;
      ptr->proc  = procs[i] ;
      LB_SET_ID (dd->gid_length, ptr->id, gid + i*dd->gid_length) ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After fill.");

   /* if count is zero, create a dummy message to force participation */
   if (count == 0)
      {
      procs[0] = dd->my_proc ;                       /* send to self */
      ptr = (DD_Find_Msg *) sbuff ;
      ptr->proc = ZOLTAN_DD_NO_PROC ;
      count = 1 ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After dummy.");

   /* create efficient communication plan */
   err = LB_Comm_Create (&plan, count, procs, dd->comm,
    ZOLTAN_DD_FIND_MSG_TAG, &nrec) ;
   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Create.");

   if (err != COMM_OK) goto fini ;

   /* allocate total receive buffer */
   rbuff = (char *) LB_MALLOC (nrec * dd->find_msg_size) ;
   if (rbuff == NULL)
      {
      err = ZOLTAN_DD_MEMORY_ERROR ;
      goto fini ;
      }

   /* send out find messages across entire system */
   err = LB_Comm_Do (plan, ZOLTAN_DD_FIND_MSG_TAG+1, sbuff,
    dd->find_msg_size, rbuff) ;
   if (err != COMM_OK)
      goto fini ;

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Do.");

   /* get find messages directed to me, fill in return information */
   errcount = 0 ;
   for (i = 0 ; i < nrec ; i++)
      {
      ptr = (DD_Find_Msg *) (rbuff + i*dd->find_msg_size) ;

      if (ptr->proc == ZOLTAN_DD_NO_PROC)
         continue ;                 /* ignore dummy message */

      err = DD_Find_Local (dd, ptr->id, ptr->id, ptr->id
       + dd->max_id_length, &ptr->partition, &ptr->proc) ;
      if (err == ZOLTAN_DD_GID_NOT_FOUND_ERROR)
         errcount++ ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After fill in return.");

   /* send return information back to requestor */
   err = LB_Comm_Do_Reverse(plan, ZOLTAN_DD_FIND_MSG_TAG+2, rbuff,
    dd->find_msg_size, NULL, sbuff) ;
   if (err != COMM_OK)
      goto fini ;

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Reverse.");

   /* fill in user supplied lists with returned information */
   for (i = 0 ; i < count ; i++)
      {
      ptr = (DD_Find_Msg *) (sbuff + i*dd->find_msg_size) ;

      if (ptr->proc == ZOLTAN_DD_NO_PROC)
         continue ;         /* ignore dummy self message */

      owner[ptr->index] = ptr->proc ;

      if (lid != NULL)
         LB_SET_ID (dd->lid_length, lid + ptr->index * dd->lid_length,
           ptr->id) ;

      if (data != NULL)
         LB_SET_ID (dd->user_data_length, data + ptr->index
          * dd->user_data_length, ptr->id + dd->max_id_length) ;

      if (partition != NULL)
         partition[ptr->index] = ptr->partition ;
      }

   /* if at least one GID was not found, notify caller of error */
   err = (errcount == 0) ? ZOLTAN_DD_NORMAL_RETURN
                         : ZOLTAN_DD_GID_NOT_FOUND_ERROR ;

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After fill return lists.");

fini:
   LB_FREE (&sbuff) ;
   LB_FREE (&rbuff) ;
   LB_FREE (&procs) ;
   LB_Comm_Destroy (&plan) ;

   if (err != ZOLTAN_DD_NORMAL_RETURN)
      ZOLTAN_PRINT_WARN (dd->my_proc, yo, "return is not normal") ;

   if (dd->debug_level > 0)
      {
      sprintf (str, "Processed %d GIDs, GIDs not found: %d", count, errcount) ;
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, str) ;
      }

   if (dd->debug_level > 1)
      ZOLTAN_TRACE_EXIT(dd->my_proc, yo, NULL);
   return err ;
   }







/******************  DD_Find_Local()  ***************************/

/* For a given gid, DD_Find_Local() provides its local ID, owner, optional
// user data, and partition from the local distributed directory. An error
// is returned if the gid is not found.
*/

static int DD_Find_Local (Zoltan_DD_Directory *dd,
 LB_ID_PTR gid,         /* incoming GID to locate (in)            */
 LB_ID_PTR lid,         /* gid's LID (out)                        */
 LB_ID_PTR user,        /* gid's user data (out)                  */
 int *partition,        /* gid's partition number (out)           */
 int *owner)            /* gid's owner (processor number) (out)   */
   {
   DD_Node *ptr ;
   int index ;
   char *yo = "DD_Find_Local" ;

   /* input sanity check */
   if (dd == NULL || owner == NULL || gid == NULL)
      {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument") ;
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   /* compute offset into hash table to find head of linked list */
   index = DD_Hash2 (gid, dd->gid_length, dd->table_length) ;

   /* walk link list until end looking for matching global ID */
   for (ptr = dd->table[index] ; ptr != NULL ; ptr = ptr->next)
      if (LB_EQ_ID (dd->gid_length, gid, ptr->gid) == TRUE)
         {
         /* matching global ID found! Return gid's information */
         LB_SET_ID (dd->lid_length, lid, ptr->gid + dd->gid_length) ;
         LB_SET_ID (dd->user_data_length, user,
          ptr->gid + (dd->gid_length + dd->lid_length)) ;

         if (owner     != NULL)   *owner     = ptr->owner ;
         if (partition != NULL)   *partition = ptr->partition ;

         return ZOLTAN_DD_NORMAL_RETURN ;
         }

   if (dd->debug_level > 0)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "Failed to find GID") ;

   return ZOLTAN_DD_GID_NOT_FOUND_ERROR ;
   }
