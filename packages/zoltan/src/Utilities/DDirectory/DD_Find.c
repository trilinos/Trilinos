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
 
   /* input sanity check */
   if (dd == NULL || count < 1 || owner == NULL || gid == NULL)
      return ZOLTAN_DD_INPUT_ERROR ;

   /* allocate memory for processors to contact for directory info */
   procs = (int *) LB_MALLOC (count * sizeof (int)) ;
   if (procs == NULL)
      return ZOLTAN_DD_MEMORY_ERROR ;

   /* allocate total memory for DD_Find_Msg send buffer */
   sbuff = (char *) LB_MALLOC (count * dd->find_msg_size) ;
   if (sbuff == NULL)
      {
      LB_FREE (&procs) ;
      return ZOLTAN_DD_MEMORY_ERROR ;
      }

   /* for each GID, fill DD_Find_Msg buffer and contact list */
   for (i = 0 ; i < count ; i++) 
      {
      procs[i] = dd->hash (gid + i*dd->gid_length, dd->gid_length,
       dd->nproc) ;

      ptr = (DD_Find_Msg *) (sbuff + i*dd->find_msg_size) ;
      ptr->index = i ;
      ptr->proc  = procs[i] ;
      LB_SET_ID (dd->gid_length, ptr->id, gid + i*dd->gid_length) ;
      }
   
   /* create efficient communication plan */
   err = LB_Comm_Create (&plan, count, procs, dd->comm,
    ZOLTAN_DD_FIND_MSG_TAG, &nrec) ;
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

   /* get find messages directed to me, fill in return information */
   errcount = 0 ;
   for (i = 0 ; i < nrec ; i++)
      {
      ptr = (DD_Find_Msg *) (rbuff + i*dd->find_msg_size) ;
      err = DD_Find_Local (dd, ptr->id, ptr->id, ptr->id
       + dd->max_id_length, &ptr->partition, &ptr->proc) ;
      if (err == ZOLTAN_DD_GID_NOT_FOUND_ERROR)
         errcount++ ;
      }

   /* send return information back to requestor */
   err = LB_Comm_Do_Reverse(plan, ZOLTAN_DD_FIND_MSG_TAG+2, rbuff,
    dd->find_msg_size, NULL, sbuff) ;
   if (err != COMM_OK)
      goto fini ;

   /* fill in user supplied lists with returned information */
   for (i = 0 ; i < count ; i++)
      {
      ptr = (DD_Find_Msg *) (sbuff + i*dd->find_msg_size) ;
      owner[ptr->index] = ptr->proc ;
      if (lid != NULL)
         LB_SET_ID (dd->lid_length, lid + ptr->index * dd->lid_length,
           ptr->id) ;

      if (data != NULL)
         LB_SET_ID (dd->user_data_length, data + ptr->index
          * dd->user_data_length, ptr->id + dd->max_id_length) ;
      }

   /* if at least one GID was not found, notify caller of error */
   err = (errcount == 0) ? ZOLTAN_DD_NORMAL_RETURN
                         : ZOLTAN_DD_GID_NOT_FOUND_ERROR ;

fini:
   LB_FREE (&sbuff) ;
   LB_FREE (&rbuff) ;
   LB_FREE (&procs) ;
   LB_Comm_Destroy (&plan) ;

   if (dd->debug_level > 0)
      printf ("ZOLTAN_DD_FIND(%d): Completed %d, Problem count %d\n",
       dd->my_proc, count, errcount) ;

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

   /* input sanity check */
   if (dd == NULL || owner == NULL || gid == NULL)
      return ZOLTAN_DD_INPUT_ERROR ;

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
      printf ("ZOLTAN_DD_FIND_LOCAL(%d): FAILED TO FIND GID\n",
       dd->my_proc) ;

   return ZOLTAN_DD_GID_NOT_FOUND_ERROR ;
   }
