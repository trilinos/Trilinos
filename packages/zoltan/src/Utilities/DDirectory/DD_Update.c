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






static int DD_Update_Local (Zoltan_DD_Directory *dd, LB_ID_PTR gid,
 LB_ID_PTR lid, LB_ID_PTR user, int partition, int owner) ;




/******************   Zoltan_DD_Update()  *************************/

int Zoltan_DD_Update (Zoltan_DD_Directory *dd,
 LB_ID_PTR gid,          /* Incoming list of GIDs to update          */
 LB_ID_PTR lid,          /* Incoming corresponding LIDs (optional)   */
 LB_ID_PTR user,         /* Incoming list of user data (optional)    */
 int *partition,         /* Optional, grouping of GIDs to partitions */
 int count)              /* Number of GIDs in update list            */
   {
   int             *procs = NULL ;   /* list of processors to contact */
   DD_Update_Msg   *ptr   = NULL ;
   struct Comm_Obj *plan  = NULL ;   /* for efficient MPI communication */
   char            *sbuff = NULL ;   /* send buffer                     */
   char            *rbuff = NULL ;   /* receive buffer                  */
   DD_Node         *ddptr = NULL ;
   int              nrec ;           /* number of receives to expect    */
   int              i ;
   int              err ;
   int              errcount ;

 
   /* input sanity checking */
   if (dd == NULL || count < 0 || gid == NULL)
      return ZOLTAN_DD_INPUT_ERROR ;

   /* for each linked list head, walk its list resetting errcheck */
   for (i = 0 ; i < dd->table_length ; i++)
      for (ddptr = dd->table[i] ; ddptr != NULL ; ddptr = ddptr->next)
         ddptr->errcheck = ZOLTAN_DD_NO_PROC ;

   /* allocate memory for list of processors to contact */
   procs = (int *) LB_MALLOC (count * sizeof (int)) ;
   if (procs == NULL)
      return ZOLTAN_DD_MEMORY_ERROR ;

   /* allocate memory for DD_Update_Msg send buffer */
   sbuff = (char *) LB_MALLOC (count * dd->update_msg_size) ;
   if (sbuff == NULL)
      {
      err = ZOLTAN_DD_MEMORY_ERROR ;
      goto fini ;
      }

   /* for each GID given, fill in contact list and then message structure */
   for (i = 0 ; i < count ; i++)
      {
      procs[i]= dd->hash(gid + i*dd->gid_length, dd->gid_length, dd->nproc);

      ptr = (DD_Update_Msg *) (sbuff + i * dd->update_msg_size) ;

      ptr->owner     = dd->my_proc ;
      ptr->partition = (partition == NULL) ? 0 : *(partition + i) ;

      LB_SET_ID (dd->gid_length, ptr->gid, gid + i * dd->gid_length) ;
      LB_SET_ID (dd->lid_length, ptr->gid + dd->gid_length,
       lid + i * dd->lid_length) ;
      LB_SET_ID (dd->user_data_length, ptr->gid +( dd->gid_length
       + dd->lid_length), user + i * dd->user_data_length) ;
      }

   /* now create efficient communication plan */
   err = LB_Comm_Create (&plan, count, procs, dd->comm,
    ZOLTAN_DD_UPDATE_MSG_TAG, &nrec) ;
   if (err != COMM_OK)
      goto fini ;

   /* allocate receive buffer for nrec DD_Update_Msg structures */
   rbuff = (char *) LB_MALLOC (nrec * dd->update_msg_size) ;
   if (rbuff == NULL)
      {
      err = ZOLTAN_DD_MEMORY_ERROR ;
      goto fini ;
      }

   /* send my update messages & receive updates directed to me */
   err = LB_Comm_Do (plan, ZOLTAN_DD_UPDATE_MSG_TAG+1, sbuff,
    dd->update_msg_size, rbuff) ;
   if (err != COMM_OK)
      goto fini ;
 
   /* for each message rec'd, update local directory information */
   errcount = 0 ;
   for (i = 0 ; i < nrec ; i++)
      {
      ptr = (DD_Update_Msg *) (rbuff + i * dd->update_msg_size) ;

      err = DD_Update_Local (dd, ptr->gid, ptr->gid + dd->gid_length,
       ptr->gid + (dd->gid_length + dd->lid_length), ptr->partition,
       ptr->owner) ;

      if (err != ZOLTAN_DD_NORMAL_RETURN && err != ZOLTAN_DD_GID_ADDED)
         errcount++ ;
      }

   /* successfully done, now free up things and return */
   err = (errcount == 0) ? ZOLTAN_DD_NORMAL_RETURN
                         : ZOLTAN_DD_GID_REDEFINED_ERROR ;

fini:
   LB_FREE (&procs) ;
   LB_FREE (&sbuff) ;
   LB_FREE (&rbuff) ;
   LB_Comm_Destroy (&plan) ;

   if (dd->debug_level > 0)
      printf ("ZOLTAN_DD_UPDATE(%d): processed %d updates, %d local, err cnt %d\n",
       dd->my_proc, count, nrec, errcount) ;

   return err ;
   }

/*****************  DD_Update_Local ()  **********************/

/* DD_Update_Local() unconditionally updates the information associated
// with a global ID. If the gid was not found, it is added to the directory.
*/

static int DD_Update_Local (Zoltan_DD_Directory *dd,
 LB_ID_PTR gid,          /* GID to update (in)                        */
 LB_ID_PTR lid,          /* gid's LID (in), NULL if not needed        */
 LB_ID_PTR user,         /* gid's user data (in), NULL if not needed  */
 int partition,          /* gid's partition (in), 0 if not used       */
 int owner)              /* gid's current owner (proc number) (in)    */
   {
   DD_Node *ptr, *old ;
   int index ;
   int head = TRUE ;       /* used to detect head of link list        */
   char err[] = "ZOLTAN_UPDATE_LOCAL(%d): multiply defined global ID" ;

   /* input sanity checking */
   if (dd == NULL || owner  < 0 || owner >= dd->nproc || gid == NULL)
      return ZOLTAN_DD_INPUT_ERROR ;

   /* compute offset into hash table to find head of linked list */
   index = DD_Hash2 (gid, dd->gid_length, dd->table_length) ;

   /* walk linked list until end looking for matching gid */
   for (ptr = dd->table[index] ; ptr != NULL ; ptr = ptr->next)
       {
       if (LB_EQ_ID (dd->gid_length, gid, ptr->gid) == TRUE)
          {
          /* found match, unconditionally update directory information */
          LB_SET_ID (dd->lid_length,        ptr->gid +  dd->gid_length,
           lid) ;
          LB_SET_ID  (dd->user_data_length, ptr->gid + (dd->gid_length
           + dd->lid_length), user) ;

          ptr->owner     = owner ;
          ptr->partition = partition ;

          /* Response to multiple updates to a gid in 1 update cycle */
          if (dd->debug_level == 0 || ptr->errcheck == ZOLTAN_DD_NO_PROC)
             {
             ptr->errcheck = owner ;
             return ZOLTAN_DD_NORMAL_RETURN ;  /* ignore all errors */
             }

          if (dd->debug_level == 1 && ptr->errcheck != owner)
             return ZOLTAN_DD_GID_REDEFINED_ERROR ; /* err return only */

          if (dd->debug_level == 2 && ptr->errcheck != owner)
             {
             fprintf (stderr, err, dd->my_proc) ;
             LB_PRINT_ID (dd->gid_length, ptr->gid) ;
             LB_PRINT_ID (dd->gid_length, gid) ;
             return ZOLTAN_DD_GID_REDEFINED_ERROR ;
             }

          fprintf (stderr, err, dd->my_proc) ;
          LB_PRINT_ID (dd->gid_length, ptr->gid) ;
          LB_PRINT_ID (dd->gid_length, gid) ;
          return ZOLTAN_DD_GID_REDEFINED_ERROR ;
          }
       head = FALSE ;     /* no longer at head of list                */
       old  = ptr ;       /* save previous in case we go one to far   */
       }

   /* gid not found. Create new DD_Node and fill it in */
   ptr = (DD_Node *) LB_MALLOC (dd->node_size)  ;
   if (ptr == NULL)
      return ZOLTAN_DD_MEMORY_ERROR ;

   if (head == TRUE)
      dd->table [index] = ptr ;  /* set pointer at head of linked list */
   else
      old->next         = ptr ;  /* set pointer in previous ptr->next  */
  
   LB_SET_ID (dd->gid_length,       ptr->gid,                   gid) ;
   LB_SET_ID (dd->lid_length,       ptr->gid + dd->gid_length,  lid) ;
   LB_SET_ID (dd->user_data_length, ptr->gid + (dd->gid_length
                                    + dd->lid_length),          user) ;

   ptr->next      = NULL ;
   ptr->owner     = owner ;
   ptr->partition = partition ;
   ptr->errcheck  = owner ;

   return ZOLTAN_DD_GID_ADDED ;
   }
