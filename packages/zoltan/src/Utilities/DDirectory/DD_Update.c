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






static int DD_Update_Local (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
 ZOLTAN_ID_PTR lid, ZOLTAN_ID_PTR user, int partition, int owner) ;




/******************   Zoltan_DD_Update()  *************************/

int Zoltan_DD_Update (
 Zoltan_DD_Directory *dd,  /* directory state information              */
 ZOLTAN_ID_PTR gid,        /* Incoming list of GIDs to update          */
 ZOLTAN_ID_PTR lid,        /* Incoming corresponding LIDs (optional)   */
 ZOLTAN_ID_PTR user,       /* Incoming list of user data (optional)    */
 int *partition,           /* Optional, grouping of GIDs to partitions */
 int count)                /* Number of GIDs in update list            */
   {
   int             *procs = NULL ;   /* list of processors to contact   */
   DD_Update_Msg   *ptr   = NULL ;
   ZOLTAN_COMM_OBJ *plan  = NULL ;   /* for efficient MPI communication */
   char            *sbuff = NULL ;   /* send buffer                     */
   char            *rbuff = NULL ;   /* receive buffer                  */
   DD_Node         *ddptr = NULL ;
   int              nrec ;           /* number of receives to expect    */
   int              i ;
   int              err ;
   int              errcount ;       /* count of GIDs not found, added  */
   char             str[100] ;       /* build error message string      */
   char            *yo = "Zoltan_DD_Update" ;

   if (dd != NULL && dd->debug_level > 1)
      ZOLTAN_TRACE_IN(dd->my_proc, yo, NULL);

   /* input sanity checking */
   if (dd == NULL || count < 0 || (gid == NULL && count > 0))
      {
      ZOLTAN_PRINT_ERROR ((dd == NULL ? ZOLTAN_DD_NO_PROC : dd->my_proc), yo,
       "Invalid input argument.") ;
      if (dd != NULL && dd->debug_level > 1)
        ZOLTAN_TRACE_OUT((dd == NULL ? ZOLTAN_DD_NO_PROC : dd->my_proc), yo,
         NULL);
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   /* part of initializing the error checking process             */
   /* for each linked list head, walk its list resetting errcheck */
   for (i = 0 ; i < dd->table_length ; i++)
      for (ddptr = dd->table[i] ; ddptr != NULL ; ddptr = ddptr->next)
         ddptr->errcheck = ZOLTAN_DD_NO_PROC ;

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After reset errcheck.");

   /* allocate memory for list of processors to contact */
   if (count > 0)
      {
      procs = (int *) ZOLTAN_MALLOC (sizeof (int) * count) ;
      if (procs == NULL)
         {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc proc list.") ;
         if (dd->debug_level > 1)
            ZOLTAN_TRACE_OUT(dd->my_proc, yo, NULL);
         return ZOLTAN_DD_MEMORY_ERROR ;
         }
      }

   /* allocate memory for DD_Update_Msg send buffer */
   if (count > 0)
      {
      sbuff = (char *) ZOLTAN_MALLOC (dd->update_msg_size * count) ;
      if (sbuff == NULL)
         {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc send buffer.") ;
         err = ZOLTAN_DD_MEMORY_ERROR ;
         goto fini ;
         }
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After mallocs.");

   /* for each GID given, fill in contact list and then message structure */
   for (i = 0 ; i < count ; i++)
      {
      procs[i] = dd->hash(gid + i*dd->gid_length, dd->gid_length, dd->nproc);
      ptr      = (DD_Update_Msg *) (sbuff + i * dd->update_msg_size) ;

      ptr->owner     = dd->my_proc ;
      ptr->partition = (partition == NULL) ? 0 : *(partition + i) ;

      ZOLTAN_SET_ID (dd->gid_length, ptr->gid, gid + i * dd->gid_length) ;

      if (lid)
         ZOLTAN_SET_ID (dd->lid_length, ptr->gid + dd->gid_length, lid
          + i * dd->lid_length) ;
      if (user)
         ZOLTAN_SET_ID (dd->user_data_length, ptr->gid + (dd->gid_length
          + dd->lid_length), user + i * dd->user_data_length) ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After fill contact list.");

   /* now create efficient communication plan */
   err = Zoltan_Comm_Create (&plan, count, procs, dd->comm,
    ZOLTAN_DD_UPDATE_MSG_TAG, &nrec) ;
   if (err != ZOLTAN_OK)
      {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "COMM Create error") ;
      goto fini ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Create.");

   /* allocate receive buffer for nrec DD_Update_Msg structures */
   if (nrec > 0)
      {
      rbuff = (char *) ZOLTAN_MALLOC (nrec * dd->update_msg_size) ;
      if (rbuff == NULL)
         {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Receive buffer malloc failed.") ;
         err = ZOLTAN_DD_MEMORY_ERROR ;
         goto fini ;
         }
      }

   /* send my update messages & receive updates directed to me */
   err = Zoltan_Comm_Do (plan, ZOLTAN_DD_UPDATE_MSG_TAG+1, sbuff,
    dd->update_msg_size, rbuff) ;
   if (err != ZOLTAN_OK)
      {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "COMM Do error.") ;
      goto fini ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Do.");

   /* for each message rec'd, update local directory information */
   errcount = 0 ;
   for (i = 0 ; i < nrec ; i++)
      {
      ptr = (DD_Update_Msg *) (rbuff + i * dd->update_msg_size) ;

      err = DD_Update_Local (dd, ptr->gid,
       (lid)   ? (ptr->gid + dd->gid_length)                    : NULL,
       (user)  ? (ptr->gid + (dd->gid_length + dd->lid_length)) : NULL,
       (partition) ? (ptr->partition) : 0,
       ptr->owner) ;

      if (err != ZOLTAN_DD_NORMAL_RETURN && err != ZOLTAN_DD_GID_ADDED)
         errcount++ ;
      }

   /* successfully done, now free up things and return */
   err = (errcount == 0) ? ZOLTAN_DD_NORMAL_RETURN
                         : ZOLTAN_DD_GID_REDEFINED_ERROR ;

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After update.");

fini:
   ZOLTAN_FREE (&procs) ;
   ZOLTAN_FREE (&sbuff) ;
   ZOLTAN_FREE (&rbuff) ;
   Zoltan_Comm_Destroy (&plan) ;

   if (dd->debug_level > 0)
      {
      sprintf (str, "Processed %d GIDs (%d local), %d GID errors.",
       count, nrec, errcount) ;
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, str) ;
      }

   if (dd->debug_level > 1)
      ZOLTAN_TRACE_OUT(dd->my_proc, yo, NULL);

   return err ;
   }







/*****************  DD_Update_Local ()  **********************/

/* DD_Update_Local() unconditionally updates the information associated
// with a global ID. If the gid was not found, it is added to the directory.
*/

static int DD_Update_Local (Zoltan_DD_Directory *dd,
 ZOLTAN_ID_PTR gid,          /* GID to update (in)                        */
 ZOLTAN_ID_PTR lid,          /* gid's LID (in), NULL if not needed        */
 ZOLTAN_ID_PTR user,         /* gid's user data (in), NULL if not needed  */
 int partition,          /* gid's partition (in), 0 if not used       */
 int owner)              /* gid's current owner (proc number) (in)    */
   {
   DD_Node **ptr ;
   int index ;
   char *yo = "DD_Update_Local" ;


   if (dd != NULL && dd->debug_level > 2)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL) ;

   /* input sanity checking */
   if (dd == NULL || owner  < 0 || owner >= dd->nproc || gid == NULL)
      {
      ZOLTAN_PRINT_ERROR ((dd == NULL) ? ZOLTAN_DD_NO_PROC : dd->my_proc,
       yo, "Invalid input parameter.") ;

      if (dd != NULL && dd->debug_level > 2)
         ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;

      return ZOLTAN_DD_INPUT_ERROR ;
      }

   /* compute offset into hash table to find head of linked list */
   index = Zoltan_DD_Hash2 (gid, dd->gid_length, dd->table_length) ;

   /* walk linked list until end looking for matching gid */
   for (ptr = dd->table+index ; *ptr != NULL ; ptr = &((*ptr)->next))
       {
       if (ZOLTAN_EQ_ID (dd->gid_length, gid, (*ptr)->gid) == TRUE)
          {
          /* found match, update directory information */
          if (lid)
             ZOLTAN_SET_ID (dd->lid_length,(*ptr)->gid + dd->gid_length, lid) ;
          if (user)
             ZOLTAN_SET_ID (dd->user_data_length, (*ptr)->gid + (dd->gid_length
              + dd->lid_length), user) ;

          (*ptr)->owner     = owner ;
          (*ptr)->partition = partition ;

          /* Response to multiple updates to a gid in 1 update cycle */
          if (dd->debug_level == 0 || (*ptr)->errcheck == ZOLTAN_DD_NO_PROC)
             {
             (*ptr)->errcheck = owner ;
             if (dd->debug_level > 2)
                ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;
             return ZOLTAN_DD_NORMAL_RETURN ;  /* ignore all errors */
             }

          if (dd->debug_level == 1 && (*ptr)->errcheck != owner)
             return ZOLTAN_DD_GID_REDEFINED_ERROR ; /* err return only */

          if (dd->debug_level == 2 && (*ptr)->errcheck != owner)
             {
             ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Multiply defined GID.") ;
             return ZOLTAN_DD_GID_REDEFINED_ERROR ;
             }

          ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Multiply defined GID.") ;
          if (dd->debug_level > 2)
             ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;

          return ZOLTAN_DD_GID_REDEFINED_ERROR ;
          }
       }

   /* gid not found. Create new DD_Node and fill it in */
   *ptr = (DD_Node *) ZOLTAN_MALLOC (dd->node_size)  ;
   if (*ptr == NULL)
      {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc new Node.") ;
      if (dd->debug_level > 2)
         ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;
      return ZOLTAN_DD_MEMORY_ERROR ;
      }

   ZOLTAN_SET_ID (dd->gid_length, (*ptr)->gid, gid) ;

   if (lid)
      ZOLTAN_SET_ID (dd->lid_length,       (*ptr)->gid + dd->gid_length, lid) ;
   if (user)
      ZOLTAN_SET_ID (dd->user_data_length, (*ptr)->gid + (dd->gid_length
       + dd->lid_length), user) ;

   (*ptr)->next      = NULL ;
   (*ptr)->owner     = owner ;
   (*ptr)->partition = partition ;
   (*ptr)->errcheck  = owner ;

   if (dd->debug_level > 2)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;

   return ZOLTAN_DD_GID_ADDED ;
   }
