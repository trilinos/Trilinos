/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"



/*  NOTE: See file, README, for associated documentation. (RTH) */



static int DD_Remove_Local (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid) ;




/********************   Zoltan_DD_Remove()  ***********************/

int Zoltan_DD_Remove (
 Zoltan_DD_Directory *dd,            /* directory state infomation      */
 ZOLTAN_ID_PTR gid,                  /* Incoming list of GIDs to remove */
 int count)                          /* Number of GIDs in removal list  */
   {
   int             *procs = NULL ;   /* list of processors to contact   */
   DD_Remove_Msg   *ptr   = NULL ;
   ZOLTAN_COMM_OBJ *plan  = NULL ;   /* efficient MPI communication     */
   char            *sbuff = NULL ;   /* send buffer                     */
   char            *rbuff = NULL ;   /* receive buffer                  */

   int              nrec ;           /* number of receives to expect    */
   int              i ;
   int              err ;            /* error condition to return       */
   int              errcount ;       /* count of GIDs not found         */
   char             str[100] ;       /* string to build error messages  */
   char            *yo = "Zoltan_DD_Remove" ;


   if (dd != NULL && dd->debug_level > 1)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL) ;

   /* input sanity checks */
   if (dd == NULL || count < 0 || (gid == NULL && count > 0))
      {
      ZOLTAN_PRINT_ERROR ((dd == NULL) ? ZOLTAN_DD_NO_PROC : dd->my_proc,
       yo, "Invalid input argument.") ;
      if (dd != NULL && dd->debug_level > 1)
         ZOLTAN_TRACE_OUT ((dd == NULL) ? ZOLTAN_DD_NO_PROC : dd->my_proc,
          yo, NULL) ;
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   /* allocate memory for processor contact list */
   if (count > 0)
      {
      procs = (int *) ZOLTAN_MALLOC (sizeof (int) * count);
      if (procs == NULL)
         {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc proc list.") ;
         if (dd->debug_level > 1)
            ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;
         return ZOLTAN_DD_MEMORY_ERROR ;
         }
      }

   /* allocate memory for DD_Remove_Msg send buffer */
   if (count > 0)
      {
      sbuff = (char *) ZOLTAN_MALLOC (dd->remove_msg_size * count) ;
      if (sbuff == NULL)
         {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc send buffer.") ;
         err = ZOLTAN_DD_MEMORY_ERROR ;
         goto fini ;
         }
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "After proc & sbuff mallocs.") ;

   /* for each GID, fill in contact list and then message structure */
   for (i = 0 ; i < count ; i++)
      {
      procs[i] = dd->hash(gid + i*dd->gid_length, dd->gid_length, dd->nproc) ;

      ptr = (DD_Remove_Msg *) (sbuff + i * dd->remove_msg_size) ;
      ptr->owner = dd->my_proc ;
      ZOLTAN_SET_ID (dd->gid_length, ptr->gid, gid + i * dd->gid_length) ;
      }

   /* now create efficient communication plan */
   err = Zoltan_Comm_Create (&plan, count, procs, dd->comm,
    ZOLTAN_DD_REMOVE_MSG_TAG, &nrec) ;
   if (err != ZOLTAN_OK)
      {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "COMM Create error.") ;
      goto fini ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "After Zoltan_Comm_Create.") ;

   /* allocate receive buffer for nrec DD_Remove_Msg structures */
   if (nrec > 0)
      {
      rbuff = (char *) ZOLTAN_MALLOC (nrec * dd->remove_msg_size) ;
      if (rbuff == NULL)
         {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Receive buffer malloc failed.") ;
         err = ZOLTAN_DD_MEMORY_ERROR ;
         goto fini ;
         }
      }

   /* send my remove messages & receive removes directed to me */
   err = Zoltan_Comm_Do (plan, ZOLTAN_DD_UPDATE_MSG_TAG+1, sbuff,
    dd->remove_msg_size, rbuff) ;
   if (err != ZOLTAN_OK)
      {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "COMM Do error.") ;
      goto fini ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "After Zoltan_Comm_Do.") ;

   /* for each message rec'd,  remove local directory info */
   errcount = 0 ;
   for (i = 0 ; i < nrec ; i++)
      {
      ptr = (DD_Remove_Msg *) (rbuff + i * dd->remove_msg_size) ;

      err = DD_Remove_Local (dd, ptr->gid) ;
      if (err == ZOLTAN_DD_GID_NOT_FOUND_ERROR)
         errcount++ ;
      }
   err = (errcount == 0) ? ZOLTAN_DD_NORMAL_RETURN
                         : ZOLTAN_DD_GID_NOT_FOUND_ERROR ;

   /* done, now free up things and return */
 fini:
   ZOLTAN_FREE (&procs) ;
   ZOLTAN_FREE (&sbuff) ;
   ZOLTAN_FREE (&rbuff) ;
   Zoltan_Comm_Destroy (&plan) ;

   if (dd->debug_level > 0)
      {
      sprintf (str, "Processed %d GIDs (%d local), %d GIDs not found.",
       count, nrec, errcount) ;
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, str) ;
      }
   if (dd->debug_level > 1)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;

   return err ;
   }







/******************  DD_Remove_Local()  **************************/

/* The given global ID, gid, is removed from the local distributed
// directory. An error is returned if the gid is not found.
*/

static int DD_Remove_Local (Zoltan_DD_Directory *dd,
 ZOLTAN_ID_PTR gid)                /* GID to be removed (in)  */
   {
   DD_Node **ptr ;
   DD_Node  *old ;
   int index ;
   char *yo = "DD_Remove_Local" ;

   /* input sanity checking */
   if (dd == NULL || gid == NULL)
      {
      ZOLTAN_PRINT_ERROR ((dd == NULL) ? ZOLTAN_DD_NO_PROC : dd->my_proc,
       yo, "Invalid input argument.") ;
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   if (dd->debug_level > 2)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL) ;

   /* compute offset into hash table to find head of linked list */
   index = Zoltan_DD_Hash2 (gid, dd->gid_length, dd->table_length) ;

   /* walk linked list until end looking for matching gid (key) */
   for (ptr = dd->table + index ; *ptr != NULL ; ptr = &((*ptr)->next))
      {
      if (ZOLTAN_EQ_ID(dd->gid_length, gid, (*ptr)->gid) == TRUE)
         {
         /* found node to remove, need to preserve its next ptr. */
          old =  *ptr ;
         *ptr = (*ptr)->next ;
         ZOLTAN_FREE (&old) ;       /* now OK to delete node */

         if (dd->debug_level > 2)
            ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;

         return ZOLTAN_DD_NORMAL_RETURN ;
         }
      }

   /* We get here only if the global ID has not been found */
   if (dd->debug_level > 2)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL) ;

   return ZOLTAN_DD_GID_NOT_FOUND_ERROR ;
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
