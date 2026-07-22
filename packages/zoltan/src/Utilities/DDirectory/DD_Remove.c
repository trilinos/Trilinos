// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include <stdio.h>
#include <stdlib.h>

#include "zoltan_dd_const.h"
#include "DD_Memory.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*  NOTE: See file, README, for associated documentation. (RTH) */



static int DD_Remove_Local (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid) ;



/********************   Zoltan_DD_Remove()  ***********************/

int Zoltan_DD_Remove (
 Zoltan_DD_Directory *dd,            /* directory state infomation      */
 ZOLTAN_ID_PTR gid,                  /* Incoming list of GIDs to remove */
 int count)                          /* Number of GIDs in removal list  */
{
   int             *procs = NULL;   /* list of processors to contact   */
   DD_Remove_Msg   *ptr   = NULL;
   ZOLTAN_COMM_OBJ *plan  = NULL;   /* efficient MPI communication     */
   char            *sbuff = NULL;   /* send buffer                     */
   char            *sbufftmp = NULL;/* pointer into send buffer        */
   char            *rbuff = NULL;   /* receive buffer                  */
   char            *rbufftmp = NULL;/* pointer into receive buffer     */

   int              nrec;           /* number of receives to expect    */
   int              i;
   int              err;            /* error condition to return       */
   int              errcount;       /* count of GIDs not found         */
   char             str[100];       /* string to build error messages  */
   char            *yo = "Zoltan_DD_Remove";



   /* input sanity checks */
   if (dd == NULL || count < 0 || (gid == NULL && count > 0)) {
      ZOLTAN_PRINT_ERROR (dd ? dd->my_proc : ZOLTAN_DD_NO_PROC, yo,
       "Invalid input argument");
      return ZOLTAN_FATAL;
   }
   if (dd->debug_level > 4)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL);


   /* allocate memory for processor contact list */
   if (count)  {
      procs = (int*) ZOLTAN_MALLOC (sizeof(int) * count);
      if (procs == NULL)  {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc proc list");
         if (dd->debug_level > 4)
            ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);
         return ZOLTAN_MEMERR;
      }
   }

   /* allocate memory for DD_Remove_Msg send buffer */
   if (count) {
      sbuff = (char*)ZOLTAN_MALLOC((size_t)(dd->remove_msg_size)*(size_t)count);
      if (sbuff == NULL)  {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc send buffer");
         err = ZOLTAN_MEMERR;
         goto fini;
      }
   }

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "After proc & sbuff mallocs");

   /* for each GID, fill in contact list and then message structure */
   sbufftmp = sbuff;
   for (i = 0; i < count; i++)  {
      procs[i] = dd->hash(gid + i*dd->gid_length, dd->gid_length, dd->nproc,
                          dd->hashdata, dd->hashfn);
      ptr = (DD_Remove_Msg*) sbufftmp;
      sbufftmp += dd->remove_msg_size;
      ptr->owner = dd->my_proc;
      ZOLTAN_SET_ID (dd->gid_length, ptr->gid, gid + i * dd->gid_length);
   }

   /* now create efficient communication plan */
   err = Zoltan_Comm_Create (&plan, count, procs, dd->comm,
    ZOLTAN_DD_REMOVE_MSG_TAG, &nrec);
   if (err != ZOLTAN_OK)  {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Comm_Create error");
      goto fini;
   }
   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "After Zoltan_Comm_Create");

   /* allocate receive buffer for nrec DD_Remove_Msg structures */
   if (nrec) {
      rbuff = (char*)ZOLTAN_MALLOC((size_t)nrec*(size_t)(dd->remove_msg_size));
      if (rbuff == NULL)  {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Receive buffer malloc failed");
         err = ZOLTAN_MEMERR;
         goto fini;
      }
   }

   /* send my remove messages & receive removes directed to me */
   err = Zoltan_Comm_Do (plan, ZOLTAN_DD_UPDATE_MSG_TAG+1, sbuff,
    dd->remove_msg_size, rbuff);
   if (err != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Comm_Do error");
      goto fini;
   }
   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "After Zoltan_Comm_Do");

   /* for each message rec'd,  remove local directory info */
   errcount = 0;
   rbufftmp = rbuff;
   for (i = 0; i < nrec; i++)  {
      ptr = (DD_Remove_Msg*) rbufftmp;
      rbufftmp += dd->remove_msg_size;

      err = DD_Remove_Local (dd, ptr->gid);
      if (err == ZOLTAN_WARN)
         ++errcount;
   }

   err = ZOLTAN_OK;
   if (dd->debug_level)  {
      sprintf (str, "Processed %d GIDs (%d local), %d GIDs not found",
       count, nrec, errcount);
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, str);
      err = (errcount) ? ZOLTAN_WARN : ZOLTAN_OK;
   }

   /* done, now free up things and return */
 fini:
   ZOLTAN_FREE (&procs);
   ZOLTAN_FREE (&sbuff);
   ZOLTAN_FREE (&rbuff);
   Zoltan_Comm_Destroy (&plan);

   if (dd->debug_level > 4)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);
   return err;
}







/******************  DD_Remove_Local()  **************************/

/* The given global ID, gid, is removed from the local distributed
 * directory. An error is returned if the gid is not found.
*/

static int DD_Remove_Local (Zoltan_DD_Directory *dd,
 ZOLTAN_ID_PTR gid)                /* GID to be removed (in)  */
{
   DD_Node *ptr;
   DD_NodeIdx nodeidx, prevnodeidx;
   int index;
   char *yo = "DD_Remove_Local";

   /* input sanity checking */
   if (dd == NULL || gid == NULL)  {
      ZOLTAN_PRINT_ERROR ((dd == NULL) ? ZOLTAN_DD_NO_PROC : dd->my_proc,
       yo, "Invalid input argument");
      return ZOLTAN_FATAL;
   }
   if (dd->debug_level > 5)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL);

   /* compute offset into hash table to find head of linked list */
   index = Zoltan_DD_Hash2 (gid, dd->gid_length, dd->table_length,
                            dd->hashdata, NULL);

   /* walk linked list until end looking for matching gid (key) */
   prevnodeidx = -1;
   for (nodeidx = dd->table[index]; nodeidx != -1;
      nodeidx = dd->nodelist[nodeidx].next) {
      ptr = dd->nodelist + nodeidx;
      if (ZOLTAN_EQ_ID(dd->gid_length, gid, ptr->gid) == TRUE)  {
         /* found node to remove, need to preserve its next ptr */
         if (prevnodeidx != -1)
            dd->nodelist[prevnodeidx].next = ptr->next;
         else
            dd->table[index] = ptr->next;
         DD_Memory_Free_Node(dd, nodeidx);       /* now OK to delete node */

         if (dd->debug_level > 5)
            ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);
         return ZOLTAN_OK;
      }
      prevnodeidx = nodeidx;
   }

   /* We get here only if the global ID has not been found */
   if (dd->debug_level > 5)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);
   return ZOLTAN_WARN;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
