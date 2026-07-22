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


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


int Zoltan_DD_GetLocalKeys(Zoltan_DD_Directory *dd,
                           ZOLTAN_ID_PTR *gid,
                           int *size)
{
  int ierr = ZOLTAN_OK;
  int i, k;
  DD_NodeIdx nodeidx;
  DD_Node *ptr;
  int gid_alloc_size;

  gid_alloc_size = dd->table_length;
  (*gid) = (ZOLTAN_ID_PTR)ZOLTAN_MALLOC(
                          gid_alloc_size*dd->gid_length*sizeof(ZOLTAN_ID_TYPE));

  k= 0;
  for (i = 0; i < dd->table_length; i++)
    for (nodeidx = dd->table[i]; nodeidx != -1;
         nodeidx = dd->nodelist[nodeidx].next) {
      ptr = dd->nodelist + nodeidx;
      if (k >= gid_alloc_size) {
        gid_alloc_size *= 2;
        (*gid) = (ZOLTAN_ID_PTR)
                  ZOLTAN_REALLOC((*gid),
                         gid_alloc_size*dd->gid_length*sizeof(ZOLTAN_ID_TYPE));
      }
      ZOLTAN_SET_ID (dd->gid_length, (*gid)+k*dd->gid_length, ptr->gid);
      k++;
    }

  (*size) = k;

  return (ierr);
}

/*  NOTE: See file, README, for associated documentation. (RTH) */


static int DD_Find_Local (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
 ZOLTAN_ID_PTR lid, char *user, int *partition, int *owner) ;



/********************  Zoltan_DD_Find()  **************************/

int Zoltan_DD_Find (
 Zoltan_DD_Directory *dd, /* contains directory state information        */
 ZOLTAN_ID_PTR gid,       /* Incoming list of GIDs to get owners proc    */
 ZOLTAN_ID_PTR lid,       /* Outgoing corresponding list of LIDs         */
 char *data,              /* Outgoing optional corresponding user data   */
 int *partition,          /* Outgoing optional partition information     */
 int  count,              /* Count of GIDs in above list (in)            */
 int *owner)              /* Outgoing optional list of data owners       */
{
   ZOLTAN_COMM_OBJ *plan  = NULL;     /* efficient MPI communication     */
   char            *rbuff = NULL;     /* receive buffer                  */
   char            *rbufftmp = NULL;  /* pointer into receive buffer     */
   char            *sbuff = NULL;     /* send buffer                     */
   char            *sbufftmp = NULL;  /* pointer into send buffer        */
   int             *procs = NULL;     /* list of processors to contact   */
   DD_Find_Msg     *ptr   = NULL;
   int              i;
   int              nrec;             /* number of messages to receive   */
   int              err = ZOLTAN_OK;  /* return error condition          */
   int              errcount;         /* count of GIDs not found         */
   char            *yo = "Zoltan_DD_Find";


   /* input sanity check */
   if (dd == NULL || count < 0 || (gid == NULL && count > 0))  {
      ZOLTAN_PRINT_ERROR (dd ? dd->my_proc : ZOLTAN_DD_NO_PROC, yo,
       "Invalid input argument");
      return ZOLTAN_FATAL;
   }
   if (dd->debug_level > 4)
      ZOLTAN_TRACE_IN(dd->my_proc, yo, NULL);

   /* allocate memory for processors to contact for directory info */
   if (count)  {
      procs = (int*) ZOLTAN_MALLOC (sizeof(int) * count);
      if (procs == NULL) {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc proc list");
         if (dd->debug_level > 4)
           ZOLTAN_TRACE_OUT(dd->my_proc, yo, NULL);
         return ZOLTAN_MEMERR;
      }
   }

   /* allocate memory for DD_Find_Msg send buffer */
   if (count)  {
      sbuff = (char*) ZOLTAN_CALLOC (count, dd->find_msg_size);
      if (sbuff == NULL)  {
         ZOLTAN_FREE (&procs);
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc send buffer");
         if (dd->debug_level > 4)
            ZOLTAN_TRACE_OUT(dd->my_proc, yo, NULL);
         return ZOLTAN_MEMERR;
      }
   }

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After mallocs");

   /* for each GID, fill DD_Find_Msg buffer and contact list */
   sbufftmp = sbuff;
   for (i = 0; i < count; i++)  {
      procs[i] = dd->hash (gid + i*dd->gid_length, dd->gid_length, dd->nproc,
                           dd->hashdata, dd->hashfn);
      ptr      = (DD_Find_Msg*) sbufftmp;
      sbufftmp += dd->find_msg_size;

      ptr->index = i;
      ptr->proc  = procs[i];
      ZOLTAN_SET_ID (dd->gid_length, ptr->id, gid + i*dd->gid_length);
   }
   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After fill");

   /* create efficient communication plan */
   err = Zoltan_Comm_Create (&plan, count, procs, dd->comm,
    ZOLTAN_DD_FIND_MSG_TAG, &nrec);
   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Create");
   if (err != ZOLTAN_OK)
      goto fini;

   /* allocate receive buffer */
   if (nrec)  {
      rbuff = (char*) ZOLTAN_MALLOC ((size_t)nrec*(size_t)(dd->find_msg_size));
      if (rbuff == NULL)  {
         err = ZOLTAN_MEMERR;
         goto fini;
      }
   }

   /* send out find messages across entire system */
   err = Zoltan_Comm_Do (plan, ZOLTAN_DD_FIND_MSG_TAG+1, sbuff,
    dd->find_msg_size, rbuff);
   if (err != ZOLTAN_OK)
      goto fini;

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Do");

   /* get find messages directed to me, fill in return information */
   errcount = 0;
   rbufftmp = rbuff;
   for (i = 0; i < nrec; i++)  {
      ptr = (DD_Find_Msg*) rbufftmp;
      rbufftmp += dd->find_msg_size;
      err = DD_Find_Local (dd, ptr->id, ptr->id,
                           (char *)(ptr->id + dd->max_id_length),
                           &ptr->partition, &ptr->proc);
      if (err == ZOLTAN_WARN)
          ++errcount;
   }
   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After fill in return info");

   /* send return information back to requester */
   err = Zoltan_Comm_Do_Reverse(plan, ZOLTAN_DD_FIND_MSG_TAG+2, rbuff,
    dd->find_msg_size, NULL, sbuff);
   if (err != ZOLTAN_OK)
      goto fini;
   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Reverse");

   /* fill in user supplied lists with returned information */
   sbufftmp = sbuff;
   for (i = 0; i < count; i++) {
      ptr = (DD_Find_Msg*) sbufftmp;
      sbufftmp += dd->find_msg_size;

      if (owner)
         owner[ptr->index] = ptr->proc;
      if (partition)
         partition[ptr->index] = ptr->partition ;
      if (lid)
         ZOLTAN_SET_ID(dd->lid_length,lid+ptr->index*dd->lid_length,ptr->id);
      if (data)
         memcpy(data + (size_t)(ptr->index) * (size_t)(dd->user_data_length),
                ptr->id + dd->max_id_length, dd->user_data_length);
   }
   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After fill return lists");
/*    err = ZOLTAN_OK;     */

   MPI_Allreduce(&errcount, &err, 1, MPI_INT, MPI_SUM, dd->comm);
   err = (err) ? ZOLTAN_WARN : ZOLTAN_OK;

   /* if at least one GID was not found, potentially notify caller of error */
   if (dd->debug_level > 0)  {
      char str[100];      /* diagnostic message string */
      sprintf (str, "Processed %d GIDs, GIDs not found: %d", count, errcount);
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, str);
   }

fini:
   ZOLTAN_FREE (&sbuff);
   ZOLTAN_FREE (&rbuff);
   ZOLTAN_FREE (&procs) ;
   Zoltan_Comm_Destroy (&plan);

   if (dd->debug_level > 4)
      ZOLTAN_TRACE_OUT(dd->my_proc, yo, NULL);
   return err;
}







/******************  DD_Find_Local()  ***************************/

/* For a given gid, DD_Find_Local() provides its local ID, owner, optional
 * user data, and partition from the local distributed directory. An error
 * is returned if the gid is not found.
*/

static int DD_Find_Local (Zoltan_DD_Directory *dd,
 ZOLTAN_ID_PTR gid,         /* incoming GID to locate (in)            */
 ZOLTAN_ID_PTR lid,         /* gid's LID (out)                        */
 char *user,                /* gid's user data (out)                  */
 int *partition,            /* gid's partition number (out)           */
 int *owner)                /* gid's owner (processor number) (out)   */
{
   DD_Node *ptr;
   DD_NodeIdx nodeidx;
   int      index;
   char    *yo = "DD_Find_Local";

   /* input sanity check */
   if (dd == NULL || owner == NULL || gid == NULL)  {
      ZOLTAN_PRINT_ERROR ((dd == NULL) ? 0 : dd->my_proc, yo, "Invalid input");
      return ZOLTAN_FATAL;
   }
   if (dd->debug_level > 5)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL);

   /* compute offset into hash table to find head of linked list */
   index = Zoltan_DD_Hash2 (gid, dd->gid_length, dd->table_length,
                            dd->hashdata, NULL);
   /* walk link list until end looking for matching global ID */
   for (nodeidx = dd->table[index]; nodeidx != -1;
        nodeidx = dd->nodelist[nodeidx].next) {
      ptr = dd->nodelist + nodeidx;
      if (ZOLTAN_EQ_ID (dd->gid_length, gid, ptr->gid) == TRUE)  {
         /* matching global ID found! Return gid's information */
         if (lid) ZOLTAN_SET_ID(dd->lid_length, lid, ptr->gid + dd->gid_length);
         if (user) memcpy(user, ptr->gid + (dd->gid_length + dd->lid_length),
                          dd->user_data_length);

         if (owner)     *owner     = ptr->owner;
         if (partition) *partition = ptr->partition;

         if (dd->debug_level > 5)
            ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);
         return ZOLTAN_OK;
      }
   }


   if (owner != NULL)
      *owner = -1;    /* JDT Added -1 owner not found */
   if (dd->debug_level > 5)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);

   if (dd->debug_level > 0)  {
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "GID not found");
      return ZOLTAN_WARN;
   }
   return ZOLTAN_WARN;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
