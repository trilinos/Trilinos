/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"
#include "DD_Memory.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


static int DD_Update_Local (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
 ZOLTAN_ID_PTR lid, char *user, int partition, int owner);


/******************   Zoltan_DD_Update()  *************************/

int Zoltan_DD_Update (
 Zoltan_DD_Directory *dd,  /* directory state information              */
 ZOLTAN_ID_PTR gid,        /* Incoming list of GIDs to update          */
 ZOLTAN_ID_PTR lid,        /* Incoming corresponding LIDs (optional)   */
 char *user,               /* Incoming list of user data (optional)    */
 int *partition,           /* Optional, grouping of GIDs to partitions */
 int count)                /* Number of GIDs in update list            */
{
   int             *procs = NULL;   /* list of processors to contact   */
   DD_Update_Msg   *ptr   = NULL;
   ZOLTAN_COMM_OBJ *plan  = NULL;   /* for efficient MPI communication */
   char            *sbuff = NULL;   /* send buffer                     */
   char            *sbufftmp = NULL;/* pointer into send buffer        */
   char            *rbuff = NULL;   /* receive buffer                  */
   char            *rbufftmp = NULL;/* pointer into receive buffer     */
   int              nrec;           /* number of receives to expect    */
   int              i;
   int              err;
   int              errcount = 0;   /* count of GIDs not found, added  */
   char             str[100];       /* build error message string      */
   char            *yo = "Zoltan_DD_Update";


   /* input sanity checking */
   if (dd == NULL || count < 0 || (gid == NULL && count > 0))  {
      ZOLTAN_PRINT_ERROR ((dd == NULL ? ZOLTAN_DD_NO_PROC : dd->my_proc), yo,
       "Invalid input argument");
      return ZOLTAN_FATAL;
   }
   if (dd->debug_level > 4)
      ZOLTAN_TRACE_IN(dd->my_proc, yo, NULL);

   /* part of initializing the error checking process             */
   /* for each linked list head, walk its list resetting errcheck */
   if (dd->debug_level)
      for (i = 0; i < dd->table_length; i++) {
         DD_NodeIdx nodeidx;
         for (nodeidx = dd->table[i]; nodeidx != -1;
              nodeidx = dd->nodelist[nodeidx].next)
            dd->nodelist[nodeidx].errcheck = ZOLTAN_DD_NO_PROC;
      }

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After reset errcheck");

   /* allocate memory for list of processors to contact */
   if (count) {
      procs = (int*) ZOLTAN_MALLOC (sizeof(int) * count);
      if (procs == NULL)  {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc proc list");
         err = ZOLTAN_MEMERR;
         goto fini;
      }
   }

   /* allocate memory for DD_Update_Msg send buffer */
   if (count)  {
      sbuff = (char*) ZOLTAN_CALLOC (count, dd->update_msg_size);
      if (sbuff == NULL)  {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to malloc send buffer");
         err = ZOLTAN_MEMERR;
         goto fini;
      }
   }

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After mallocs");

   /* for each GID given, fill in contact list and then message structure */
   sbufftmp = sbuff;
   for (i = 0; i < count; i++)  {
      procs[i] = dd->hash(gid + i*dd->gid_length, dd->gid_length, dd->nproc,
                          dd->hashdata, dd->hashfn);
      ptr      = (DD_Update_Msg*) sbufftmp;
      sbufftmp += dd->update_msg_size;

      ptr->lid_flag       = (lid)  ? 1 : 0;
      ptr->user_flag      = (user) ? 1 : 0;
      ptr->partition_flag = (partition) ? 1 : 0;
      ptr->partition      = (partition) ? *(partition + i) :  -1;
      ptr->owner          = dd->my_proc;

      ZOLTAN_SET_ID (dd->gid_length, ptr->gid, gid + i * dd->gid_length);
      if (lid) {
         ZOLTAN_SET_ID(dd->lid_length, ptr->gid + dd->gid_length,
                       lid + i * dd->lid_length);
      }
      else {
         memset(ptr->gid + dd->gid_length, 0, dd->lid_length);
      }
      if (user) {
         memcpy(ptr->gid + (dd->gid_length + dd->lid_length),
                user + (size_t)i * (size_t)(dd->user_data_length),
                dd->user_data_length);
      }
      else {
         memset(ptr->gid + (dd->gid_length + dd->lid_length), 0,
                dd->user_data_length);
      }
   }

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After fill contact list");

   /* now create efficient communication plan */
   err = Zoltan_Comm_Create (&plan, count, procs, dd->comm,
    ZOLTAN_DD_UPDATE_MSG_TAG, &nrec);
   if (err != ZOLTAN_OK)  {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Comm_Create error");
      goto fini;
   }

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Create");

   /* If dd has no nodes allocated (e.g., first call to DD_Update; 
    * create the nodelist and freelist 
    */
   if (nrec && dd->nodelistlen == 0) {
      DD_Memory_Alloc_Nodelist(dd, (DD_NodeIdx) nrec, 0.); 
                               /* TODO Add overalloc parameter */
   }

   /* allocate receive buffer for nrec DD_Update_Msg structures */
   if (nrec)  {
      rbuff = (char*)ZOLTAN_MALLOC((size_t)nrec*(size_t)(dd->update_msg_size));
      if (rbuff == NULL)  {
         ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Receive buffer malloc failed");
         err = ZOLTAN_MEMERR;
         goto fini;
      }
   }

   /* send my update messages & receive updates directed to me */
   err = Zoltan_Comm_Do (plan, ZOLTAN_DD_UPDATE_MSG_TAG+1, sbuff,
    dd->update_msg_size, rbuff);
   if (err != ZOLTAN_OK)  {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Comm_Do error");
      goto fini;
   }

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Comm_Do");

   /* for each message rec'd, update local directory information */
   errcount = 0;
   rbufftmp = rbuff;
   for (i = 0; i < nrec; i++)  {
      ptr = (DD_Update_Msg *) rbufftmp;
      rbufftmp += dd->update_msg_size;

      err = DD_Update_Local (dd, ptr->gid,
       (ptr->lid_flag)  ? (ptr->gid + dd->gid_length) : NULL,
       (char*)((ptr->user_flag)?(ptr->gid+(dd->gid_length+dd->lid_length)):NULL),
       (ptr->partition_flag) ? (ptr->partition) : -1,  /* illegal partition */
       ptr->owner);

      if (err != ZOLTAN_OK)
         ++errcount;
   }
   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO(dd->my_proc, yo, "After Local update");

   err = ZOLTAN_OK;
   if (dd->debug_level)  /* overwrite error return if extra checking is on */
      err = (errcount) ? ZOLTAN_WARN : ZOLTAN_OK;

fini:
   ZOLTAN_FREE (&procs);
   ZOLTAN_FREE (&sbuff);
   ZOLTAN_FREE (&rbuff);
   Zoltan_Comm_Destroy (&plan);

   if (dd->debug_level)  {
      sprintf (str, "Processed %d GIDs (%d local), %d GID errors", count,
       nrec, errcount);
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, str);
   }

   if (dd->debug_level > 4)
      ZOLTAN_TRACE_OUT(dd->my_proc, yo, NULL);
   return err;
}







/*****************  DD_Update_Local ()  **********************/

/* DD_Update_Local() unconditionally updates the information associated
 * with a global ID. If the gid was not found, it is added to the directory.
*/

static int DD_Update_Local (Zoltan_DD_Directory *dd,
 ZOLTAN_ID_PTR gid,         /* GID to update (in)                        */
 ZOLTAN_ID_PTR lid,         /* gid's LID (in), NULL if not needed        */
 char *user,                /* gid's user data (in), NULL if not needed  */
 int partition,             /* gid's partition (in), -1 if not used      */
 int owner)                 /* gid's current owner (proc number) (in)    */
{
   int index;
   char *yo = "DD_Update_Local";
   DD_NodeIdx nodeidx;
   DD_Node *ptr;

   /* input sanity checking */
   if (dd == NULL || owner  < 0 || owner >= dd->nproc || gid == NULL)  {
      ZOLTAN_PRINT_ERROR (dd ? dd->my_proc : ZOLTAN_DD_NO_PROC, yo,
       "Invalid input parameter");
      return ZOLTAN_FATAL;
   }
   if (dd->debug_level > 5)
      ZOLTAN_TRACE_IN (dd->my_proc, yo, NULL);

   /* compute offset into hash table to find head of linked list */
   index = Zoltan_DD_Hash2 (gid, dd->gid_length, dd->table_length,
                            dd->hashdata, NULL);

   /* walk linked list until end looking for matching gid */
   for (nodeidx = dd->table[index]; nodeidx != -1;
        nodeidx = dd->nodelist[nodeidx].next) {
       ptr = dd->nodelist + nodeidx;
       if (ZOLTAN_EQ_ID (dd->gid_length, gid, ptr->gid) == TRUE)  {
          /* found match, update directory information */
          if (lid)
             ZOLTAN_SET_ID (dd->lid_length,ptr->gid + dd->gid_length, lid);
          if (user)
             memcpy(ptr->gid + (dd->gid_length + dd->lid_length), user,
                    dd->user_data_length);

          ptr->owner = owner;
          if (partition != -1)
             ptr->partition = partition;

          /* Response to multiple updates to a gid in 1 update cycle */
          if (dd->debug_level > 0 && ptr->errcheck != owner)  {
             ZOLTAN_PRINT_INFO (dd->my_proc, yo, "Multiply defined GID");
             if (dd->debug_level > 4)
                ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);
             return ZOLTAN_WARN;
          }

          ptr->errcheck = owner;
          if (dd->debug_level > 5)
             ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);
          return ZOLTAN_OK;          /* ignore all errors */
       }
   }

   /* gid not found. Create new DD_Node and fill it in */
   nodeidx = DD_Memory_Alloc_Node(dd);
   ptr = dd->nodelist + nodeidx;

   ZOLTAN_SET_ID (dd->gid_length, ptr->gid, gid);

   if (lid) {
      ZOLTAN_SET_ID(dd->lid_length,ptr->gid + dd->gid_length, lid);
   }
   else  {
      memset(ptr->gid + dd->gid_length, 0,
             dd->lid_length*sizeof(ZOLTAN_ID_TYPE));
   }
   if (user) {
      memcpy(ptr->gid + (dd->gid_length + dd->lid_length), user,
             dd->user_data_length);
   }
   else {
      memset(ptr->gid + (dd->gid_length+dd->lid_length), 0,
             dd->user_data_length);
   }
   ptr->partition = partition;
   ptr->owner = owner;
   ptr->errcheck = owner;

   /* Add node to the linked list */
   ptr->next = dd->table[index];
   dd->table[index] = nodeidx;

   if (dd->debug_level > 6)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "Created new directory item");
   if (dd->debug_level > 5)
      ZOLTAN_TRACE_OUT (dd->my_proc, yo, NULL);

   return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
