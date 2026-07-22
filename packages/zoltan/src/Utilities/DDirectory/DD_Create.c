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
#include "zoltan_align.h"
#include "zz_hash.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*  NOTE: See file, README, for associated documentation. (RTH) */

/* Zoltan_DD_Create assumes the global object being managed by the
 * directory is a Zoltan global ID, which is a ZOLTAN_ID_TYPE-tuple.  
 * The "num_gid" parameter is where we specify how many ZOLTAN_ID_TYPEs 
 * are in the object (zz->Num_GID).
 *
 * However, some Zoltan code uses a data directory to manage other
 * global integer values.  When the ZOLTAN_ID_TYPE was always be an
 * unsigned integer, that worked.  
 *
 * But now that the ZOLTAN_ID_TYPE can be specified at compile time,
 * we need to be more careful.
 *
 * If the global value is not a Zoltan global ID, then the "num_gid"
 * parameter should be factor which, when multiplied by a ZOLTAN_ID_TYPE,
 * give an object of the same length as the global value.
 *
 * So if ZOLTAN_ID_TYPE is a 32-bit int, and the global value being
 * managed by the data directory is a 64-bit int, "num_gid" should
 * be "2".
 *
 * The "num_lid" parameter specifies the number ZOLTAN_ID_TYPEs in 
 * the local ID.
 *
 * The "user_length" parameter specifies the number of chars in the
 * user data.
 */

/*******************  Zoltan_DD_Create()  ***************************/

int Zoltan_DD_Create (
 Zoltan_DD_Directory **dd,    /* contains directory state and pointers */
 MPI_Comm comm,               /* Dup'ed and saved for future use       */
 int num_gid,                 /* Number of entries in a global ID.     */
 int num_lid,                 /* Number of entries in a local ID.      
                                 If zero, ignore LIDs                  */
 int user_length,             /* Optional user data length in chars, 0 ignore */
 int table_length,            /* sizeof hash table, use default if 0   */
 int debug_level              /* control actions to errors, normally 0 */
)
{
   int size, i;
   int my_proc;
   int array[3], max_array[3], min_array[3];
   char *yo = "Zoltan_DD_Create";

   if (MPI_Comm_rank(comm, &my_proc) != MPI_SUCCESS)  {
      ZOLTAN_PRINT_ERROR (-1, yo, "MPI_Comm_rank failed");
      return ZOLTAN_FATAL;
   }

   /* input sanity check */
   if (dd == NULL || num_gid < 1 || table_length < 0 || num_lid < 0)  {
      ZOLTAN_PRINT_ERROR (my_proc, yo, "Invalid input argument");
      return ZOLTAN_FATAL;
   }
   if (debug_level > 4)
      ZOLTAN_TRACE_IN (my_proc, yo, NULL);

   /* insure all processors are using the same GID, LID, USER lengths */
   array[0] = num_gid;
   array[1] = num_lid;
   array[2] = user_length;
   MPI_Allreduce (array, max_array, 3, MPI_INT, MPI_MAX, comm);
   MPI_Allreduce (array, min_array, 3, MPI_INT, MPI_MIN, comm);
   if (max_array[0] != min_array[0] || max_array[1] != min_array[1]
    || max_array[2] != min_array[2])  {
       ZOLTAN_PRINT_ERROR(-1,yo,"LID, GID, USER data lengths differ globally");
       return ZOLTAN_FATAL;
   }

   /* malloc memory for the directory structure + hash table */
   size = (table_length) ? table_length: ZOLTAN_DD_HASH_TABLE_COUNT;
   size = Zoltan_Recommended_Hash_Size(size);
   *dd  = (Zoltan_DD_Directory*) ZOLTAN_MALLOC (sizeof (Zoltan_DD_Directory)
        + size * sizeof(DD_NodeIdx));
   if (*dd == NULL)  {
      ZOLTAN_PRINT_ERROR (my_proc, yo, "Can not malloc hash table");
      if (debug_level > 4)
        ZOLTAN_TRACE_OUT(my_proc, yo, NULL);
      return ZOLTAN_MEMERR;
   }

   /* NULL heads of link list in hash table */
   for (i = 0; i < size; i++) (*dd)->table[i] = -1;  /* NULL values */
   (*dd)->nodecnt = 0;
   (*dd)->nodelist = NULL;
   (*dd)->nodedata = NULL;
   (*dd)->nodelistlen = 0;
   (*dd)->nextfreenode = -1;

   /* save useful constants into directory for convenience */
   (*dd)->debug_level      = debug_level;  /* [0,3], default 0          */
   (*dd)->gid_length       = num_gid;      /* saved input Num_GID       */
   (*dd)->lid_length       = num_lid;      /* saved input Num_LIB       */
   (*dd)->table_length     = size;         /* # of linked list heads    */
   (*dd)->user_data_length = user_length;  /* optional user data length */
   (*dd)->hash             = Zoltan_DD_Hash2;/* default hash algorithm   */
   (*dd)->hashdata         = NULL;         /* no hash data */
   (*dd)->hashfn           = NULL;         /* no hash function */
   (*dd)->cleanup          = NULL;         /* user registered cleanup   */
   (*dd)->max_id_length    = (num_gid > num_lid) ? num_gid : num_lid;

   /* frequently used dynamic allocation computed sizes */
   size = ((num_gid + num_lid) * sizeof(ZOLTAN_ID_TYPE)) + user_length;
   (*dd)->nodedata_size   = size;
   (*dd)->update_msg_size = size + sizeof(DD_Update_Msg);

   size = num_gid * sizeof(ZOLTAN_ID_TYPE);
   (*dd)->remove_msg_size = size + sizeof(DD_Remove_Msg);

   size = user_length + ((*dd)->max_id_length * sizeof(ZOLTAN_ID_TYPE));
   (*dd)->find_msg_size   = size + sizeof(DD_Find_Msg);

   /* force alignment */
   (*dd)->nodedata_size   = Zoltan_Align_size_t((*dd)->nodedata_size);
   (*dd)->update_msg_size = Zoltan_Align_size_t((*dd)->update_msg_size);
   (*dd)->remove_msg_size = Zoltan_Align_size_t((*dd)->remove_msg_size);
   (*dd)->find_msg_size   = Zoltan_Align_size_t((*dd)->find_msg_size);

   /* duplicate MPI comm to prevent future comm changes from disrupting  */
   /* directory communications & save the associated comm size & rank    */
   if (MPI_Comm_dup  (comm,  &((*dd)->comm))    != MPI_SUCCESS
    || MPI_Comm_size (comm,  &((*dd)->nproc))   != MPI_SUCCESS
    || MPI_Comm_rank (comm,  &((*dd)->my_proc)) != MPI_SUCCESS)  {
         ZOLTAN_PRINT_ERROR (my_proc, yo, "MPI Problem, unable to continue");
         return ZOLTAN_FATAL;
   }

   if (debug_level > 4)
      ZOLTAN_TRACE_OUT (my_proc, yo, NULL);
   return ZOLTAN_OK;
}

/*******************  Copy functions  ***************************/
    

Zoltan_DD_Directory *Zoltan_DD_Copy(Zoltan_DD_Directory *from)
{
  Zoltan_DD_Directory *to = NULL;

  Zoltan_DD_Copy_To(&to, from);

  return to;
}


int Zoltan_DD_Copy_To(Zoltan_DD_Directory **toptr, Zoltan_DD_Directory *from)
{
  static char *yo = "Zoltan_DD_Copy_To";
  Zoltan_DD_Directory *to= NULL;

  if (!toptr) {
    return ZOLTAN_FATAL;
  }

  if (*toptr) {
    Zoltan_DD_Destroy(toptr);
  }

  if (from) {
    DD_NodeIdx i;

    to = *toptr = 
      (Zoltan_DD_Directory *)ZOLTAN_MALLOC(
        sizeof (Zoltan_DD_Directory) + 
        (from->table_length * sizeof(DD_NodeIdx)));

    if (!to) {
      ZOLTAN_PRINT_ERROR(from->my_proc, yo, "Insufficient memory."); 
      return ZOLTAN_MEMERR;
    }
  
    *to = *from;
    memcpy(to->table, from->table, to->table_length * sizeof(DD_NodeIdx));

    MPI_Comm_dup(from->comm, &(to->comm));

    if (to->nodelistlen) {
      to->nodelist = (DD_Node *) ZOLTAN_MALLOC(to->nodelistlen * sizeof(DD_Node));
      memcpy(to->nodelist, from->nodelist, to->nodelistlen * sizeof(DD_Node));

      to->nodedata = (char *) ZOLTAN_MALLOC(to->nodelistlen * to->nodedata_size);
      memcpy(to->nodedata, from->nodedata, to->nodelistlen * to->nodedata_size);

      for (i = 0; i < to->nodelistlen; i++) {
        to->nodelist[i].gid = (ZOLTAN_ID_PTR)(to->nodedata + i*to->nodedata_size);
      }
    }
  }

  return ZOLTAN_OK;
}
  
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
