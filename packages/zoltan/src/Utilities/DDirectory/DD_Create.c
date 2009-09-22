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


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"
#include "zoltan_align.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

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
   int size;
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
   *dd  = (Zoltan_DD_Directory*) ZOLTAN_MALLOC (sizeof (Zoltan_DD_Directory)
        + size * sizeof(DD_Node*));
   if (*dd == NULL)  {
      ZOLTAN_PRINT_ERROR (my_proc, yo, "Can not malloc hash table");
      if (debug_level > 4)
        ZOLTAN_TRACE_OUT(my_proc, yo, NULL);
      return ZOLTAN_MEMERR;
   }

   /* NULL heads of link list in hash table */
   memset ((char*) (*dd)->table, '\0', size * sizeof(DD_Node*));

   /* save useful constants into directory for convenience */
   (*dd)->debug_level      = debug_level;  /* [0,3], default 0          */
   (*dd)->gid_length       = num_gid;      /* saved input Num_GID       */
   (*dd)->lid_length       = num_lid;      /* saved input Num_LIB       */
   (*dd)->table_length     = size;         /* # of linked list heads    */
   (*dd)->user_data_length = user_length;  /* optional user data length */
   (*dd)->hash             = Zoltan_DD_Hash2;/* default hash algorithm   */
   (*dd)->data             = NULL;         /* no data */
   (*dd)->cleanup          = NULL;         /* user registered cleanup   */
   (*dd)->max_id_length    = (num_gid > num_lid) ? num_gid : num_lid;

   /* frequently used dynamic allocation computed sizes */
   size = (num_gid + num_lid + user_length) * sizeof(ZOLTAN_ID_TYPE);
   (*dd)->node_size       = size + sizeof(DD_Node);
   (*dd)->update_msg_size = size + sizeof(DD_Update_Msg);

   size = num_gid * sizeof(ZOLTAN_ID_TYPE);
   (*dd)->remove_msg_size = size + sizeof(DD_Remove_Msg);

   size = (user_length + (*dd)->max_id_length) * sizeof(ZOLTAN_ID_TYPE);
   (*dd)->find_msg_size   = size + sizeof(DD_Find_Msg);

   /* force alignment */
   (*dd)->update_msg_size = Zoltan_Align((*dd)->update_msg_size);
   (*dd)->remove_msg_size = Zoltan_Align((*dd)->remove_msg_size);
   (*dd)->find_msg_size   = Zoltan_Align((*dd)->find_msg_size);

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
    
static void allocate_copy_list(DD_Node **new, DD_Node *l, int len);

Zoltan_DD_Directory *Zoltan_DD_Copy(Zoltan_DD_Directory *from)
{
  Zoltan_DD_Directory *to = NULL;

  Zoltan_DD_Copy_To(&to, from);

  return to;
}
int Zoltan_DD_Copy_To(Zoltan_DD_Directory **toptr, Zoltan_DD_Directory *from)
{
  static char *yo = "Zoltan_DD_Copy_To";
  int i, proc = 0;
  Zoltan_DD_Directory *to= NULL;

  if (!toptr){
    return ZOLTAN_FATAL;
  }

  if (*toptr){
    Zoltan_DD_Destroy(toptr);
  }

  if (from){
    proc = from->my_proc;

    to = *toptr = 
      (Zoltan_DD_Directory *)ZOLTAN_MALLOC(
        sizeof (Zoltan_DD_Directory) + 
        (from->table_length * sizeof(DD_Node*)));

    if (!to){
      ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory."); 
      return ZOLTAN_MEMERR;
    }
  
    *to = *from;

    MPI_Comm_dup(from->comm, &(to->comm));

    for (i=0; i<to->table_length; i++){
      allocate_copy_list(&(to->table[i]), from->table[i], 
                                to->node_size);
    }
  }

  return ZOLTAN_OK;
}
static void allocate_copy_list(DD_Node **new, DD_Node *l, int len)
{
  DD_Node *next = l;
  DD_Node *node = NULL, *prev = NULL;

  *new = NULL;

  if (len == 0){
    return ;
  }

  while (next){
    node = (DD_Node *)ZOLTAN_MALLOC(len);
    memcpy(node, next, len);
    if (prev){
      prev->next = node;
    }
    else{
      *new = node;
    }
    node->next = NULL;
    prev = node;
    next = next->next;
  }
}
  
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
