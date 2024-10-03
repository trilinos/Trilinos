// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ZOLTAN_DD_H
#define ZOLTAN_DD_H

#include "zoltan_mem.h"
#include "zoltan_comm.h"
#include "zoltan_types.h"
#include "zoltan_id.h"
#include "zoltan_util.h"
#include "zoltan_dd.h"
#ifdef __cplusplus

/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#define ZOLTAN_DD_HASH_TABLE_COUNT 100000 /* default # of linked list heads */
#define ZOLTAN_DD_NO_PROC           -1    /* not a possible processor #     */

/* Tags for MPI communications.  These need unique values. Arbitrary */
#define ZOLTAN_DD_FIND_MSG_TAG     29137  /* needs 3 consecutive values */
#define ZOLTAN_DD_UPDATE_MSG_TAG   29140  /* needs 2 consecutive values */
#define ZOLTAN_DD_REMOVE_MSG_TAG   29142  /* needs 2 consecutive values */

#ifndef TRUE
#define FALSE (0)
#define TRUE (1)
#endif /* !TRUE */

/***********  Distributed Directory Function Prototypes ************/

unsigned int Zoltan_DD_Hash2(ZOLTAN_ID_PTR key, int num_id_entries,
 unsigned int n, void *hashdata, ZOLTAN_HASH_FN *fn);

void Zoltan_DD_default_cleanup(void *hashdata);

typedef unsigned int DD_Hash_fn(ZOLTAN_ID_PTR, int, unsigned int, void *,
                                ZOLTAN_HASH_FN *);
typedef void DD_Cleanup_fn(void*);



/************  Zoltan_DD_Directory, DD_Node  **********************/

typedef int DD_NodeIdx;   /* Index into dd->nodelist; 
                             must be a signed type as -1 indicates NULL */

/* The following structure, DD_Node, is the basic unit of the hash table's
 * linked list.  DD_Node contains the global ID(used as the table lookup
 * key) and other necessary information.  NOTE: DD_Node is malloc'd to
 * store gid, lid & user data beyond the struct's end.
*/

typedef struct DD_Node  {
  int              owner;      /* processor hosting global ID object    */
  int              partition;  /* Optional data                         */
  int              errcheck;   /* Error checking(inconsistent updates) */
  DD_NodeIdx       next;       /* index in nodelist of next DD_Node in 
                                  linked list or free node list */
  ZOLTAN_ID_TYPE  *gid;        /* gid used as key for update & lookup   */
                               /* lid starts at gid + dd->gid_length    */
                               /*(user) data starts at                 */
                               /* gid + dd->gid_length + dd->lid_length */
  int              free;       /* flag indicating whether node is
                                  free or used */
} DD_Node;


/* The directory structure, Zoltan_DD_Struct, is created by the call
 * to Zoltan_DD_Create(). It maintains the state information and storage
 * allocation for the distributed directory. Other state information may
 * be added in the future. This structure must be passed back to all other
 * distributed directory calls: Zoltan_DD_Update(), Zoltan_DD_Find(),
 * Zoltan_DD_Destroy(), Zoltan_DD_Set_Hash_Fn(), DD_Update_Local(),
 * DD_Find_Local(), DD_Remove_Local().  NOTE: Zoltan_DD_Struct is
 * malloc'd for storage beyond the struct's end to hold hash table.
*/

struct Zoltan_DD_Struct {
  int my_proc;            /* My identity in MPI Comm                */
  int nproc;              /* Number of processors in MPI Comm       */
  int gid_length;         /* = zz->Num_GID -- avoid needing Zoltan_Struct */
  int lid_length;         /* = zz->Num_LID -- avoid needing Zoltan_Struct */
  int max_id_length;      /* max(gid_length, lid_length)           */
  int user_data_length;   /* Optional user data length in chars */
  int table_length;       /* # of heads of linked lists             */
  size_t nodedata_size;   /* Malloc for GID & LID & user storage    */
  size_t find_msg_size;   /* Total allocation for DD_FIND_MSG       */
  size_t update_msg_size; /* Total allocation for DD_UPDATE_MSG     */
  size_t remove_msg_size; /* Total allocation for DD_REMOVE_MSG     */
  int debug_level;        /* Determines actions to multiple updates */

  DD_Hash_fn *hash;       /* Hash function used by this DD         */
  void *hashdata;         /* Either our hash data ...              */
  ZOLTAN_HASH_FN *hashfn; /* ... Or user's hash function, not both */
  DD_Cleanup_fn *cleanup; /* Functioned to free our hash data      */

  MPI_Comm comm;          /* Dup of original MPI Comm(KDD)         */
  DD_NodeIdx nodecnt;     /* Number of nodes used in nodelist      */
                          /* AKA Number of local entries in directory      */
  DD_Node *nodelist;      /* Memory for storing all nodes in the directory */
  char *nodedata;         /* Memory for storing all data in the directory  */
  DD_NodeIdx nodelistlen; /* Length of the nodelist. */
  DD_NodeIdx nextfreenode;/* Index of first free node in nodelist; 
                             -1 if no nodes are free */
  DD_NodeIdx table[1];    /* Hash table heads of the link lists     */
};


/*************** DD Communication Messages *********************/

/* Note: These message structures should become MPI datatypes(KDD)   */


typedef struct {           /* Only used by Zoltan_DD_Update()              */
   char lid_flag;          /* indicates if LID data are present            */
   char user_flag;         /* indicates if USER data are present           */
   char partition_flag;    /* indicates if optional partition data present */
   int owner;              /* range [0, nproc-1] or -1                     */
   int partition;
   ZOLTAN_ID_TYPE gid[1];   /* struct malloc'd to include gid & lid & user */
                            /* LID found at gid[dd->gid_length]            */
                            /* USER found at gid[dd->gid_length            */
                            /*  + dd->lid_length] if used                  */
} DD_Update_Msg;



/* A single structure can serve for both the request and its answer(in
 * DD_Find_Msg if its memory is sized to hold either a global ID or a
 * local ID.  On the request direction, the proc field is the
 * destination and the id field holds the global ID being located.  In
 * the return direction, the proc field holds the objects location and
 * the id field holds its corresponding local ID. The index field is
 * untouched & unused.
*/

typedef struct {           /* Only used by Zoltan_DD_Find()         */
   int        proc;        /* destination or location               */
   int        partition;
   int        index;       /* to put things back in order afterward */
   ZOLTAN_ID_TYPE id[1];   /* allocated as max(Num_GID, Num_LID)    */
                           /* + user data length                    */
} DD_Find_Msg;


typedef struct  {          /* Only used by Zoltan_DD_Remove()      */
   int        owner;       /* range [0, nproc-1] or -1             */
   ZOLTAN_ID_TYPE gid[1];  /* structure malloc'd to include gid    */
}  DD_Remove_Msg;



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif


#endif
