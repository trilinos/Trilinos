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

#ifndef ZOLTAN_DD_DDIRECTORY_H
#define ZOLTAN_DD_DDIRECTORY_H

#include "mem_const.h"
#include "comm_const.h"
#include "zoltan_types.h"
#include "zoltan_id.h"
#include "zoltan_util.h"

#define ZOLTAN_DD_HASH_TABLE_COUNT  503   /* default # of linked list heads */
#define ZOLTAN_DD_NO_PROC           -1    /* not a possible processor #     */

/* The following are used as return value error codes */
#define ZOLTAN_DD_NORMAL_RETURN         0
#define ZOLTAN_DD_INPUT_ERROR           137   /* arbitrary, but distinct */
#define ZOLTAN_DD_MEMORY_ERROR          138
#define ZOLTAN_DD_MPI_ERROR             139
#define ZOLTAN_DD_COMM_ERROR            140
#define ZOLTAN_DD_GID_ADDED             141
#define ZOLTAN_DD_GID_NOT_FOUND_ERROR   142
#define ZOLTAN_DD_GID_REDEFINED_ERROR   143



/* Tags for MPI communications.  These need unique values. Arbitrary */
#define ZOLTAN_DD_FIND_MSG_TAG     29137  /* needs 3 consecutive values */
#define ZOLTAN_DD_UPDATE_MSG_TAG   29140  /* needs 2 consecutive values */
#define ZOLTAN_DD_REMOVE_MSG_TAG   29142  /* needs 2 consecutive values */


#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */




/************  Zoltan_DD_Directory, DD_Node  **********************/


/* The following structure, DD_Node, is the basic unit of the hash table's
// linked list.  DD_Node contains the global ID (used as the table lookup
// key) and other necessary information.  NOTE: DD_Node is malloc'd to
// store gid, lid & user data beyond the struct's end.
*/

typedef struct DD_Node
   {
   int              owner ;      /* processor hosting global ID object    */
   int              partition ;  /* Optional data                         */
   int              errcheck ;   /* Error checking (inconsistent updates) */
   struct DD_Node  *next ;       /* Next DD_Node in linked list or NULL   */
   ZOLTAN_ID_TYPE   gid[1] ;     /* gid used as key for update & lookup   */
                                 /* lid starts at gid + dd->gid_length    */
                                 /* (user) data starts at                 */
                                 /* gid + dd->gid_length + dd->lid_length */
   } DD_Node ;


/* The directory structure, Zoltan_DD_Directory, is created by the call
// to Zoltan_DD_Create(). It maintains the state information and storage
// allocation for the distributed directory. Other state information may
// be added in the future. This structure must be passed back to all other
// distributed directory calls: Zoltan_DD_Update(), Zoltan_DD_Find(),
// Zoltan_DD_Destroy(), Zoltan_DD_Set_Hash_Fn(), DD_Update_Local(),
// DD_Find_Local(), DD_Remove_Local().  NOTE: Zoltan_DD_Directory is
// malloc'd for storage beyond the struct's end to hold hash table.
*/

typedef struct
   {
   int my_proc ;            /* My identity in MPI Comm                */
   int nproc ;              /* Number of processors in MPI Comm       */
   int gid_length ;         /* = zz->Num_GID -- avoid needing Zoltan_Struct     */
   int lid_length ;         /* = zz->Num_LID -- avoid needing Zoltan_Struct     */
   int max_id_length ;      /* max (gid_length, lid_length)           */
   int user_data_length ;   /* Optional user data stored as ZOLTAN_ID_PTR */
   int table_length ;       /* # of heads of linked lists             */
   int node_size ;          /* Malloc'd to include GID & LID storage  */
   int find_msg_size ;      /* Total allocation for DD_FIND_MSG       */
   int update_msg_size ;    /* Total allocation for DD_UPDATE_MSG     */
   int remove_msg_size ;    /* Total allocation for DD_REMOVE_MSG     */
   int debug_level ;        /* Determines actions to multiple updates */

   unsigned int (*hash)(ZOLTAN_ID_PTR, int, unsigned int) ;
   void (*cleanup) (void) ;

   MPI_Comm comm ;          /* Dup of original MPI Comm (KDD)         */
   DD_Node *table[1] ;      /* Hash table, heads of the link lists    */
   } Zoltan_DD_Directory ;







/*************** DD Communication Messages *********************/

/* Note: These message structures should become MPI datatypes (KDD)   */


typedef struct           /* Only used by Zoltan_DD_Update()           */
   {
   int owner ;           /* range [0, nproc-1]                        */
   int partition ;
   ZOLTAN_ID_TYPE gid[1] ;   /* struct malloc'd to include gid & lid & user */
                             /* LID found at gid[dd->gid_length]            */
                             /* USER found at gid[dd->gid_length            */
                             /*  + dd->lid_length] if used                  */
   } DD_Update_Msg ;






/* A single structure can serve for both the request and its answer (in
// DD_Find_Msg if its memory is sized to hold either a global ID or a
// local ID.  On the request direction, the proc field is the
// destination and the id field holds the global ID being located.  In
// the return direction, the proc field holds the objects location and
// the id field holds its corresponding local ID. The index field is
// untouched & unused.
*/

typedef struct            /* Only used by Zoltan_DD_Find()         */
   {
   int        proc ;      /* destination or location               */
   int        partition ;
   int        index ;     /* to put things back in order afterward */
   ZOLTAN_ID_TYPE id[1] ; /* allocated as max(Num_GID, Num_LID)    */
                          /* + user data length                    */
   } DD_Find_Msg ;






typedef struct             /* Only used by Zoltan_DD_Remove()      */
   {
   int        owner ;      /* range [0, nproc-1]                   */
   ZOLTAN_ID_TYPE gid[1] ; /* structure malloc'd to include gid    */
   }  DD_Remove_Msg ;







/***********  Distributed Directory Function Prototypes ************/

int Zoltan_DD_Create (Zoltan_DD_Directory **dd, MPI_Comm comm, int num_gid,
 int num_lid, int user_length,  int table_length, int debug_level) ;


void Zoltan_DD_Destroy (Zoltan_DD_Directory **dd) ;

int Zoltan_DD_Update (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
 ZOLTAN_ID_PTR lid, ZOLTAN_ID_PTR user, int *partition, int count) ;

int Zoltan_DD_Find (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid, 
 ZOLTAN_ID_PTR lid, ZOLTAN_ID_PTR data, int *partition, int count,
 int *owner) ;

int Zoltan_DD_Remove (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
 int count) ;

int Zoltan_DD_Set_Hash_Fn (Zoltan_DD_Directory *dd,
 unsigned int (*hash) (ZOLTAN_ID_PTR, int, unsigned int)) ;

unsigned int DD_Hash2(ZOLTAN_ID_PTR key, int num_id_entries,
 unsigned int n) ;

void Zoltan_DD_Stats (Zoltan_DD_Directory *dd) ;


int Zoltan_DD_Set_Neighbor_Hash_Fn1 (Zoltan_DD_Directory *dd, int size) ;


int Zoltan_DD_Set_Neighbor_Hash_Fn2 (Zoltan_DD_Directory *dd, int *proc,
 int *low, int *high, int count) ;

int Zoltan_DD_Print (Zoltan_DD_Directory *dd) ;

#endif
