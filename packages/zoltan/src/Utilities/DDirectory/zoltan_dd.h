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

#include "zoltan_types.h"
#include "mpi.h"

struct Zoltan_DD_Struct;
typedef struct Zoltan_DD_Struct Zoltan_DD_Directory;

/* The following are used as return value error codes */
/* These values should be made to conform with Zoltan Error Codes. KDD */
#define ZOLTAN_DD_NORMAL_RETURN         0
#define ZOLTAN_DD_INPUT_ERROR           137   /* arbitrary, but distinct */
#define ZOLTAN_DD_MEMORY_ERROR          138
#define ZOLTAN_DD_MPI_ERROR             139
#define ZOLTAN_DD_COMM_ERROR            140
#define ZOLTAN_DD_GID_ADDED             141
#define ZOLTAN_DD_GID_NOT_FOUND_ERROR   142
#define ZOLTAN_DD_GID_REDEFINED_ERROR   143


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

void Zoltan_DD_Stats (Zoltan_DD_Directory *dd) ;

int Zoltan_DD_Set_Neighbor_Hash_Fn1 (Zoltan_DD_Directory *dd, int size) ;

int Zoltan_DD_Set_Neighbor_Hash_Fn2 (Zoltan_DD_Directory *dd, int *proc,
 int *low, int *high, int count) ;

int Zoltan_DD_Print (Zoltan_DD_Directory *dd) ;

#endif
