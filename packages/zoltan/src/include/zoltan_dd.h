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
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING
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

#ifndef ZOLTAN_DD_DDIRECTORY_H
#define ZOLTAN_DD_DDIRECTORY_H

#include "zoltan_types.h"
#include <mpi.h>

/*
** Must define this function prototype before #ifdef __cplusplus
** to avoid warning when compiling with C++ on solaris
*/
typedef unsigned int ZOLTAN_HASH_FN(ZOLTAN_ID_PTR, int, unsigned int);

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

struct Zoltan_DD_Struct;

typedef struct Zoltan_DD_Struct Zoltan_DD_Directory;


/***********  Distributed Directory Function Prototypes ************/

int Zoltan_DD_Create(Zoltan_DD_Directory **dd, MPI_Comm comm, 
                     int num_gid, int num_lid, int user_length,
                     int table_length, int debug_level);

int Zoltan_DD_Copy_To(Zoltan_DD_Directory **toptr, Zoltan_DD_Directory *from);

Zoltan_DD_Directory *Zoltan_DD_Copy(Zoltan_DD_Directory *from);

void Zoltan_DD_Destroy(Zoltan_DD_Directory **dd);

int Zoltan_DD_Update(Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
                     ZOLTAN_ID_PTR lid, char *user, int *partition, int count);

int Zoltan_DD_Find(Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
                   ZOLTAN_ID_PTR lid, char *data, int *partition, int count,
                   int *owner);

int Zoltan_DD_Remove(Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
                     int count);

int Zoltan_DD_Set_Hash_Fn(Zoltan_DD_Directory *dd, ZOLTAN_HASH_FN *hash);

void Zoltan_DD_Stats(Zoltan_DD_Directory *dd);

int Zoltan_DD_Set_Neighbor_Hash_Fn1(Zoltan_DD_Directory *dd, int size);

int Zoltan_DD_Set_Neighbor_Hash_Fn2(Zoltan_DD_Directory *dd, int *proc,
                                    int *low, int *high, int count);

int Zoltan_DD_Set_Neighbor_Hash_Fn3(Zoltan_DD_Directory *dd, int total);

int Zoltan_DD_Print(Zoltan_DD_Directory *dd);

int Zoltan_DD_GetLocalKeys(Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR* gid, 
                           int* size);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
