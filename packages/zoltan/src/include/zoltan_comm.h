// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __COMM_CONST_H
#define __COMM_CONST_H

#include <mpi.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Interface to the Zoltan Communication Package. */
/* This file should be included by user's of the  */
/* Communication package.                         */

struct Zoltan_Comm_Obj;
typedef struct Zoltan_Comm_Obj ZOLTAN_COMM_OBJ;

/* function prototypes */

MPI_Comm zoltan_get_global_comm();

int Zoltan_Comm_Create(ZOLTAN_COMM_OBJ**, int, int*, MPI_Comm, int, int*);

int Zoltan_Comm_Copy_To(ZOLTAN_COMM_OBJ **toptr, ZOLTAN_COMM_OBJ *from);

ZOLTAN_COMM_OBJ *Zoltan_Comm_Copy(ZOLTAN_COMM_OBJ *from);

int Zoltan_Comm_Destroy(ZOLTAN_COMM_OBJ**);

int Zoltan_Comm_Invert_Map(int*, int*, int, int, int**, int**, int*, int, int,
 int, int, MPI_Comm);

int Zoltan_Comm_Sort_Ints(int*, int*, int);

int Zoltan_Comm_Exchange_Sizes(int*, int*, int, int, int*, int*, int, int*, int,
 int, MPI_Comm);

int Zoltan_Comm_Resize(ZOLTAN_COMM_OBJ*, int*, int, int*);

int Zoltan_Comm_Do     (ZOLTAN_COMM_OBJ*, int, char*, int, char*);
int Zoltan_Comm_Do_Post(ZOLTAN_COMM_OBJ*, int, char*, int, char*);
int Zoltan_Comm_Do_Wait(ZOLTAN_COMM_OBJ*, int, char*, int, char*);
int Zoltan_Comm_Do_AlltoAll(ZOLTAN_COMM_OBJ*, char*, int, char*);

int Zoltan_Comm_Do_Reverse     (ZOLTAN_COMM_OBJ*, int, char*, int, int*, char*);
int Zoltan_Comm_Do_Reverse_Post(ZOLTAN_COMM_OBJ*, int, char*, int, int*, char*);
int Zoltan_Comm_Do_Reverse_Wait(ZOLTAN_COMM_OBJ*, int, char*, int, int*, char*);

int Zoltan_Comm_Info(ZOLTAN_COMM_OBJ*, int*, int*, int*, int*, int*, int*, int*,
 int*, int*, int*, int*, int*, int*);

int Zoltan_Comm_Invert_Plan(ZOLTAN_COMM_OBJ**);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
