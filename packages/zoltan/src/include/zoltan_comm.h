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
