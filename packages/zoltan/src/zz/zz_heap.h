// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ZOLTAN_HEAP_H
#define ZOLTAN_HEAP_H

#include "zz_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Heap datastructure */
typedef struct {
   int    space;
   int    n;
   int   *ele;
   int   *pos;
   float *value;
   } HEAP;

#define Zoltan_Heap_Empty(H)         (((H)->n)==0)
#define Zoltan_Heap_Not_Empty(H)     (((H)->n)!=0)
#define Zoltan_Heap_Has_Elem(H,e)    (((H)->pos[e])!=-1)
#define Zoltan_Heap_Value(H,e)       ((H)->value[e])    
#define Zoltan_Heap_Max_Value(H)     ((H)->value[(H)->ele[0]])
#define Zoltan_Heap_Peek_Max(H)      ((H)->ele[0])
#define Zoltan_Heap_Count(H)         ((H)->n)

int  Zoltan_Heap_Init         (ZZ*, HEAP*, int);
void Zoltan_Heap_Clear        (HEAP*);
void Zoltan_Heap_Free         (HEAP*);
int  Zoltan_Heap_Check        (HEAP*);
int  Zoltan_Heap_Input        (HEAP*, int, float);
int  Zoltan_Heap_Make         (HEAP*);
int  Zoltan_Heap_Change_Value (HEAP*, ZOLTAN_GNO_TYPE, float);
int  Zoltan_Heap_Extract_Max  (HEAP*);
int  Zoltan_Heap_Extract      (HEAP*, int);



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HEAP_H_ */
