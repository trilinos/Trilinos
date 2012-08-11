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
