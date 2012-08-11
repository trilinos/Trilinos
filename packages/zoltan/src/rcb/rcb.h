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


#ifndef __RCB_H
#define __RCB_H

#include "shared.h"
#include "rcb_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Data structures for parallel RCB */

struct rcb_tree {	     /* tree of RCB cuts */
  double    cut;        	/* position of cut */
  int       dim;	        /* dimension (012) of cut */
  int       parent;             /* parent of this node in cut tree */
  int       left_leaf;          /* left child of this node in cut tree */
  int       right_leaf;         /* right child of this node in cut tree */
};

struct rcb_median {          /* RCB cut info */
  double    totallo, totalhi;   /* weight in each half of active partition */
  double    valuelo, valuehi;	/* position of dot(s) nearest to cut */
  double    wtlo, wthi;         /* total weight of dot(s) at that position */
  ZOLTAN_GNO_TYPE countlo, counthi;   /* # of dots at that position */
  int       proclo, prochi;	/* unique proc who owns a nearest dot */
};

struct rcb_box {       	     /* bounding box */
  double    lo[3], hi[3];	/* xyz lo/hi bounds */
};

typedef struct RCB_Struct {
  ZOLTAN_ID_PTR Global_IDs;  /* Pointer to array of global IDs; global ID of 
                                Dots[i] starts in Global_IDs[i*zz->Num_GID].
                                Because zz->Num_GID is determined at runtime,
                                this info is most easily stored, allocated and
                                reallocated separately from Dots. 
                                This array is NOT used if Zoltan_RB_Use_IDs 
                                returns FALSE.   */
  ZOLTAN_ID_PTR Local_IDs;   /* Pointer to array of local IDs; local ID of 
                                Dots[i] starts in Local_IDs[i*zz->Num_LID].
                                Because zz->Num_LID is determined at runtime,
                                this info is most easily stored, allocated and
                                reallocated separately from Dots. 
                                This array is NOT used if Zoltan_RB_Use_IDs 
                                returns FALSE.   */
  struct Dot_Struct Dots;       /* coordinates, weights, etc */
  struct rcb_tree *Tree_Ptr;
  struct rcb_box *Box;
  int Num_Dim;    /* Number of dimensions in the input geometry. */
  ZZ_Transform Tran;        /* transformation for degenerate geometry */
} RCB_STRUCT;

extern int Zoltan_RCB_Build_Structure(ZZ *, int *, int *, int, double, int,int);
extern void Zoltan_RCB_Print_Structure(ZZ *zz, int howMany);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
