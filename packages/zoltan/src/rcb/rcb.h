// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


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

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
