/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __RCB_CONST_H
#define __RCB_CONST_H

#include "shared_const.h"

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
  int       countlo, counthi;   /* # of dots at that position */
  int       proclo, prochi;	/* unique proc who owns a nearest dot */
};

struct rcb_box {       	     /* bounding box */
  double    lo[3], hi[3];	/* xyz lo/hi bounds */
};

typedef struct RCB_Struct {
  LB_ID_PTR Global_IDs;      /* Pointer to array of global IDs; global ID of 
                                Dots[i] starts in Global_IDs[i*lb->Num_GID].
                                Because lb->Num_GID is determined at runtime,
                                this info is most easily stored, allocated and
                                reallocated separately from Dots. 
                                This array is NOT used if LB_Use_IDs returns
                                FALSE.   */
  LB_ID_PTR Local_IDs;       /* Pointer to array of local IDs; local ID of 
                                Dots[i] starts in Local_IDs[i*lb->Num_LID].
                                Because lb->Num_LID is determined at runtime,
                                this info is most easily stored, allocated and
                                reallocated separately from Dots. 
                                This array is NOT used if LB_Use_IDs returns
                                FALSE.   */
  struct Dot_Struct *Dots;     
  struct rcb_tree *Tree_Ptr;
  struct rcb_box *Box;
} RCB_STRUCT;

extern int LB_RCB_Build_Structure(LB *, int *, int *, int, int);
extern void LB_RCB_Free_Structure(LB *);
extern int LB_Set_RCB_Param(char *, char *);

#endif
