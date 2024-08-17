// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __ORDER_CONST_H
#define __ORDER_CONST_H

#include "zoltan.h"
#ifndef __PARAMS_CONST_H
#include "params_const.h" /* needed for MAX_PARAM_STRING_LEN */
#endif
#include "zoltan_util.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*
 * Definition of the Zoltan Ordering Struct (ZOS).
 * This structure contains information about one particular ordering.
 */

struct Zoltan_Order_Struct {
  int needfree;
  int nbr_objects;              /* # of objects (local) */
  ZOLTAN_ID_PTR gids;           /* ptr to list of global ids */
  ZOLTAN_ID_PTR lids;           /* ptr to list of local ids */
  int *rank;		/* rank[i] is the rank of gids[i] */
  ZOLTAN_ID_PTR gidrank;
  int *iperm;
  int  start_index;
  char method[MAX_PARAM_STRING_LEN]; /* Ordering method used */
  char order_type[MAX_PARAM_STRING_LEN]; /* Ordering method used */

  /* Elimination Tree */
  int nbr_blocks;               /* Out: number of ordering blocks */
  int *start;                   /* Out: start[i] is the first vertex of block i */
  int *ancestor;                /* Out: father of block i */
  int *leaves;                  /* Out: list of all leaves */
  int nbr_leaves;               /* Number of leaves */

  int *vtxdist;                 /* How vertices are distributed accross processors */

  /* Deprecated */
  int  num_separators;          /* Optional: # of separators. */
  int *sep_sizes;               /* Optional: Separator sizes. */
};

typedef struct Zoltan_Order_Struct ZOS;

/*
 * Definition of Zoltan Order Option struct.
 * This structure contains options that are passed on to the ordering method.
 */

struct Zoltan_Order_Options {
  char method[MAX_PARAM_STRING_LEN];	   /* In: Ordering method. */
  int start_index;		/* In: Permutations start at 0 or 1? */
  int use_order_info;		/* In: Put order info into ZOS? */
  int return_args;		/* Out: What return arguments were computed? */
};

typedef struct Zoltan_Order_Options ZOOS;

/*
 * Type definitions for functions that depend on
 * ordering method or uses the ordering struct.
 */

typedef int ZOLTAN_ORDER_FN(  struct Zoltan_Struct *zz, int,
			      ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
			      int *, ZOOS *);

/*****************************************************************************/
/* PROTOTYPES */

/* Ordering functions */
extern ZOLTAN_ORDER_FN Zoltan_ParMetis_Order;

#ifdef ZOLTAN_SCOTCH
extern ZOLTAN_ORDER_FN Zoltan_Scotch_Order;
#endif /* ZOLTAN_SCOTCH */

#ifdef CEDRIC_2D_PARTITIONS
int Zoltan_HUND(
  struct Zoltan_Struct *zz,               /* Zoltan structure */
  int num_gid_entries, /* # of entries for a global id */
  int num_obj,		/* Number of objects to order */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
                        /* The application must allocate enough space */
  ZOLTAN_ID_PTR rank            /* rank[i] is the rank of gids[i] */
  /* int *iperm            /\* iperm[rank[i]]=i, only for sequential ordering *\/ */
  );
#endif /* CEDRIC_2D_PARTITIONS */

/* Parameter routine */
extern int Zoltan_Order_Set_Param(char *, char *);

/* Utility routines for permutations */
extern int Zoltan_Get_Distribution(  struct Zoltan_Struct *zz, int **);
extern int Zoltan_Inverse_Perm(  struct Zoltan_Struct *zz, int *, int *, int *, char *, int);
extern int Zoltan_Get_Processor_Graph(int *vtxdist, int Num_Proc, int i);


/* Utility routines for memory management */
extern int  Zoltan_Order_Init_Tree (struct Zoltan_Order_Struct *order, int blocknbr, int leavesnbr);
extern void Zoltan_Order_Free_Struct(struct Zoltan_Order_Struct *order);



/*****************************************************************************/
/* Misc. constants */
#define RETURN_RANK  1
#define RETURN_IPERM 2

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
