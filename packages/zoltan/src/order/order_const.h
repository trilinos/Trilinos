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

#ifndef __ORDER_CONST_H
#define __ORDER_CONST_H

#include "zoltan.h"

/*
 * Definition of the Zoltan Ordering Struct (ZOS).
 */

struct Zoltan_Order_Struct {
  ZZ *zz;                       /* ptr to Zoltan struct */
  int num_objects;              /* # of objects (local) */
  ZOLTAN_ID_PTR gids;           /* ptr to list of global ids */
  ZOLTAN_ID_PTR lids;           /* ptr to list of local ids */
  int *rank;        		/* rank[i] is the rank of gigs[i] */
  char *method;   		/* Ordering method used */
  int  num_separators;          /* Optional: # of separators. */
  int *sep_sizes;               /* Optional: Separator sizes. */
};

typedef struct Zoltan_Order_Struct ZOS;

/*
 * Type definitions for functions that depend on 
 * ordering method or uses the ordering struct.
 */

typedef int ZOLTAN_ORDER_FN(ZZ *, 
                         ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, 
                         int *, int *, int *,
                         ZOS *);

/*****************************************************************************/
/* PROTOTYPES */

/* Ordering functions */
extern ZOLTAN_ORDER_FN Zoltan_ParMetis_Order;

/* Parameter routine */
extern int Zoltan_Order_Set_Param(char *, char *);

/* Utility routines for permutations */
extern int Zoltan_Get_Distribution(ZZ *, int **);
extern int Zoltan_Inverse_Perm(ZZ *, int *, int *, int *, char *, int);

/*****************************************************************************/
/* Misc. constants */
#define RETURN_RANK  1
#define RETURN_IPERM 2

#endif
