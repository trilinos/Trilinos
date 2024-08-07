// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"


/****************************************************************************/
int Zoltan_LB_Point_Assign (
 ZZ *zz,
 double *x,
 int *proc)
{
/* Returns processor to which a point should be assigned. */
  char *yo = "Zoltan_LB_Point_Assign";
  if (zz->LB.Point_Assign == NULL) {
    /* function not supported by current decomposition method */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Point_Assign not supported by chosen partitioning method.");
    return ZOLTAN_FATAL;  
  }

  if (zz->LB.PartDist != NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "Non-uniform distribution of partitions over processors is specified; "
      "use Zoltan_LB_Point_PP_Assign.");
    return ZOLTAN_FATAL;
  }

  /* call appropriate method; pass proc in partition argument for greater
   * efficiency within LB.Point_Assign (Zoltan is partition-based). */
  return zz->LB.Point_Assign(zz, x, NULL, proc); 
}

/****************************************************************************/
int Zoltan_LB_Point_PP_Assign (
 ZZ *zz,
 double *x,
 int *proc,
 int *part)
{
/* Returns processor and partition to which a point should be assigned. */
  char *yo = "Zoltan_LB_Point_PP_Assign";
  if (zz->LB.Point_Assign == NULL) {
    /* function not supported by current decomposition method */
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Point_Assign not supported by chosen partitioning method.");
     return ZOLTAN_FATAL ;   
  }

  return zz->LB.Point_Assign(zz, x, proc, part);  /* call appropriate method */
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
