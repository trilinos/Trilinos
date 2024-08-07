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
int Zoltan_LB_Box_Assign (
 ZZ *zz,
 double xlo,
 double ylo,
 double zlo,
 double xhi,
 double yhi,
 double zhi,
 int *procs,
 int *count)
{
  char *yo = "Zoltan_LB_Box_Assign";
  int tmp = 0;

  if (zz->LB.Box_Assign == NULL) {
    /* function not supported by current decomposition method */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Box_Assign not supported by chosen partitioning method.");
    return ZOLTAN_FATAL;  
  }

  if (zz->LB.PartDist != NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Non-uniform distribution of partitions over processors is specified; "
      "use Zoltan_LB_Box_PP_Assign.");
    return ZOLTAN_FATAL;
  }

  /* Call appropriate method.  Pass procs and count in partition arguments
   * for greater efficiency in LB.Box_Assign (Zoltan is partition-based.) */
  return zz->LB.Box_Assign(zz, xlo, ylo, zlo, xhi, yhi, zhi, NULL, &tmp, 
                           procs, count);
}

/****************************************************************************/
int Zoltan_LB_Box_PP_Assign (
 ZZ *zz,
 double xlo,
 double ylo,
 double zlo,
 double xhi,
 double yhi,
 double zhi,
 int *procs,
 int *proc_count,
 int *parts,
 int *part_count)
{
  char *yo = "Zoltan_LB_Box_PP_Assign";

  if (zz->LB.Box_Assign == NULL) {
    /* function not supported by current decomposition method */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Box_Assign not supported by chosen partitioning method.");
    return ZOLTAN_FATAL;  
  }

  /* Call appropriate method.  Pass procs and count in partition arguments
   * for greater efficiency in LB.Box_Assign (Zoltan is partition-based.) */
  return zz->LB.Box_Assign(zz, xlo, ylo, zlo, xhi, yhi, zhi, procs, proc_count,
                           parts, part_count);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
