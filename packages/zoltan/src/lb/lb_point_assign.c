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

#include "zz_const.h"


int Zoltan_LB_Point_Assign (
 ZZ *zz,
 double *x,
 int *proc)
  {
  char *yo = "Zoltan_LB_Point_Assign";
  if (zz->LB.Point_Assign == NULL) {
    /* function not supported by current decomposition method */
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Point_Assign not supported by chosen partitioning method.");
     return ZOLTAN_FATAL ;   
  }

  return zz->LB.Point_Assign (zz, x, proc) ;    /* call appropriate method */
  }
