/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
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
  if (zz->LB.Box_Assign == NULL) {
    /* function not supported by current decomposition method */
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                   "Box_Assign not supported by chosen partitioning method.");
     return ZOLTAN_FATAL ;  
  }

  return zz->LB.Box_Assign (zz, xlo, ylo, zlo, xhi, yhi, zhi, procs, count) ;
  }
