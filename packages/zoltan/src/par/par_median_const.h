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


#ifndef __PAR_MEDIAN_CONST_H
#define __PAR_MEDIAN_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <mpi.h>

extern int Zoltan_RB_find_median(int, double *, double *, int *,
  int, int, double, MPI_Comm,
  double *, int, int *, int, int,
  int, int, int, double, double, double,
  double *, double *, int *, int, int);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
