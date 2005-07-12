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


#ifndef __ZOLTAN_UTIL_CONST_H
#define __ZOLTAN_UTIL_CONST_H

#include "zoltan_types.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

extern unsigned int Zoltan_Hash(ZOLTAN_ID_PTR, int, unsigned int);
extern int Zoltan_Clean_String(char *, char **);
void Zoltan_Transform_Box(double *lo, double *hi, double (*m)[3], int ndims);
void Zoltan_Transform_Box_Points(
  double *lo, double *hi, double (*m)[3], int ndims, double (*v)[3]);
void Zoltan_Box_Box_Intersection( double *lo1, double *hi1,
  double *lo2, double *hi2);
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
