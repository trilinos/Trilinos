/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/

#ifndef __INERTIAL_H
#define __INERTIAL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* macros for routines */
#define max(a, b) ((a) < (b) ? (b) : (a))
#define min(a, b) ((a) > (b) ? (b) : (a))
#define sign(x)   ((x) >= 0 ? 1.0 : -1.0)

/* function prototypes */

extern void evals3(double[3][3], double *, double *, double *);
extern double determinant(double[3][3]);
extern void eigenvec3(double[3][3], double, double *, double *);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
