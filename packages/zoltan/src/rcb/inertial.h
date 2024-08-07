// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __INERTIAL_H
#define __INERTIAL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* function prototypes */

extern void Zoltan_evals3(double[3][3], double *, double *, double *);
extern double Zoltan_determinant(double[3][3]);
extern void Zoltan_eigenvec3(double[3][3], double, double *, double *);

extern void Zoltan_evals2(double[2][2], double *, double *);
extern void Zoltan_eigenvec2(double[2][2], double, double *, double *);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
