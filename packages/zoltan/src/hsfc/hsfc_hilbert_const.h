// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __HSFC_HILBERT_CONST_H
#define __HSFC_HILBERT_CONST_H

#include "zz_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/***************************************************************************
*  The modules and data arrays in hsfc_hilbert.c and hsfc_box_assign.c are 
*  influenced by a series of papers:
*  "Load balancing for the parallel adaptive solution of partial differential
*    equations", 1994, deCougny, Devine, Flaherty, Loy, Ozturan, Shephard
*  "Adaptive Local Refinement with Octree Load-Balancing for the parallel
*    solution of three-dimensional conservation laws", 1998, Loy, (Ph.D.)
*  "The Performance of an Octree Load Balancer for Parallel Adaptive Finite
*   Element Computation", 2001, Campbell
*  "Using State Diagrams for Hilbert Curve Mappings", Lawder, 2000
*  "Querying Multi-demensional Data Indexed Using the Hilbert Space-Filling
*    Curve", Lawder, King, 2000
*  "Calculation of Mappings between One and n-dimensional Values Using the
*    Hilbert Space-filling Curve", Lawder, 2000
*  "Using Space-filling Curves for Multi-dimensional Indexing", Lawder, King,
*    2000
*
*  A useful reference discussing the Hilbert curve for spatial data indexing,
*   data queries, etc. is "Fundamentals of Spatial Information Systems, 1992,
*   Laurini, Thompson.
*
*  Useful code examples for the generation of Hilbert and Inverse Hilbert 
*   coordinates came from Octree (in above papers) also using state tables,
*   H. Carter Edwards (1997), (Ph.D. dissertation, code copyrighted 1997),
*   Doug Moore, Rice University (copyrighted 1998-2000) whose code also
*    includes an implimenation of a query search for box like regions in the
*    Hilbert space as well as other database oriented functions.
*
*  This code is closest in spirit to Lawder, but with significant differences:
*    octree state tables modified to match existing HSFC direction
*     (Note, Lawder tables are incorrect!),
*    no binary search, greatly simplifing the test for intersection,
*    backtracking does not explicitly work backwards up the tree,
*    query region and search region intersection calculated fresh at each level,
*    partitioning rather than database application.
*    (Note, no code examples were available in Lawder's publications, so this
*     is an original adaptation.)
*
****************************************************************************/


static unsigned const int IMAX = ~(0U);



static unsigned const int idata2d[] =  /* 2 dimension to nkey conversion */
 {0, 3, 1, 2,
  0, 1, 3, 2,
  2, 3, 1, 0,
  2, 1, 3, 0};

static unsigned const int istate2d[] = /* 2 dimension to nkey state transitions */
 {1, 2, 0, 0,
  0, 1, 3, 1,
  2, 0, 2, 3,
  3, 3, 1, 2};

static unsigned const int data2d[] = /* nkey to 2 dimension conversion */
 {0, 2, 3, 1,
  0, 1, 3, 2,
  3, 2, 0, 1,
  3, 1, 0, 2};

static unsigned const int state2d[] = /* nkey to 2 dimension state transitions */
 {1, 0, 0, 2,
  0, 1, 1, 3,
  3, 2, 2, 0,
  2, 3, 3, 1};

static unsigned const data3d [] = {  /* nkey to 3 dimension conversion */
 0,  4,  6,  2,  3,  7,  5,  1,
 0,  1,  3,  2,  6,  7,  5,  4,
 0,  4,  5,  1,  3,  7,  6,  2,
 5,  4,  0,  1,  3,  2,  6,  7,
 6,  7,  3,  2,  0,  1,  5,  4,
 3,  7,  6,  2,  0,  4,  5,  1,
 5,  4,  6,  7,  3,  2,  0,  1,
 0,  1,  5,  4,  6,  7,  3,  2,
 5,  1,  0,  4,  6,  2,  3,  7,
 5,  1,  3,  7,  6,  2,  0,  4,
 0,  2,  6,  4,  5,  7,  3,  1,
 3,  1,  0,  2,  6,  4,  5,  7,
 5,  7,  6,  4,  0,  2,  3,  1,
 6,  7,  5,  4,  0,  1,  3,  2,
 3,  1,  5,  7,  6,  4,  0,  2,
 0,  2,  3,  1,  5,  7,  6,  4,
 3,  2,  0,  1,  5,  4,  6,  7,
 3,  2,  6,  7,  5,  4,  0,  1,
 6,  2,  0,  4,  5,  1,  3,  7,
 3,  7,  5,  1,  0,  4,  6,  2,
 5,  7,  3,  1,  0,  2,  6,  4,
 6,  2,  3,  7,  5,  1,  0,  4,
 6,  4,  0,  2,  3,  1,  5,  7,
 6,  4,  5,  7,  3,  1,  0,  2};

static unsigned const state3d [] = { /* nkey to 3 dimension state transitions */
    1,  2,  0,  3,  4,  0,  5,  6,
    0,  7,  1,  8,  5,  1,  4,  9,
   15,  0,  2, 22, 20,  2, 19, 23,
   20,  6,  3, 23, 15,  3, 16, 22,
   22, 13,  4, 12, 11,  4,  1, 20,
   11, 19,  5, 20, 22,  5,  0, 12,
    9,  3,  6,  2, 21,  6, 17,  0,
   10,  1,  7, 11, 12,  7, 13, 14,
   12,  9,  8, 14, 10,  8, 18, 11,
    6,  8,  9,  7, 17,  9, 21,  1,
    7, 15, 10, 16, 13, 10, 12, 17,
    5, 14, 11,  9,  0, 11, 22,  8,
    8, 20, 12, 19, 18, 12, 10,  5,
   18,  4, 13,  5,  8, 13,  7, 19,
   17, 11, 14,  1,  6, 14, 23,  7,
    2, 10, 15, 18, 19, 15, 20, 21,
   19, 17, 16, 21,  2, 16,  3, 18,
   14, 16, 17, 15, 23, 17,  6, 10,
   13, 21, 18, 17,  7, 18,  8, 16,
   16,  5, 19,  4,  3, 19,  2, 13,
    3, 12, 20, 13, 16, 20, 15,  4,
   23, 18, 21, 10, 14, 21,  9, 15,
    4, 23, 22,  6,  1, 22, 11,  3,
   21, 22, 23,  0,  9, 23, 14,  2};

static unsigned const idata3d [] = {   /* 3 dimension to nkey conversion */
 0,  7,  3,  4,  1,  6,  2,  5,
 0,  1,  3,  2,  7,  6,  4,  5,
 0,  3,  7,  4,  1,  2,  6,  5,
 2,  3,  5,  4,  1,  0,  6,  7,
 4,  5,  3,  2,  7,  6,  0,  1,
 4,  7,  3,  0,  5,  6,  2,  1,
 6,  7,  5,  4,  1,  0,  2,  3,
 0,  1,  7,  6,  3,  2,  4,  5,
 2,  1,  5,  6,  3,  0,  4,  7,
 6,  1,  5,  2,  7,  0,  4,  3,
 0,  7,  1,  6,  3,  4,  2,  5,
 2,  1,  3,  0,  5,  6,  4,  7,
 4,  7,  5,  6,  3,  0,  2,  1,
 4,  5,  7,  6,  3,  2,  0,  1,
 6,  1,  7,  0,  5,  2,  4,  3,
 0,  3,  1,  2,  7,  4,  6,  5,
 2,  3,  1,  0,  5,  4,  6,  7,
 6,  7,  1,  0,  5,  4,  2,  3,
 2,  5,  1,  6,  3,  4,  0,  7,
 4,  3,  7,  0,  5,  2,  6,  1,
 4,  3,  5,  2,  7,  0,  6,  1,
 6,  5,  1,  2,  7,  4,  0,  3,
 2,  5,  3,  4,  1,  6,  0,  7,
 6,  5,  7,  4,  1,  2,  0,  3};

static unsigned const istate3d [] ={ /* 3 dimension to nkey state transitions */
 1,  6,  3,  4,  2,  5,  0,  0,
 0,  7,  8,  1,  9,  4,  5,  1,
15, 22, 23, 20,  0,  2, 19,  2,
 3, 23,  3, 15,  6, 20, 16, 22,
11,  4, 12,  4, 20,  1, 22, 13,
22, 12, 20, 11,  5,  0,  5, 19,
17,  0,  6, 21,  3,  9,  6,  2,
10,  1, 14, 13, 11,  7, 12,  7,
 8,  9,  8, 18, 14, 12, 10, 11,
21,  8,  9,  9,  1,  6, 17,  7,
 7, 17, 15, 12, 16, 13, 10, 10,
11, 14,  9,  5, 11, 22,  0,  8,
18,  5, 12, 10, 19,  8, 12, 20,
 8, 13, 19,  7,  5, 13, 18,  4,
23, 11,  7, 17, 14, 14,  6,  1,
 2, 18, 10, 15, 21, 19, 20, 15,
16, 21, 17, 19, 16,  2,  3, 18,
 6, 10, 16, 14, 17, 23, 17, 15,
18, 18, 21,  8, 17,  7, 13, 16,
 3,  4, 13, 16, 19, 19,  2,  5,
16, 13, 20, 20,  4,  3, 15, 12,
 9, 21, 18, 21, 15, 14, 23, 10,
22, 22,  6,  1, 23, 11,  4,  3,
14, 23,  2,  9, 22, 23, 21,  0};



double Zoltan_HSFC_InvHilbert1d (ZZ*, double *coord);
double Zoltan_HSFC_InvHilbert2d (ZZ*, double *coord);
double Zoltan_HSFC_InvHilbert3d (ZZ*, double *coord);

void Zoltan_HSFC_Hilbert1d (ZZ*, double *coord, double key);
void Zoltan_HSFC_Hilbert2d (ZZ*, double *coord, double key);
void Zoltan_HSFC_Hilbert3d (ZZ*, double *coord, double key);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
