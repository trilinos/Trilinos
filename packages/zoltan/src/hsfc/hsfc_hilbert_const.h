/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __HSFC_HILBERT_CONST_H
#define __HSFC_HILBERT_CONST_H

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
*    state tables modified to match existing HSFC direction (2-d),
*    no binary search, greatly simplifing the test for intersection,
*    backtracking does not explicitly work backwards up the tree,
*    query region and search region intersection calculated fresh at each level,
*    partitioning rather than database application.
*    (Note, no code examples were available in Lawder's publications, so this
*     is an original adaptation.)
*
****************************************************************************/

#include "zz_const.h"

static unsigned const int IMAX = ~(0U);



static unsigned const int data2d[] =  /* 2 dimension to nkey conversion */
 {0, 3, 1, 2,
  0, 1, 3, 2,
  2, 3, 1, 0,
  2, 1, 3, 0};

static unsigned const int state2d[] = /* 2 dimension to nkey state transitions */
 {1, 2, 0, 0,
  0, 1, 3, 1,
  2, 0, 2, 3,
  3, 3, 1, 2};

static unsigned const int idata2d[] = /* nkey to 2 dimension conversion */
 {0, 2, 3, 1,
  0, 1, 3, 2,
  3, 2, 0, 1,
  3, 1, 0, 2};

static unsigned const int istate2d[] = /* nkey to 2 dimension state transitions */
 {1, 0, 0, 2,
  0, 1, 1, 3,
  3, 2, 2, 0,
  2, 3, 3, 1};

  


static unsigned const int data3d[] =  /* 3d to nkey conversion */
 {0, 1, 3, 2, 7, 6, 4, 5,
  0, 7, 1, 6, 3, 4, 2, 5,
  0, 3, 7, 4, 1, 2, 6, 5,
  2, 3, 1, 0, 5, 4, 6, 7,
  4, 3, 5, 2, 7, 0, 6, 1,
  6, 5, 1, 2, 7, 4, 0, 3,
  4, 7, 3, 0, 5, 6, 2, 1,
  6, 7, 5, 4, 1, 0, 2, 3,
  2, 5, 3, 4, 1, 6, 0, 7,
  2, 1, 5, 6, 3, 0, 4, 7,
  4, 5, 7, 6, 3, 2, 0, 1,
  6, 1, 7, 0, 5, 2, 4, 3};

static unsigned const int state3d[] = /* 3d to nkey state transitions */
 {1,  2,  3,  2,  4, 5,  3,  5,
  2,  6,  0,  7,  8, 8,  0,  7,
  0,  9, 10,  9,  1, 1, 11, 11,
  6,  0,  6, 11,  9, 0,  9,  8,
 11, 11,  0,  7,  5, 9,  0,  7,
  4,  4,  8,  8,  0, 6, 10,  6,
  5,  7,  5,  3,  1, 1, 11, 11,
  6,  1,  6, 10,  9, 4,  9, 10,
 10,  3,  1,  1, 10, 3,  5,  9,
  4,  4,  8,  8,  2, 7,  2,  3,
  7,  2, 11,  2,  7, 5,  8,  5,
 10,  3,  2,  6, 10, 3,  4,  4};


static unsigned const int idata3d[] =  /* nkey to 3d conversion */
 {0, 1, 3, 2, 6, 7, 5, 4,
  0, 2, 6, 4, 5, 7, 3, 1,
  0, 4, 5, 1, 3, 7, 6, 2,
  3, 2, 0, 1, 5, 4, 6, 7,
  5, 7, 3, 1, 0, 2, 6, 4,
  6, 2, 3, 7, 5, 1, 0, 4,
  3, 7, 6, 2, 0, 4, 5, 1,
  5, 4, 6, 7, 3, 2, 0, 1,
  6, 4, 0, 2, 3, 1, 5, 7,
  5, 1, 0, 4, 6, 2, 3, 7,
  6, 7, 5, 4, 0, 1, 3, 2,
  3, 1, 5, 7, 6, 4, 0, 2};

static signed const int istate3d[] =   /* nkey to 3d state transitions */
 {1,  2,  2,  3,  3,  5,  5,  4,
  2,  0,  0,  8,  8,  7,  7,  6,
  0,  1,  1,  9,  9, 11, 11, 10,
 11,  6,  6,  0,  0,  9,  9,  8,
  9,  7,  7, 11, 11,  0,  0,  5,
 10,  8,  8,  6,  6,  4,  4,  0,
  3, 11, 11,  5,  5,  1,  1,  7,
  4,  9,  9, 10, 10,  6,  6,  1,
  5, 10, 10,  1,  1,  3,  3,  9,
  7,  4,  4,  2,  2,  8,  8,  3,
  8,  5,  5,  7,  7,  2,  2, 11,
  6,  3,  3,  4,  4, 10, 10,  2};

  

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
