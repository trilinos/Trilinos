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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"

int Zoltan_PHG_Gno_To_Proc_Block(
  int gno,
  int *dist_dim,
  int nProc_dim
)
{
/* Function that locates a given global number gno within a distribution 
 * vector dist.
 * Works for both vtx and edges.
 * Takes an initial guess based on equal distribution of gno's.
 * Modifies guess based on actual distribution.
 */

int idx;
int maxgno = dist_dim[nProc_dim];

  idx = gno * nProc_dim / maxgno;

  while (gno < dist_dim[idx]) idx--;
  while (gno >= dist_dim[idx+1]) idx++;

  return idx;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
