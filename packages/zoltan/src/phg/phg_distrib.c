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

  if (maxgno<=0) {
     printf("ERROR: maxgno=%d nProc_dim=%d dist_dim=[%d,%d,%d]\n", maxgno, nProc_dim, dist_dim[0], dist_dim[1], dist_dim[2]);
     return 0;
  }
   
  idx = gno * nProc_dim / maxgno;

   if (!dist_dim) {
      printf("ERRROR: distdim is NULL\n");
      return 0;
   }
  // printf("UMIT: gno=%d idx = %d\n", gno, idx);
  while (gno < dist_dim[idx]) {
    idx--;
    if (idx<0) {
       printf("ERROR: idx<0");
       return 0;
    }
       
    }
  if (idx+1>nProc_dim) {
      printf("ERRROR idx+1(%d) > nProc_dim(%d)\n", idx+1, nProc_dim);
      return 0;
  }
  while (gno >= dist_dim[idx+1]) { 
      idx++;
   if (idx+1>nProc_dim) {
      printf("ERRROR idx+1(%d) > nProc_dim(%d)\n", idx+1, nProc_dim);
      return 0;
   } 
    }

  return idx;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
