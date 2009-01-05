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

#include "zz_const.h"
#include "phg_tree.h"
#include "phg.h"
#include <limits.h>

  /* The tree is a list of couple (-min, max) but only declare as an array of int */

/* Create tree structure */
int
Zoltan_PHG_create_tree(int **ptr, int part_number, int* tree_size)
{
  int part2; /* Power of 2 parts */
  int i;

  if (part_number == 0)
    return ZOLTAN_OK;

  /* Round up to the next highest power of 2 */
  /* From http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2 */
  part2 = part_number;
  part2--;
  part2 |= part2 >> 1;
  part2 |= part2 >> 2;
  part2 |= part2 >> 4;
  part2 |= part2 >> 8;
  part2 |= part2 >> 16;
  part2 |= part2 >> 32; /* On 64 bits */
  part2++;

  *tree_size = 2*part2-1;
  *ptr = (int*) ZOLTAN_MALLOC(sizeof(int)*2*(*tree_size));
  if (*ptr == NULL)
    return ZOLTAN_MEMERR;
  /* TRICK: we store -low, thus we can use MPI_MAX for low and high */
  for (i = 0 ; i < *tree_size ; ++i) {
    SET_MIN_NODE(*ptr, i, *tree_size + 1);
    SET_MAX_NODE(*ptr, i, -1);
  }

  *ptr -=2; /* Begin at offset 1 */

  return ZOLTAN_OK;
}


/* Build a centralized tree */
int
Zoltan_PHG_centralize_tree(ZZ *zz, int p, int tree_size)
{
  /* TRICK: we store -low, thus we can use MPI_MAX for low and high */
  MPI_Allreduce(MPI_IN_PLACE, zz->LB.Tree + 2, 2*tree_size, MPI_INT, MPI_MAX, zz->Communicator);
  return ZOLTAN_OK;
}

/******************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
