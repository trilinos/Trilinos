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

#include <stdio.h>
#include "zoltan.h"

#undef DEBUG_PRINT

int Zoltan_Order_Test(
  struct Zoltan_Struct *zz,
  int *num_gid_entries,
  int *num_lid_entries,
  int num_obj,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int *rank,
  int *iperm
)
{
  int i;
  int nbr;
  int first, last;
  int tmp;
#ifdef DEBUG_PRINT
  int ancestor;
#endif /* DEBUG_PRINT */
  ZOLTAN_ID_PTR perm_gids;

  nbr = Zoltan_Order_Get_Num_Blocks(zz);
#ifdef DEBUG_PRINT
  fprintf(stderr, "Nbr block : %d\n", nbr);
  fprintf(stderr,"details : \n");
#endif /* DEBUG_PRINT */

  tmp = 0;
  for (i = 0 ; i<nbr ; ++i) {
    int size;

    Zoltan_Order_Get_Block_Bounds(zz, i, &first, &last);
#ifdef DEBUG_PRINT
    fprintf(stderr, "block %i : %d to %d\n", i, first, last);
#endif /* DEBUG_PRINT */
    if (first != tmp) {
      fprintf(stderr, "Error, non consecutive numbering\n");
      return (ZOLTAN_FATAL);
    }
    tmp = last;
    size = Zoltan_Order_Get_Block_Size(zz, i);
    if (size != last - first) {
      fprintf(stderr, "Error, size doesn't match\n");
      return (ZOLTAN_FATAL);
    }

#ifdef DEBUG_PRINT
    ancestor = Zoltan_Order_Get_Block_Parent(zz, i);
    fprintf(stderr, "Father of %d : %d\n", i, ancestor);
#endif /* DEBUG_PRINT */
  }

  {
    int nbrleaves;
    int *leaves;
    int *blocks;

    nbrleaves = Zoltan_Order_Get_Num_Leaves(zz);
#ifdef DEBUG_PRINT
    fprintf(stderr, "Number of leaves : %d\n", nbrleaves);
#endif /* DEBUG_PRINT */

    blocks = (int*)ZOLTAN_MALLOC(nbr*sizeof(int));
    for (i = 0 ; i < nbr ; ++i)
      blocks[i] = 0;

    leaves = (int*)ZOLTAN_MALLOC((nbrleaves+1)*sizeof(int));
    Zoltan_Order_Get_Block_Leaves(zz, leaves);
    for (i = 0 ; i<nbrleaves ; ++i) {
      int parent;
#ifdef DEBUG_PRINT
      fprintf(stderr, "Leaf %d : %d\n", i, leaves[i]);
#endif /* DEBUG_PRINT */
      if (leaves[i] == -1) {
	fprintf(stderr, "Error, leavessize doesn't match\n");
	return (ZOLTAN_FATAL);
      }

      parent = leaves[i];
      do {
	blocks[parent] = 1;
	parent = Zoltan_Order_Get_Block_Parent(zz, parent);
      } while (parent != -1);
    }
    if (leaves[nbrleaves] != -1) {
      fprintf(stderr, "Error, leaves array is not valid\n");
      return (ZOLTAN_FATAL);
    }

    for (i = 0 ; i <nbr ; ++i) {
      if (blocks[i] != 1) {
	fprintf(stderr, "Error, not all nodes in array are in the tree\n");
	return (ZOLTAN_FATAL);
      }
    }
    ZOLTAN_FREE(&leaves);
    ZOLTAN_FREE(&blocks);
  }

  perm_gids = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC(num_obj *sizeof(int) * (*num_gid_entries));

  Zoltan_Order_Get_GID_Order(zz,global_ids, perm_gids);
/*   for (i = 0 ; i <num_obj ; ++i) { */
/*     fprintf(stderr, "%d(%d) becomes %d (%d)\n", global_ids[i], i, perm_gids[i], rank[i]); */
/*   } */


  ZOLTAN_FREE(&perm_gids);
  return (ZOLTAN_OK);
}


#ifdef __cplusplus
}
#endif
