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
#ifndef __ZOLTAN_PHG_TREE_H
#define __ZOLTAN_PHG_TREE_H

#include "zz_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#define SET_MIN_NODE(ptr, offset, val) (ptr)[2*(offset)]=-(val)
#define SET_MAX_NODE(ptr, offset, val) (ptr)[2*(offset)+1]=(val)

/* Create tree structure */
int
Zoltan_PHG_create_tree(int **ptr, int part_number, int* tree_size);

int
Zoltan_PHG_centralize_tree(ZZ *zz, int p, int tree_size);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_PHG_COMM_H */

