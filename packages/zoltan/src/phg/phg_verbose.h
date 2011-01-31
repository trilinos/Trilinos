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
#ifndef __ZOLTAN_PHG_VERBOSE_H
#define __ZOLTAN_PHG_VERBOSE_H

#include "phg.h"
#include "phg_lookup.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


void print_zoltan_pins(zoltan_pins *z, int me, int ewgt_dim);
void print_hypergraph(ZZ *zz, ZHG *zhg, int sumWeight);
void show_edges(char *s, ZZ *zz, int num_lists, int num_pins,
                ZOLTAN_ID_TYPE *edg_GID, int *row_ptr, ZOLTAN_ID_TYPE  *vtx_GID);
void debug_graph_to_hg(
  int nedges, ZOLTAN_ID_PTR egids, ZOLTAN_ID_PTR elids,
  int *esizes, float *ewgts, int npins,
  ZOLTAN_ID_PTR pins, int *pin_procs, int ewgtdim, int lenGID, int lenLID);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
