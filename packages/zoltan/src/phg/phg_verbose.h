// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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
