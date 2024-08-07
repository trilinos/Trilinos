// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef CH_INIT_DIST_CONST_H
#define CH_INIT_DIST_CONST_H

#include "dr_const.h"
#include "dr_input_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern void ch_dist_init(int, int, PARIO_INFO_PTR, short **, int, MPI_Comm);
extern int ch_dist_num_vtx(int, short *);
extern int ch_dist_max_num_vtx(short *);
extern void ch_dist_vtx_list(int *, int*, int, short *);
extern int ch_dist_proc(int, short *, int);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif  /* CH_INIT_DIST_CONST_H */
