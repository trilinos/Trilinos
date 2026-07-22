// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef _DR_PAR_UTIL_CONST_H_
#define _DR_PAR_UTIL_CONST_H_

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


extern
void print_sync_start (
  int proc,
  int do_print_line
);

extern
void print_sync_end (
  int proc,
  int nprocs,
  int do_print_line
);

extern
void boundary_exchange(
  MESH_INFO_PTR mesh,
  int vec_len,
  int *send_vec,
  int *recv_vec
);
/* Function prototypes */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
#endif /* _DR_PAR_UTIL_CONST_H_ */
