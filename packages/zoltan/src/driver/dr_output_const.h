// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef _DR_OUTPUT_CONST_H_
#define _DR_OUTPUT_CONST_H_

#include "dr_input_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern void print_distributed_mesh(
  int Proc,
  int Num_Proc,
  MESH_INFO_PTR mesh
);

extern int output_results(
  const char *cmd_file,
  const char *tag,
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh);

extern int output_gnu(
  const char *cmd_file,
  const char *tag,
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* _DR_OUTPUT_CONST_H_ */
