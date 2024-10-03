// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __GRAPH_PARAMS_H
#define __GRAPH_PARAMS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "zz_util_const.h"
#include "params_const.h"

/* Parameters for how to build the graph */
static PARAM_VARS ZG_params[] = {
	{ "GRAPH_SYMMETRIZE", NULL, "STRING", 0 },
	{ "GRAPH_SYM_WEIGHT", NULL, "STRING", 0 },
	{ "GRAPH_BIPARTITE_TYPE", NULL, "STRING", 0},
	{ "GRAPH_BUILD_TYPE", NULL, "STRING", 0},
	{ "GRAPH_FAST_BUILD_BASE", NULL, "INTEGER", 0},
	{ NULL, NULL, NULL, 0 } };


#ifdef __cplusplus
}
#endif

#endif
