// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __RCB_PARAMS_H
#define __RCB_PARAMS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "params_const.h"

/*  Parameters structure for RCB method.  Used in  */
/*  Zoltan_RCB_Set_Param and Zoltan_RCB.                   */
static PARAM_VARS RCB_params[] = {
                  { "RCB_OVERALLOC", NULL, "DOUBLE", 0 },
                  { "RCB_REUSE", NULL, "INT", 0 },
                  { "CHECK_GEOM", NULL, "INT", 0 },
                  { "RCB_OUTPUT_LEVEL", NULL, "INT", 0 },
                  { "KEEP_CUTS", NULL, "INT", 0 },
                  { "RCB_LOCK_DIRECTIONS", NULL, "INT", 0 },
                  { "RCB_SET_DIRECTIONS", NULL, "INT", 0 },
                  { "RCB_RECTILINEAR_BLOCKS", NULL, "INT", 0 },
                  { "OBJ_WEIGHTS_COMPARABLE", NULL, "INT", 0 },
                  { "RCB_MULTICRITERIA_NORM", NULL, "INT", 0 },
                  { "RCB_MAX_ASPECT_RATIO", NULL, "DOUBLE", 0 },
                  { "AVERAGE_CUTS", NULL, "INT", 0 },
                  { "RANDOM_PIVOTS", NULL, "INT", 0 },
                  { "RCB_RECOMPUTE_BOX", NULL, "INT", 0 },
                  { "REDUCE_DIMENSIONS", NULL, "INT", 0 },
                  { "DEGENERATE_RATIO", NULL, "DOUBLE", 0 },
                  {"FINAL_OUTPUT",      NULL,  "INT",    0},
                  { NULL, NULL, NULL, 0 } };

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
