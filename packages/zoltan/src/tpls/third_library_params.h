// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __THIRD_LIBRARY_PARAMS_H
#define __THIRD_LIBRARY_PARAMS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "zz_util_const.h"
#include "params_const.h"

/**********  parameters structure used by PHG, ParMetis and Jostle **********/
static PARAM_VARS Graph_Package_params[] = {
        { "GRAPH_PACKAGE", NULL, "STRING", 0 },
        { "ORDER_TYPE", NULL, "STRING", 0 },
        { NULL, NULL, NULL, 0 } };


#ifdef __cplusplus
}
#endif

#endif
