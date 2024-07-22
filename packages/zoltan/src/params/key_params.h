// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __KEY_PARAMS_H
#define __KEY_PARAMS_H

#include "zz_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern int Zoltan_Set_Key_Param(ZZ *, const char *, const char *, int);
extern void Zoltan_Print_Key_Params(ZZ const *);
extern void Zoltan_Print_Configuration(char *indent);
extern int Zoltan_Filter_Params(ZZ *, ZZ *, PARAM_VARS *, int , int, int);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
