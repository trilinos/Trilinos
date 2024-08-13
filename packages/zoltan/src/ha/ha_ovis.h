// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __HA_OVIS_H
#define __HA_OVIS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#ifdef ZOLTAN_OVIS

#include "ovis.h"

/* Structure for storing OVIS parameters */

struct OVIS_parameters {
  int outputLevel;
  char hello[MAX_PARAM_STRING_LEN];
  char dll[MAX_PARAM_STRING_LEN];
  double minVersion;
};

extern int Zoltan_OVIS_Setup(ZZ *, struct OVIS_parameters *);
extern int Zoltan_OVIS_Set_Param(char *, char *);

#endif /* ZOLTAN_OVIS */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
