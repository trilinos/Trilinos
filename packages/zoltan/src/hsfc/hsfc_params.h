
#ifndef __HSFC_PARAMS_H
#define __HSFC_PARAMS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"

/* This structure is the Zoltan convention for user settable parameters */
static PARAM_VARS HSFC_params[] =
   {{"KEEP_CUTS", NULL, "INT", 0},
    { "REDUCE_DIMENSIONS", NULL, "INT", 0 },
    { "DEGENERATE_RATIO", NULL, "DOUBLE", 0 },
    {"FINAL_OUTPUT",  NULL,  "INT",    0},
    {NULL,        NULL,  NULL, 0}};


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
