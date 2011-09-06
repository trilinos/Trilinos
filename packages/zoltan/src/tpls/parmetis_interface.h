/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifndef __PARMETIS_INTERFACE_H
#define __PARMETIS_INTERFACE_H

#include <limits.h>
#include "zoltan_comm.h"
#include "third_library_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Guess the version number of ParMetis if not defined */
/* PARMETIS_MAJOR_VERSION is only defined in version 3.0 and higher */
#if (!defined(PARMETIS_MAJOR_VERSION))
#define PARMETIS_MAJOR_VERSION 0
#define PARMETIS_MINOR_VERSION 0
#endif

/* ParMetis option defs. These must be identical to the defs
 * in defs.h in the version of ParMetis you are using!
 * Both ParMetis 2.0 and 3.0 defs are included below.
 */
#define OPTION_IPART            1
#define OPTION_FOLDF            2
#define OPTION_DBGLVL           3
#define PMV3_OPTION_DBGLVL      1
#define PMV3_OPTION_SEED        2
#define PMV3_OPTION_IPART       3
#define PMV3_OPTION_PSR         3
#define PMV3_OPT_USE_OBJ_SIZE   9  /* Added by EB, not in ParMetis */
#define MAX_OPTIONS             10 /* Max number of options +1 */
/* Other ParMetis constants we may need */
#define GLOBAL_DBGLVL		0  /* Default debug level */
#define GLOBAL_SEED		15 /* Default random seed */
#define COUPLED                 1  /* Processors coupled to partitions? */
#define DISCOUPLED              2  /* Processors coupled to partitions? */


int Zoltan_ParMetis(
  ZZ *, float *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  int **, int **, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  int **, int **);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
