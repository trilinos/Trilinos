/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __PHG_UTIL_H
#define __PHG_UTIL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdarg.h>
#include "phg_comm.h"

    /* UVC: some utility functions not particularly related to hypergraph */
extern char *uMe(PHGComm *);
extern void uprintf(PHGComm *, char *,...);
extern void errexit(char *,...);
    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __HG_UTIL_H */
