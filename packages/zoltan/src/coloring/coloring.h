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

#ifndef __COLORING_H
#define __COLORING_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "coloring_const.h"

/* Metis also has a swap */
#ifdef SWAP
#undef SWAP
#endif

#define SWAP(a,b) tmp=(a);(a)=(b);(b)=tmp;

/* Macros for error handling */
#define ZOLTAN_COLOR_ERROR(error,str) {ierr = error ; \
 ZOLTAN_PRINT_ERROR(zz->Proc, yo, str) ; goto End ;}

#ifdef __cplusplus
}
#endif
#endif
