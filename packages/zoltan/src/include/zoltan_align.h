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

#ifndef __ZOLTAN_ALIGN_H
#define __ZOLTAN_ALIGN_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* 
 * Plauger alignment algorithm, The Standard C Library.
 * Forces malloc'ed variable size struct alignment.
 * ZOLTAN_ALIGN_VAL is defined in Zoltan/include/zoltan_align.h;
 * values are 0,1,3,7U depending upon machine.
 * (E.g., 7U == 8 byte alignment.)
 */

#ifndef ZOLTAN_ALIGN_VAL
#define ZOLTAN_ALIGN_VAL 7U
#endif

extern int Zoltan_Align(int);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
