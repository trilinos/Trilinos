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


#ifndef _ZOLTAN_HG_READFILE_CONST_H_
#define _ZOLTAN_HG_READFILE_CONST_H_

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdio.h>

/* Hypergraph read from file */
int Zoltan_HG_Readfile (int, FILE*, int*, int*, int*, int**, int**, int**, int*,
 float**, int*, float**, int*);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif 
