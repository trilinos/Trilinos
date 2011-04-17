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
