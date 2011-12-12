/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000-2012, Sandia National Laboratories.                    *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/

#ifndef __HA_OVIS_H
#define __HA_OVIS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#ifdef ZOLTAN_OVIS

extern int Zoltan_OVIS_Setup(ZZ *);
extern int Zoltan_OVIS_Set_Param(char *, char *);

#endif /* ZOLTAN_OVIS */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
