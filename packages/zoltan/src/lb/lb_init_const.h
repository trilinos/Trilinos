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


#ifndef __LB_INIT_CONST_H
#define __LB_INIT_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


extern void Zoltan_Migrate_Init(struct Zoltan_Migrate_Struct *);

extern void Zoltan_LB_Init(struct Zoltan_LB_Struct *, int);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
