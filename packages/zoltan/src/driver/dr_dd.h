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


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#ifndef __DR_DD_H
#define __DR_DD_H

extern int build_elem_dd(MESH_INFO_PTR);
extern int update_elem_dd(MESH_INFO_PTR);
extern int update_hvertex_proc(MESH_INFO_PTR);

#endif  /* __DR_DD_H */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
