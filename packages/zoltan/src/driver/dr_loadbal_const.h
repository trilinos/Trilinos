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


#ifndef _DR_LOADBAL_CONST_H_
#define _DR_LOADBAL_CONST_H_

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "dr_input_const.h"


extern int setup_zoltan(struct Zoltan_Struct *, int, PROB_INFO_PTR, MESH_INFO_PTR); 
extern int run_zoltan(struct Zoltan_Struct *, int, PROB_INFO_PTR, MESH_INFO_PTR,
                      PARIO_INFO_PTR); 
extern int migrate_elements(int, MESH_INFO_PTR, struct Zoltan_Struct *, 
                            int, int, 
                            int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *,
                            int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *);

extern ELEM_INFO *search_by_global_id(MESH_INFO *, int, int *);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* _DR_LOADBAL_CONST_H_ */
