// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef _DR_LOADBAL_CONST_H_
#define _DR_LOADBAL_CONST_H_

#include "dr_input_const.h"

#ifdef __cplusplus

  #include "zoltan_cpp.h"

  #define ZOLTAN_STRUCT Zoltan &

  extern "C" {

#else

#define ZOLTAN_STRUCT struct Zoltan_Struct *

#endif

extern int setup_zoltan(ZOLTAN_STRUCT, int, PROB_INFO_PTR, MESH_INFO_PTR,
                                PARIO_INFO_PTR); 
extern void setup_fixed_obj(MESH_INFO_PTR, int); 

extern int run_zoltan(ZOLTAN_STRUCT, int, PROB_INFO_PTR, MESH_INFO_PTR,
                      PARIO_INFO_PTR); 


extern int migrate_elements(int, MESH_INFO_PTR, ZOLTAN_STRUCT,
                            int, int, 
                            int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *,
                            int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, int *);

extern ELEM_INFO *search_by_global_id(MESH_INFO *, ZOLTAN_ID_TYPE, int *);

extern ZOLTAN_OBJ_SIZE_FN migrate_elem_size;
extern ZOLTAN_OBJ_SIZE_MULTI_FN migrate_elem_size_multi;

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* _DR_LOADBAL_CONST_H_ */
