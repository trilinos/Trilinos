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

#ifndef __LBI_CONST_H
#define __LBI_CONST_H

#include <mpi.h>
#include "zoltan.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#warning "WARNING:  Use of zoltan include file lbi_const.h is deprecated and will not be supported in Trilinos v11.  Update your code to include zoltan.h instead. "
#warning "WARNING:  Use of zoltan include file lbi_const.h is deprecated and will not be supported in Trilinos v11.  Update your code to include zoltan.h instead. "
#warning "WARNING:  Use of zoltan include file lbi_const.h is deprecated and will not be supported in Trilinos v11.  Update your code to include zoltan.h instead. "

/****************************************************************************
 *  This file is maintained for backward compatability with previous versions
 *  of Zoltan that used LB_* for function and variable names.
 *  New Zoltan users should include zoltan.h instead of this file and use
 *  the currently supported Zoltan interface described there.
 */

/****************************************************************************
 *  Data type for global and local identifiers used in Zoltan.
 */

typedef ZOLTAN_ID_TYPE  LB_ID_TYPE;
typedef ZOLTAN_ID_TYPE *LB_ID_PTR;
#define LB_ID_MPI_TYPE  ZOLTAN_ID_MPI_TYPE

/****************************************************************************
 *  Function types for callback functions; used in LB_Set_Fn.
 */

#define LB_Fn_Type                         Zoltan_Fn_Type
#define LB_NUM_EDGES_FN_TYPE               ZOLTAN_NUM_EDGES_FN_TYPE
#define LB_EDGE_LIST_FN_TYPE               ZOLTAN_EDGE_LIST_FN_TYPE
#define LB_NUM_GEOM_FN_TYPE                ZOLTAN_NUM_GEOM_FN_TYPE
#define LB_GEOM_FN_TYPE                    ZOLTAN_GEOM_FN_TYPE
#define LB_NUM_OBJ_FN_TYPE                 ZOLTAN_NUM_OBJ_FN_TYPE
#define LB_OBJ_LIST_FN_TYPE                ZOLTAN_OBJ_LIST_FN_TYPE
#define LB_FIRST_OBJ_FN_TYPE               ZOLTAN_FIRST_OBJ_FN_TYPE
#define LB_NEXT_OBJ_FN_TYPE                ZOLTAN_NEXT_OBJ_FN_TYPE
#define LB_NUM_BORDER_OBJ_FN_TYPE          ZOLTAN_NUM_BORDER_OBJ_FN_TYPE
#define LB_BORDER_OBJ_LIST_FN_TYPE         ZOLTAN_BORDER_OBJ_LIST_FN_TYPE
#define LB_FIRST_BORDER_OBJ_FN_TYPE        ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE
#define LB_NEXT_BORDER_OBJ_FN_TYPE         ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE
#define LB_PRE_MIGRATE_FN_TYPE             ZOLTAN_PRE_MIGRATE_FN_TYPE
#define LB_MID_MIGRATE_FN_TYPE             ZOLTAN_MID_MIGRATE_FN_TYPE
#define LB_POST_MIGRATE_FN_TYPE            ZOLTAN_POST_MIGRATE_FN_TYPE
#define LB_OBJ_SIZE_FN_TYPE                ZOLTAN_OBJ_SIZE_FN_TYPE
#define LB_PACK_OBJ_FN_TYPE                ZOLTAN_PACK_OBJ_FN_TYPE
#define LB_UNPACK_OBJ_FN_TYPE              ZOLTAN_UNPACK_OBJ_FN_TYPE
#define LB_NUM_COARSE_OBJ_FN_TYPE          ZOLTAN_NUM_COARSE_OBJ_FN_TYPE
#define LB_COARSE_OBJ_LIST_FN_TYPE         ZOLTAN_COARSE_OBJ_LIST_FN_TYPE
#define LB_FIRST_COARSE_OBJ_FN_TYPE        ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE
#define LB_NEXT_COARSE_OBJ_FN_TYPE         ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE
#define LB_NUM_CHILD_FN_TYPE               ZOLTAN_NUM_CHILD_FN_TYPE
#define LB_CHILD_LIST_FN_TYPE              ZOLTAN_CHILD_LIST_FN_TYPE
#define LB_CHILD_WEIGHT_FN_TYPE            ZOLTAN_CHILD_WEIGHT_FN_TYPE
#define LB_GET_PROCESSOR_NAME_FN_TYPE      ZOLTAN_PROC_NAME_FN_TYPE
#define LB_MAX_FN_TYPES                    ZOLTAN_MAX_FN_TYPES

typedef enum Zoltan_Fn_Type LB_FN_TYPE;

/****************************************************************************
 * Types of refinement for building a refinement tree.
 */

#define LB_Ref_Type      Zoltan_Ref_Type
#define LB_OTHER_REF     ZOLTAN_OTHER_REF
#define LB_IN_ORDER      ZOLTAN_IN_ORDER
#define LB_TRI_BISECT    ZOLTAN_TRI_BISECT
#define LB_QUAD_QUAD     ZOLTAN_QUAD_QUAD
#define LB_HEX3D_OCT     ZOLTAN_HEX3D_OCT

typedef enum Zoltan_Ref_Type LB_REF_TYPE;

/****************************************************************************
 *  Other common definitions:
 */

#define LB_Struct Zoltan_Struct

/****************************************************************************
 * Error codes for Zoltan library; defined in "zoltan_types.h".
 */

#define LB_OK     ZOLTAN_OK
#define LB_WARN   ZOLTAN_WARN
#define LB_FATAL  ZOLTAN_FATAL
#define LB_MEMERR ZOLTAN_MEMERR

/****************************************************************************
 * Callback functions.
 */

#define LB_NUM_EDGES_FN             ZOLTAN_NUM_EDGES_FN
#define LB_NUM_EDGES_FORT_FN        ZOLTAN_NUM_EDGES_FORT_FN

#define LB_EDGE_LIST_FN             ZOLTAN_EDGE_LIST_FN
#define LB_EDGE_LIST_FORT_FN        ZOLTAN_EDGE_LIST_FORT_FN

#define LB_NUM_GEOM_FN              ZOLTAN_NUM_GEOM_FN
#define LB_NUM_GEOM_FORT_FN         ZOLTAN_NUM_GEOM_FORT_FN

#define LB_GEOM_FN                  ZOLTAN_GEOM_FN
#define LB_GEOM_FORT_FN             ZOLTAN_GEOM_FORT_FN

#define LB_NUM_OBJ_FN               ZOLTAN_NUM_OBJ_FN
#define LB_NUM_OBJ_FORT_FN          ZOLTAN_NUM_OBJ_FORT_FN

#define LB_OBJ_LIST_FN              ZOLTAN_OBJ_LIST_FN
#define LB_OBJ_LIST_FORT_FN         ZOLTAN_OBJ_LIST_FORT_FN

#define LB_FIRST_OBJ_FN             ZOLTAN_FIRST_OBJ_FN
#define LB_FIRST_OBJ_FORT_FN        ZOLTAN_FIRST_OBJ_FORT_FN

#define LB_NEXT_OBJ_FN              ZOLTAN_NEXT_OBJ_FN
#define LB_NEXT_OBJ_FORT_FN         ZOLTAN_NEXT_OBJ_FORT_FN

#define LB_OBJ_SIZE_FN              ZOLTAN_OBJ_SIZE_FN
#define LB_OBJ_SIZE_FORT_FN         ZOLTAN_OBJ_SIZE_FORT_FN

#define LB_PRE_MIGRATE_FN           ZOLTAN_PRE_MIGRATE_FN
#define LB_PRE_MIGRATE_FORT_FN      ZOLTAN_PRE_MIGRATE_FORT_FN

#define LB_MID_MIGRATE_FN           ZOLTAN_MID_MIGRATE_FN
#define LB_MID_MIGRATE_FORT_FN      ZOLTAN_MID_MIGRATE_FORT_FN

#define LB_POST_MIGRATE_FN          ZOLTAN_POST_MIGRATE_FN
#define LB_POST_MIGRATE_FORT_FN     ZOLTAN_POST_MIGRATE_FORT_FN

#define LB_PACK_OBJ_FN              ZOLTAN_PACK_OBJ_FN
#define LB_PACK_OBJ_FORT_FN         ZOLTAN_PACK_OBJ_FORT_FN

#define LB_UNPACK_OBJ_FN            ZOLTAN_UNPACK_OBJ_FN
#define LB_UNPACK_OBJ_FORT_FN       ZOLTAN_UNPACK_OBJ_FORT_FN

#define LB_GET_PROCESSOR_NAME_FN    ZOLTAN_PROC_NAME_FN

#define LB_NUM_COARSE_OBJ_FN        ZOLTAN_NUM_COARSE_OBJ_FN
#define LB_NUM_COARSE_OBJ_FORT_FN   ZOLTAN_NUM_COARSE_OBJ_FORT_FN

#define LB_COARSE_OBJ_LIST_FN       ZOLTAN_COARSE_OBJ_LIST_FN
#define LB_COARSE_OBJ_LIST_FORT_FN  ZOLTAN_COARSE_OBJ_LIST_FORT_FN

#define LB_FIRST_COARSE_OBJ_FN      ZOLTAN_FIRST_COARSE_OBJ_FN
#define LB_FIRST_COARSE_OBJ_FORT_FN ZOLTAN_FIRST_COARSE_OBJ_FORT_FN

#define LB_NEXT_COARSE_OBJ_FN       ZOLTAN_NEXT_COARSE_OBJ_FN
#define LB_NEXT_COARSE_OBJ_FORT_FN  ZOLTAN_NEXT_COARSE_OBJ_FORT_FN

#define LB_NUM_CHILD_FN             ZOLTAN_NUM_CHILD_FN
#define LB_NUM_CHILD_FORT_FN        ZOLTAN_NUM_CHILD_FORT_FN

#define LB_CHILD_LIST_FN            ZOLTAN_CHILD_LIST_FN
#define LB_CHILD_LIST_FORT_FN       ZOLTAN_CHILD_LIST_FORT_FN

#define LB_CHILD_WEIGHT_FN          ZOLTAN_CHILD_WEIGHT_FN
#define LB_CHILD_WEIGHT_FORT_FN     ZOLTAN_CHILD_WEIGHT_FORT_FN

#define LB_NUM_BORDER_OBJ_FN        ZOLTAN_NUM_BORDER_OBJ_FN
#define LB_NUM_BORDER_OBJ_FORT_FN   ZOLTAN_NUM_BORDER_OBJ_FORT_FN

#define LB_BORDER_OBJ_LIST_FN       ZOLTAN_BORDER_OBJ_LIST_FN
#define LB_BORDER_OBJ_LIST_FORT_FN  ZOLTAN_BORDER_OBJ_LIST_FORT_FN

#define LB_FIRST_BORDER_OBJ_FN      ZOLTAN_FIRST_BORDER_OBJ_FN
#define LB_FIRST_BORDER_OBJ_FORT_FN ZOLTAN_FIRST_BORDER_OBJ_FORT_FN

#define LB_NEXT_BORDER_OBJ_FN       ZOLTAN_NEXT_BORDER_OBJ_FN
#define LB_NEXT_BORDER_OBJ_FORT_FN  ZOLTAN_NEXT_BORDER_OBJ_FORT_FN


/****************************************************************************
 *  Zoltan functions 
 */

#define LB_Initialize               Zoltan_Initialize
#define LB_Create                   Zoltan_Create
#define LB_Destroy                  Zoltan_Destroy
#define LB_Set_Param                Zoltan_Set_Param

#define LB_Set_Fn                   Zoltan_Set_Fn
#define LB_Set_Num_Edges_Fn         Zoltan_Set_Num_Edges_Fn
#define LB_Set_Edge_List_Fn         Zoltan_Set_Edge_List_Fn
#define LB_Set_Num_Geom_Fn          Zoltan_Set_Num_Geom_Fn
#define LB_Set_Geom_Fn              Zoltan_Set_Geom_Fn
#define LB_Set_Num_Obj_Fn           Zoltan_Set_Num_Obj_Fn
#define LB_Set_Obj_List_Fn          Zoltan_Set_Obj_List_Fn
#define LB_Set_First_Obj_Fn         Zoltan_Set_First_Obj_Fn
#define LB_Set_Next_Obj_Fn          Zoltan_Set_Next_Obj_Fn
#define LB_Set_Num_Border_Obj_Fn    Zoltan_Set_Num_Border_Obj_Fn
#define LB_Set_Border_Obj_List_Fn   Zoltan_Set_Border_Obj_List_Fn
#define LB_Set_First_Border_Obj_Fn  Zoltan_Set_First_Border_Obj_Fn
#define LB_Set_Next_Border_Obj_Fn   Zoltan_Set_Next_Border_Obj_Fn
#define LB_Set_Pre_Migrate_Fn       Zoltan_Set_Pre_Migrate_Fn
#define LB_Set_Mid_Migrate_Fn       Zoltan_Set_Mid_Migrate_Fn
#define LB_Set_Post_Migrate_Fn      Zoltan_Set_Post_Migrate_Fn
#define LB_Set_Obj_Size_Fn          Zoltan_Set_Obj_Size_Fn
#define LB_Set_Pack_Obj_Fn          Zoltan_Set_Pack_Obj_Fn
#define LB_Set_Unpack_Obj_Fn        Zoltan_Set_Unpack_Obj_Fn
#define LB_Set_Num_Coarse_Obj_Fn    Zoltan_Set_Num_Coarse_Obj_Fn
#define LB_Set_Coarse_Obj_List_Fn   Zoltan_Set_Coarse_Obj_List_Fn
#define LB_Set_First_Coarse_Obj_Fn  Zoltan_Set_First_Coarse_Obj_Fn
#define LB_Set_Next_Coarse_Obj_Fn   Zoltan_Set_Next_Coarse_Obj_Fn
#define LB_Set_Num_Child_Fn         Zoltan_Set_Num_Child_Fn
#define LB_Set_Child_List_Fn        Zoltan_Set_Child_List_Fn
#define LB_Set_Child_Weight_Fn      Zoltan_Set_Child_Weight_Fn

#define LB_Set_Method(a,b)          Zoltan_Set_Param(a,"LB_METHOD",b)
#define LB_Balance                  Zoltan_LB_Balance
#define LB_Free_Data                Zoltan_LB_Free_Data
#define LB_Eval                     Zoltan_LB_Eval 

#define LB_Compute_Destinations     Zoltan_Compute_Destinations
#define LB_Help_Migrate             Zoltan_Help_Migrate


#define LB_Point_Assign             Zoltan_LB_Point_Assign
#define LB_Box_Assign               Zoltan_LB_Box_Assign

/*******  Memory management package ********/

#define LB_Memory_Stats             Zoltan_Memory_Stats
#define LB_Malloc                   Zoltan_Malloc
#define LB_Realloc                  Zoltan_Realloc
#define LB_Array_Alloc              Zoltan_Array_Alloc
#define LB_Free                     Zoltan_Free
#define LB_Multifree                Zoltan_Multifree
#define LB_Set_Memory_Debug         Zoltan_Memory_Debug
#define LB_Memory_Num               Zoltan_Memory_Num
#define LB_MALLOC                   ZOLTAN_MALLOC
#define LB_REALLOC                  ZOLTAN_REALLOC
#define LB_FREE                     ZOLTAN_FREE

/*******  Unstructured communication package ********/

#define Comm_Obj                Zoltan_Comm_Obj
#define LB_Comm_Create          Zoltan_Comm_Create
#define LB_Comm_Destroy         Zoltan_Comm_Destroy
#define LB_Comm_Invert_Map      Zoltan_Comm_Invert_Map
#define LB_Comm_Sort_Ints       Zoltan_Comm_Sort_Ints
#define LB_Comm_Exchange_Sizes  Zoltan_Comm_Exchange_Sizes
#define LB_Comm_Resize          Zoltan_Comm_Resize
#define LB_Comm_Do              Zoltan_Comm_Do
#define LB_Comm_Do_Reverse      Zoltan_Comm_Do_Reverse

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* !__LBI_CONST_H */
