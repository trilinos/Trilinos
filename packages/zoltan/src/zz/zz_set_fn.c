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


#include "zz_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines used to register callback functions.
 *  A generic routine (Zoltan_Set_Fn) is provided, as well as callback-function
 *  specific registration routines that enable tighter type-checking by
 *  the compiler.
 *
 *  When new callback functions are added to Zoltan, they should be 
 *  added to the case statement in Zoltan_Set_Fn and should have separate
 *  registration functions specific to the new callback types.
 */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_Set_Fn(ZZ *zz, ZOLTAN_FN_TYPE fn_type, void (*fn)(void), void *data)
{
/*
 *  Function to initialize a given LB interface function.
 *  Input:
 *    zz                --  Pointer to a Zoltan structure.
 *    fn_type           --  Enum type indicating the function to be set.
 *    fn                --  Pointer to the function to be used in the
 *                          assignment.
 *    data              --  Pointer to data that the LB library will
 *                          pass as an argument to fn(). May be NULL.
 *  Output:
 *    zz                --  Appropriate field set to value in void *().
 */

char *yo = "Zoltan_Set_Fn";
char msg[256];
int ierr;

  switch (fn_type) {
  case ZOLTAN_PARTITION_FN_TYPE:
    ierr = Zoltan_Set_Partition_Fn(zz, 
                  (ZOLTAN_PARTITION_FN *) fn, data);
    break;
  case ZOLTAN_PARTITION_MULTI_FN_TYPE:
    ierr = Zoltan_Set_Partition_Multi_Fn(zz, 
                  (ZOLTAN_PARTITION_MULTI_FN *) fn, data);
    break;
  case ZOLTAN_NUM_EDGES_FN_TYPE:
    ierr = Zoltan_Set_Num_Edges_Fn(zz, 
                  (ZOLTAN_NUM_EDGES_FN *) fn, data);
    break;
  case ZOLTAN_NUM_EDGES_MULTI_FN_TYPE:
    ierr = Zoltan_Set_Num_Edges_Multi_Fn(zz, 
                  (ZOLTAN_NUM_EDGES_MULTI_FN *) fn, data);
    break;
  case ZOLTAN_EDGE_LIST_FN_TYPE:
    ierr = Zoltan_Set_Edge_List_Fn(zz, 
                  (ZOLTAN_EDGE_LIST_FN *) fn, data);
    break;
  case ZOLTAN_EDGE_LIST_MULTI_FN_TYPE:
    ierr = Zoltan_Set_Edge_List_Multi_Fn(zz, 
                  (ZOLTAN_EDGE_LIST_MULTI_FN *) fn, data);
    break;
  case ZOLTAN_NUM_GEOM_FN_TYPE:
    ierr = Zoltan_Set_Num_Geom_Fn(zz, 
                  (ZOLTAN_NUM_GEOM_FN *) fn, data);
    break;
  case ZOLTAN_GEOM_MULTI_FN_TYPE:
    ierr = Zoltan_Set_Geom_Multi_Fn(zz, 
                  (ZOLTAN_GEOM_MULTI_FN *) fn, data);
    break;
  case ZOLTAN_GEOM_FN_TYPE:
    ierr = Zoltan_Set_Geom_Fn(zz, 
                  (ZOLTAN_GEOM_FN *) fn, data);
    break;
  case ZOLTAN_NUM_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Num_Obj_Fn(zz, 
                  (ZOLTAN_NUM_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_OBJ_LIST_FN_TYPE:
    ierr = Zoltan_Set_Obj_List_Fn(zz, 
                  (ZOLTAN_OBJ_LIST_FN *) fn, data);
    break;
  case ZOLTAN_FIRST_OBJ_FN_TYPE:
    ierr = Zoltan_Set_First_Obj_Fn(zz, 
                  (ZOLTAN_FIRST_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_NEXT_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Next_Obj_Fn(zz, 
                  (ZOLTAN_NEXT_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_NUM_BORDER_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Num_Border_Obj_Fn(zz, 
                  (ZOLTAN_NUM_BORDER_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_BORDER_OBJ_LIST_FN_TYPE:
    ierr = Zoltan_Set_Border_Obj_List_Fn(zz, 
                  (ZOLTAN_BORDER_OBJ_LIST_FN *) fn, data);
    break;
  case ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE:
    ierr = Zoltan_Set_First_Border_Obj_Fn(zz, 
                  (ZOLTAN_FIRST_BORDER_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Next_Border_Obj_Fn(zz, 
                  (ZOLTAN_NEXT_BORDER_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_PRE_MIGRATE_PP_FN_TYPE:
    ierr = Zoltan_Set_Pre_Migrate_PP_Fn(zz, 
                  (ZOLTAN_PRE_MIGRATE_PP_FN *) fn, data);
    break;
  case ZOLTAN_MID_MIGRATE_PP_FN_TYPE:
    ierr = Zoltan_Set_Mid_Migrate_PP_Fn(zz, 
                  (ZOLTAN_MID_MIGRATE_PP_FN *) fn, data);
    break;
  case ZOLTAN_POST_MIGRATE_PP_FN_TYPE:
    ierr = Zoltan_Set_Post_Migrate_PP_Fn(zz, 
                  (ZOLTAN_POST_MIGRATE_PP_FN *) fn, data);
    break;
  case ZOLTAN_PRE_MIGRATE_FN_TYPE:
    ierr = Zoltan_Set_Pre_Migrate_Fn(zz, 
                  (ZOLTAN_PRE_MIGRATE_FN *) fn, data);
    break;
  case ZOLTAN_MID_MIGRATE_FN_TYPE:
    ierr = Zoltan_Set_Mid_Migrate_Fn(zz, 
                  (ZOLTAN_MID_MIGRATE_FN *) fn, data);
    break;
  case ZOLTAN_POST_MIGRATE_FN_TYPE:
    ierr = Zoltan_Set_Post_Migrate_Fn(zz, 
                  (ZOLTAN_POST_MIGRATE_FN *) fn, data);
    break;
  case ZOLTAN_OBJ_SIZE_FN_TYPE:
    ierr = Zoltan_Set_Obj_Size_Fn(zz, 
                  (ZOLTAN_OBJ_SIZE_FN *) fn, data);
    break;
  case ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE:
    ierr = Zoltan_Set_Obj_Size_Multi_Fn(zz, 
                  (ZOLTAN_OBJ_SIZE_MULTI_FN *) fn, data);
    break;
  case ZOLTAN_PACK_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Pack_Obj_Fn(zz, 
                  (ZOLTAN_PACK_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_PACK_OBJ_MULTI_FN_TYPE:
    ierr = Zoltan_Set_Pack_Obj_Multi_Fn(zz, 
                  (ZOLTAN_PACK_OBJ_MULTI_FN *) fn, data);
    break;
  case ZOLTAN_UNPACK_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Unpack_Obj_Fn(zz, 
                  (ZOLTAN_UNPACK_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE:
    ierr = Zoltan_Set_Unpack_Obj_Multi_Fn(zz, 
                  (ZOLTAN_UNPACK_OBJ_MULTI_FN *) fn, data);
    break;
  case ZOLTAN_NUM_COARSE_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Num_Coarse_Obj_Fn(zz, 
                  (ZOLTAN_NUM_COARSE_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_COARSE_OBJ_LIST_FN_TYPE:
    ierr = Zoltan_Set_Coarse_Obj_List_Fn(zz, 
                  (ZOLTAN_COARSE_OBJ_LIST_FN *) fn, data);
    break;
  case ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE:
    ierr = Zoltan_Set_First_Coarse_Obj_Fn(zz, 
                  (ZOLTAN_FIRST_COARSE_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Next_Coarse_Obj_Fn(zz, 
                  (ZOLTAN_NEXT_COARSE_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_NUM_CHILD_FN_TYPE:
    ierr = Zoltan_Set_Num_Child_Fn(zz, 
                  (ZOLTAN_NUM_CHILD_FN *) fn, data);
    break;
  case ZOLTAN_CHILD_LIST_FN_TYPE:
    ierr = Zoltan_Set_Child_List_Fn(zz, 
                  (ZOLTAN_CHILD_LIST_FN *) fn, data);
    break;
  case ZOLTAN_CHILD_WEIGHT_FN_TYPE:
    ierr = Zoltan_Set_Child_Weight_Fn(zz, 
                  (ZOLTAN_CHILD_WEIGHT_FN *) fn, data);
    break;
  case ZOLTAN_HG_SIZE_CS_FN_TYPE:
    ierr = Zoltan_Set_HG_Size_CS_Fn(zz, 
                  (ZOLTAN_HG_SIZE_CS_FN *) fn, data);
    break;
  case ZOLTAN_HG_CS_FN_TYPE:
    ierr = Zoltan_Set_HG_CS_Fn(zz, 
                  (ZOLTAN_HG_CS_FN *) fn, data);
    break;
  case ZOLTAN_HG_SIZE_EDGE_WTS_FN_TYPE:
    ierr = Zoltan_Set_HG_Size_Edge_Wts_Fn(zz, 
                  (ZOLTAN_HG_SIZE_EDGE_WTS_FN *) fn, data);
    break;
  case ZOLTAN_HG_EDGE_WTS_FN_TYPE:
    ierr = Zoltan_Set_HG_Edge_Wts_Fn(zz, 
                  (ZOLTAN_HG_EDGE_WTS_FN *) fn, data);
    break;
  case ZOLTAN_NUM_FIXED_OBJ_FN_TYPE:
    ierr = Zoltan_Set_Num_Fixed_Obj_Fn(zz, 
                  (ZOLTAN_NUM_FIXED_OBJ_FN *) fn, data);
    break;
  case ZOLTAN_FIXED_OBJ_LIST_FN_TYPE:
    ierr = Zoltan_Set_Fixed_Obj_List_Fn(zz, 
                  (ZOLTAN_FIXED_OBJ_LIST_FN *) fn, data);
    break;
  case ZOLTAN_HIER_NUM_LEVELS_FN_TYPE:
    ierr = Zoltan_Set_Hier_Num_Levels_Fn(zz,
		  (ZOLTAN_HIER_NUM_LEVELS_FN *) fn, data);
    break;
  case ZOLTAN_HIER_PARTITION_FN_TYPE:
    ierr = Zoltan_Set_Hier_Partition_Fn(zz,
		  (ZOLTAN_HIER_PARTITION_FN *) fn, data);
    break;
  case ZOLTAN_HIER_METHOD_FN_TYPE:
    ierr = Zoltan_Set_Hier_Method_Fn(zz,
		  (ZOLTAN_HIER_METHOD_FN *) fn, data);
    break;
  default:
    sprintf(msg, "ZOLTAN_FN_TYPE %d is invalid.\n"
            "Value must be in range 0 to %d.", fn_type, ZOLTAN_MAX_FN_TYPES);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ierr = ZOLTAN_WARN;
  }

  return (ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * Callback registration functions that are specific to a given function
 * type.  Each callback function type must have a related registration 
 * function.  Using these functions enables more type-checking between the
 * application and Zoltan.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Set_Partition_Multi_Fn(
  ZZ *zz, 
  ZOLTAN_PARTITION_MULTI_FN *fn, 
  void *data
)
{
  zz->Get_Partition_Multi = fn;
  zz->Get_Partition_Multi_Data = data;
  zz->Get_Partition = NULL;
  zz->Get_Partition_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Partition_Fn(
  ZZ *zz, 
  ZOLTAN_PARTITION_FN *fn, 
  void *data
)
{
  zz->Get_Partition = fn;
  zz->Get_Partition_Data = data;
  zz->Get_Partition_Multi = NULL;
  zz->Get_Partition_Multi_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Num_Edges_Multi_Fn(
  ZZ *zz, 
  ZOLTAN_NUM_EDGES_MULTI_FN *fn, 
  void *data
)
{
  zz->Get_Num_Edges_Multi = fn;
  zz->Get_Num_Edges_Multi_Data = data;
  zz->Get_Num_Edges = NULL;
  zz->Get_Num_Edges_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Num_Edges_Fn(
  ZZ *zz, 
  ZOLTAN_NUM_EDGES_FN *fn, 
  void *data
)
{
  zz->Get_Num_Edges = fn;
  zz->Get_Num_Edges_Data = data;
  zz->Get_Num_Edges_Multi = NULL;
  zz->Get_Num_Edges_Multi_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Edge_List_Multi_Fn(
  ZZ *zz, 
  ZOLTAN_EDGE_LIST_MULTI_FN *fn, 
  void *data
)
{
  zz->Get_Edge_List_Multi = fn;
  zz->Get_Edge_List_Multi_Data = data;
  zz->Get_Edge_List = NULL;
  zz->Get_Edge_List_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Edge_List_Fn(
  ZZ *zz, 
  ZOLTAN_EDGE_LIST_FN *fn, 
  void *data
)
{
  zz->Get_Edge_List = fn;
  zz->Get_Edge_List_Data = data;
  zz->Get_Edge_List_Multi = NULL;
  zz->Get_Edge_List_Multi_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Num_Geom_Fn(
  ZZ *zz, 
  ZOLTAN_NUM_GEOM_FN *fn, 
  void *data
)
{
  zz->Get_Num_Geom = fn;
  zz->Get_Num_Geom_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Geom_Multi_Fn(
  ZZ *zz, 
  ZOLTAN_GEOM_MULTI_FN *fn, 
  void *data
)
{
  zz->Get_Geom_Multi = fn;
  zz->Get_Geom_Multi_Data = data;
  zz->Get_Geom = NULL;
  zz->Get_Geom_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Geom_Fn(
  ZZ *zz, 
  ZOLTAN_GEOM_FN *fn, 
  void *data
)
{
  zz->Get_Geom = fn;
  zz->Get_Geom_Data = data;
  zz->Get_Geom_Multi = NULL;
  zz->Get_Geom_Multi_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Num_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_NUM_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_Num_Obj = fn;
  zz->Get_Num_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Obj_List_Fn(
  ZZ *zz, 
  ZOLTAN_OBJ_LIST_FN *fn, 
  void *data
)
{
  zz->Get_Obj_List = fn;
  zz->Get_Obj_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_First_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_FIRST_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_First_Obj = fn;
  zz->Get_First_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Next_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_NEXT_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_Next_Obj = fn;
  zz->Get_Next_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Num_Border_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_NUM_BORDER_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_Num_Border_Obj = fn;
  zz->Get_Num_Border_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Border_Obj_List_Fn(
  ZZ *zz, 
  ZOLTAN_BORDER_OBJ_LIST_FN *fn, 
  void *data
)
{
  zz->Get_Border_Obj_List = fn;
  zz->Get_Border_Obj_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_First_Border_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_FIRST_BORDER_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_First_Border_Obj = fn;
  zz->Get_First_Border_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Next_Border_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_NEXT_BORDER_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_Next_Border_Obj = fn;
  zz->Get_Next_Border_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Obj_Size_Fn(
  ZZ *zz, 
  ZOLTAN_OBJ_SIZE_FN *fn, 
  void *data
)
{
  zz->Get_Obj_Size = fn;
  zz->Get_Obj_Size_Data = data;
  zz->Get_Obj_Size_Multi = NULL;
  zz->Get_Obj_Size_Multi_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Obj_Size_Multi_Fn(
  ZZ *zz, 
  ZOLTAN_OBJ_SIZE_MULTI_FN *fn, 
  void *data
)
{
  zz->Get_Obj_Size_Multi = fn;
  zz->Get_Obj_Size_Multi_Data = data;
  zz->Get_Obj_Size = NULL;
  zz->Get_Obj_Size_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Pack_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_PACK_OBJ_FN *fn, 
  void *data
)
{
  zz->Pack_Obj = fn;
  zz->Pack_Obj_Data = data;
  zz->Pack_Obj_Multi = NULL;
  zz->Pack_Obj_Multi_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Pack_Obj_Multi_Fn(
  ZZ *zz, 
  ZOLTAN_PACK_OBJ_MULTI_FN *fn, 
  void *data
)
{
  zz->Pack_Obj_Multi = fn;
  zz->Pack_Obj_Multi_Data = data;
  zz->Pack_Obj = NULL;
  zz->Pack_Obj_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Unpack_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_UNPACK_OBJ_FN *fn, 
  void *data
)
{
  zz->Unpack_Obj = fn;
  zz->Unpack_Obj_Data = data;
  zz->Unpack_Obj_Multi = NULL;
  zz->Unpack_Obj_Multi_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Unpack_Obj_Multi_Fn(
  ZZ *zz, 
  ZOLTAN_UNPACK_OBJ_MULTI_FN *fn, 
  void *data
)
{
  zz->Unpack_Obj_Multi = fn;
  zz->Unpack_Obj_Multi_Data = data;
  zz->Unpack_Obj = NULL;
  zz->Unpack_Obj_Data = NULL;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Num_Coarse_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_NUM_COARSE_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_Num_Coarse_Obj = fn;
  zz->Get_Num_Coarse_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Coarse_Obj_List_Fn(
  ZZ *zz, 
  ZOLTAN_COARSE_OBJ_LIST_FN *fn, 
  void *data
)
{
  zz->Get_Coarse_Obj_List = fn;
  zz->Get_Coarse_Obj_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_First_Coarse_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_FIRST_COARSE_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_First_Coarse_Obj = fn;
  zz->Get_First_Coarse_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Next_Coarse_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_NEXT_COARSE_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_Next_Coarse_Obj = fn;
  zz->Get_Next_Coarse_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Num_Child_Fn(
  ZZ *zz, 
  ZOLTAN_NUM_CHILD_FN *fn, 
  void *data
)
{
  zz->Get_Num_Child = fn;
  zz->Get_Num_Child_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Child_List_Fn(
  ZZ *zz, 
  ZOLTAN_CHILD_LIST_FN *fn, 
  void *data
)
{
  zz->Get_Child_List = fn;
  zz->Get_Child_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Child_Weight_Fn(
  ZZ *zz, 
  ZOLTAN_CHILD_WEIGHT_FN *fn, 
  void *data
)
{
  zz->Get_Child_Weight = fn;
  zz->Get_Child_Weight_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_HG_Size_CS_Fn(
  ZZ *zz, 
  ZOLTAN_HG_SIZE_CS_FN *fn, 
  void *data
)
{
  zz->Get_HG_Size_CS = fn;
  zz->Get_HG_Size_CS_Data = data;
  return ZOLTAN_OK;
}
/*****************************************************************************/

int Zoltan_Set_HG_CS_Fn(
  ZZ *zz, 
  ZOLTAN_HG_CS_FN *fn, 
  void *data
)
{
  zz->Get_HG_CS = fn;
  zz->Get_HG_CS_Data = data;
  return ZOLTAN_OK;
}
/*****************************************************************************/

int Zoltan_Set_HG_Size_Edge_Wts_Fn(
  ZZ *zz, 
  ZOLTAN_HG_SIZE_EDGE_WTS_FN *fn, 
  void *data
)
{
  zz->Get_HG_Size_Edge_Wts = fn;
  zz->Get_HG_Size_Edge_Wts_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_HG_Edge_Wts_Fn(
  ZZ *zz, 
  ZOLTAN_HG_EDGE_WTS_FN *fn, 
  void *data
)
{
  zz->Get_HG_Edge_Wts = fn;
  zz->Get_HG_Edge_Wts_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Num_Fixed_Obj_Fn(
  ZZ *zz, 
  ZOLTAN_NUM_FIXED_OBJ_FN *fn, 
  void *data
)
{
  zz->Get_Num_Fixed_Obj = fn;
  zz->Get_Num_Fixed_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Fixed_Obj_List_Fn(
  ZZ *zz, 
  ZOLTAN_FIXED_OBJ_LIST_FN *fn, 
  void *data
)
{
  zz->Get_Fixed_Obj_List = fn;
  zz->Get_Fixed_Obj_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Hier_Num_Levels_Fn(
  ZZ *zz, 
  ZOLTAN_HIER_NUM_LEVELS_FN *fn, 
  void *data
)
{
  zz->Get_Hier_Num_Levels = fn;
  zz->Get_Hier_Num_Levels_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Hier_Partition_Fn(
  ZZ *zz, 
  ZOLTAN_HIER_PARTITION_FN *fn, 
  void *data
)
{
  zz->Get_Hier_Partition = fn;
  zz->Get_Hier_Partition_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Hier_Method_Fn(
  ZZ *zz, 
  ZOLTAN_HIER_METHOD_FN *fn, 
  void *data
)
{
  zz->Get_Hier_Method = fn;
  zz->Get_Hier_Method_Data = data;
  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
