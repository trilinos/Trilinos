/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines used to register callback functions.
 *  A generic routine (LB_Set_Fn) is provided, as well as callback-function
 *  specific registration routines that enable tighter type-checking by
 *  the compiler.
 *
 *  When new callback functions are added to Zoltan, they should be 
 *  added to the case statement in LB_Set_Fn and should have separate
 *  registration functions specific to the new callback types.
 */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Set_Fn(LB *lb, LB_FN_TYPE fn_type, void (*fn)(), void *data)
{
/*
 *  Function to initialize a given LB interface function.
 *  Input:
 *    lb                --  Pointer to a LB structure.
 *    fn_type           --  Enum type indicating the function to be set.
 *    fn                --  Pointer to the function to be used in the
 *                          assignment.
 *    data              --  Pointer to data that the LB library will
 *                          pass as an argument to fn(). May be NULL.
 *  Output:
 *    lb                --  Appropriate field set to value in void *().
 */

char *yo = "LB_Set_Fn";
char msg[256];
int ierr;

  switch (fn_type) {
  case LB_NUM_EDGES_FN_TYPE:
    ierr = LB_Set_Num_Edges_Fn(lb, (LB_NUM_EDGES_FN *) fn, data);
    break;
  case LB_EDGE_LIST_FN_TYPE:
    ierr = LB_Set_Edge_List_Fn(lb, (LB_EDGE_LIST_FN *) fn, data);
    break;
  case LB_NUM_GEOM_FN_TYPE:
    ierr = LB_Set_Num_Geom_Fn(lb, (LB_NUM_GEOM_FN *) fn, data);
    break;
  case LB_GEOM_FN_TYPE:
    ierr = LB_Set_Geom_Fn(lb, (LB_GEOM_FN *) fn, data);
    break;
  case LB_NUM_OBJ_FN_TYPE:
    ierr = LB_Set_Num_Obj_Fn(lb, (LB_NUM_OBJ_FN *) fn, data);
    break;
  case LB_OBJ_LIST_FN_TYPE:
    ierr = LB_Set_Obj_List_Fn(lb, (LB_OBJ_LIST_FN *) fn, data);
    break;
  case LB_FIRST_OBJ_FN_TYPE:
    ierr = LB_Set_First_Obj_Fn(lb, (LB_FIRST_OBJ_FN *) fn, data);
    break;
  case LB_NEXT_OBJ_FN_TYPE:
    ierr = LB_Set_Next_Obj_Fn(lb, (LB_NEXT_OBJ_FN *) fn, data);
    break;
  case LB_NUM_BORDER_OBJ_FN_TYPE:
    ierr = LB_Set_Num_Border_Obj_Fn(lb, (LB_NUM_BORDER_OBJ_FN *) fn, data);
    break;
  case LB_BORDER_OBJ_LIST_FN_TYPE:
    ierr = LB_Set_Border_Obj_List_Fn(lb, (LB_BORDER_OBJ_LIST_FN *) fn, data);
    break;
  case LB_FIRST_BORDER_OBJ_FN_TYPE:
    ierr = LB_Set_First_Border_Obj_Fn(lb, (LB_FIRST_BORDER_OBJ_FN *) fn, data);
    break;
  case LB_NEXT_BORDER_OBJ_FN_TYPE:
    ierr = LB_Set_Next_Border_Obj_Fn(lb, (LB_NEXT_BORDER_OBJ_FN *) fn, data);
    break;
  case LB_PRE_MIGRATE_FN_TYPE:
    ierr = LB_Set_Pre_Migrate_Fn(lb, (LB_PRE_MIGRATE_FN *) fn, data);
    break;
  case LB_MID_MIGRATE_FN_TYPE:
    ierr = LB_Set_Mid_Migrate_Fn(lb, (LB_MID_MIGRATE_FN *) fn, data);
    break;
  case LB_POST_MIGRATE_FN_TYPE:
    ierr = LB_Set_Post_Migrate_Fn(lb, (LB_POST_MIGRATE_FN *) fn, data);
    break;
  case LB_OBJ_SIZE_FN_TYPE:
    ierr = LB_Set_Obj_Size_Fn(lb, (LB_OBJ_SIZE_FN *) fn, data);
    break;
  case LB_PACK_OBJ_FN_TYPE:
    ierr = LB_Set_Pack_Obj_Fn(lb, (LB_PACK_OBJ_FN *) fn, data);
    break;
  case LB_UNPACK_OBJ_FN_TYPE:
    ierr = LB_Set_Unpack_Obj_Fn(lb, (LB_UNPACK_OBJ_FN *) fn, data);
    break;
  case LB_NUM_COARSE_OBJ_FN_TYPE:
    ierr = LB_Set_Num_Coarse_Obj_Fn(lb, (LB_NUM_COARSE_OBJ_FN *) fn, data);
    break;
  case LB_COARSE_OBJ_LIST_FN_TYPE:
    ierr = LB_Set_Coarse_Obj_List_Fn(lb, (LB_COARSE_OBJ_LIST_FN *) fn, data);
    break;
  case LB_FIRST_COARSE_OBJ_FN_TYPE:
    ierr = LB_Set_First_Coarse_Obj_Fn(lb, (LB_FIRST_COARSE_OBJ_FN *) fn, data);
    break;
  case LB_NEXT_COARSE_OBJ_FN_TYPE:
    ierr = LB_Set_Next_Coarse_Obj_Fn(lb, (LB_NEXT_COARSE_OBJ_FN *) fn, data);
    break;
  case LB_NUM_CHILD_FN_TYPE:
    ierr = LB_Set_Num_Child_Fn(lb, (LB_NUM_CHILD_FN *) fn, data);
    break;
  case LB_CHILD_LIST_FN_TYPE:
    ierr = LB_Set_Child_List_Fn(lb, (LB_CHILD_LIST_FN *) fn, data);
    break;
  case LB_CHILD_WEIGHT_FN_TYPE:
    ierr = LB_Set_Child_Weight_Fn(lb, (LB_CHILD_WEIGHT_FN *) fn, data);
    break;
  default:
    sprintf(msg, "LB_FN_TYPE %d is invalid.\n"
                 "Value must be in range 0 to %d.", fn_type, LB_MAX_FN_TYPES);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
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

int LB_Set_Num_Edges_Fn(LB *lb, LB_NUM_EDGES_FN *fn, void *data)
{
  lb->Get_Num_Edges = fn;
  lb->Get_Num_Edges_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Edge_List_Fn(LB *lb, LB_EDGE_LIST_FN *fn, void *data)
{
  lb->Get_Edge_List = fn;
  lb->Get_Edge_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Num_Geom_Fn(LB *lb, LB_NUM_GEOM_FN *fn, void *data)
{
  lb->Get_Num_Geom = fn;
  lb->Get_Num_Geom_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Geom_Fn(LB *lb, LB_GEOM_FN *fn, void *data)
{
  lb->Get_Geom = fn;
  lb->Get_Geom_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Num_Obj_Fn(LB *lb, LB_NUM_OBJ_FN *fn, void *data)
{
  lb->Get_Num_Obj = fn;
  lb->Get_Num_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Obj_List_Fn(LB *lb, LB_OBJ_LIST_FN *fn, void *data)
{
  lb->Get_Obj_List = fn;
  lb->Get_Obj_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_First_Obj_Fn(LB *lb, LB_FIRST_OBJ_FN *fn, void *data)
{
  lb->Get_First_Obj = fn;
  lb->Get_First_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Next_Obj_Fn(LB *lb, LB_NEXT_OBJ_FN *fn, void *data)
{
  lb->Get_Next_Obj = fn;
  lb->Get_Next_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Num_Border_Obj_Fn(LB *lb, LB_NUM_BORDER_OBJ_FN *fn, void *data)
{
  lb->Get_Num_Border_Obj = fn;
  lb->Get_Num_Border_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Border_Obj_List_Fn(LB *lb, LB_BORDER_OBJ_LIST_FN *fn, void *data)
{
  lb->Get_Border_Obj_List = fn;
  lb->Get_Border_Obj_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_First_Border_Obj_Fn(LB *lb, LB_FIRST_BORDER_OBJ_FN *fn, void *data)
{
  lb->Get_First_Border_Obj = fn;
  lb->Get_First_Border_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Next_Border_Obj_Fn(LB *lb, LB_NEXT_BORDER_OBJ_FN *fn, void *data)
{
  lb->Get_Next_Border_Obj = fn;
  lb->Get_Next_Border_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Pre_Migrate_Fn(LB *lb, LB_PRE_MIGRATE_FN *fn, void *data)
{
  lb->Migrate.Pre_Migrate = fn;
  lb->Migrate.Pre_Migrate_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Mid_Migrate_Fn(LB *lb, LB_MID_MIGRATE_FN *fn, void *data)
{
  lb->Migrate.Mid_Migrate = fn;
  lb->Migrate.Mid_Migrate_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Post_Migrate_Fn(LB *lb, LB_POST_MIGRATE_FN *fn, void *data)
{
  lb->Migrate.Post_Migrate = fn;
  lb->Migrate.Post_Migrate_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Obj_Size_Fn(LB *lb, LB_OBJ_SIZE_FN *fn, void *data)
{
  lb->Migrate.Get_Obj_Size = fn;
  lb->Migrate.Get_Obj_Size_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Pack_Obj_Fn(LB *lb, LB_PACK_OBJ_FN *fn, void *data)
{
  lb->Migrate.Pack_Obj = fn;
  lb->Migrate.Pack_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Unpack_Obj_Fn(LB *lb, LB_UNPACK_OBJ_FN *fn, void *data)
{
  lb->Migrate.Unpack_Obj = fn;
  lb->Migrate.Unpack_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Num_Coarse_Obj_Fn(LB *lb, LB_NUM_COARSE_OBJ_FN *fn, void *data)
{
  lb->Get_Num_Coarse_Obj = fn;
  lb->Get_Num_Coarse_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Coarse_Obj_List_Fn(LB *lb, LB_COARSE_OBJ_LIST_FN *fn, void *data)
{
  lb->Get_Coarse_Obj_List = fn;
  lb->Get_Coarse_Obj_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_First_Coarse_Obj_Fn(LB *lb, LB_FIRST_COARSE_OBJ_FN *fn, void *data)
{
  lb->Get_First_Coarse_Obj = fn;
  lb->Get_First_Coarse_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Next_Coarse_Obj_Fn(LB *lb, LB_NEXT_COARSE_OBJ_FN *fn, void *data)
{
  lb->Get_Next_Coarse_Obj = fn;
  lb->Get_Next_Coarse_Obj_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Num_Child_Fn(LB *lb, LB_NUM_CHILD_FN *fn, void *data)
{
  lb->Get_Num_Child = fn;
  lb->Get_Num_Child_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Child_List_Fn(LB *lb, LB_CHILD_LIST_FN *fn, void *data)
{
  lb->Get_Child_List = fn;
  lb->Get_Child_List_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int LB_Set_Child_Weight_Fn(LB *lb, LB_CHILD_WEIGHT_FN *fn, void *data)
{
  lb->Get_Child_Weight = fn;
  lb->Get_Child_Weight_Data = data;
  return ZOLTAN_OK;
}
