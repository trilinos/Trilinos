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
int Zoltan_Set_Pre_Migrate_PP_Fn(
  ZZ *zz,
  ZOLTAN_PRE_MIGRATE_PP_FN *fn,
  void *data
)
{
  zz->Migrate.Pre_Migrate_PP = fn;
  zz->Migrate.Pre_Migrate_PP_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Mid_Migrate_PP_Fn(
  ZZ *zz,
  ZOLTAN_MID_MIGRATE_PP_FN *fn,
  void *data
)
{
  zz->Migrate.Mid_Migrate_PP = fn;
  zz->Migrate.Mid_Migrate_PP_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Post_Migrate_PP_Fn(
  ZZ *zz,
  ZOLTAN_POST_MIGRATE_PP_FN *fn,
  void *data
)
{
  zz->Migrate.Post_Migrate_PP = fn;
  zz->Migrate.Post_Migrate_PP_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Pre_Migrate_Fn(
  ZZ *zz,
  ZOLTAN_PRE_MIGRATE_FN *fn,
  void *data
)
{
  zz->Migrate.Pre_Migrate = fn;
  zz->Migrate.Pre_Migrate_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Mid_Migrate_Fn(
  ZZ *zz,
  ZOLTAN_MID_MIGRATE_FN *fn,
  void *data
)
{
  zz->Migrate.Mid_Migrate = fn;
  zz->Migrate.Mid_Migrate_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/

int Zoltan_Set_Post_Migrate_Fn(
  ZZ *zz,
  ZOLTAN_POST_MIGRATE_FN *fn,
  void *data
)
{
  zz->Migrate.Post_Migrate = fn;
  zz->Migrate.Post_Migrate_Data = data;
  return ZOLTAN_OK;
}

/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
