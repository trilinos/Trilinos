// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


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
