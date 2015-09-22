/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


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
