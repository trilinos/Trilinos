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
#include "order_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines to create, access, and destroy 
 *  Zoltan Ordering Structs (ZOS).
 *  These functions are all callable by the application.  
 *
 *  NOTE: These routines cannot be used in any useful way with Zoltan 1.5!
 *  The ordering struct is currently not supported by Zoltan_Order,
 *  but this may change in future versions.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Order_Create(ZOS **order_info, ZZ *zz)
{
  int ierr = ZOLTAN_OK; /* Error code to return */
  static char *yo = "Zoltan_Order_Create";

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (zz == NULL){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Pointer to Zoltan struct is NULL");
    ierr = ZOLTAN_FATAL;
    return (ierr);
  }

  *order_info = (ZOS *) ZOLTAN_MALLOC (sizeof(ZOS));
  if (!(*order_info)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Not enough memory");
    ierr = ZOLTAN_MEMERR;
    return (ierr);
  }

  /* Initialize ordering struct */
/*   (*order_info)->zz = zz; */
  (*order_info)->nbr_objects = 0;
  (*order_info)->gids = NULL;
  (*order_info)->lids = NULL;
  (*order_info)->method[0] = '\0';
  (*order_info)->num_separators = 0;
  (*order_info)->sep_sizes = NULL;
  (*order_info)->rank = NULL;
  (*order_info)->vtxdist = NULL;

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

int Zoltan_Order_Destroy(ZOS **order_info)
{
  int ierr = ZOLTAN_OK; /* Error code to return */
  /* static char *yo = "Zoltan_Order_Destroy"; */

  /* ZOLTAN_TRACE_ENTER(zz, yo); */

  if ((*order_info)->gids)      ZOLTAN_FREE(&((*order_info)->gids));
  if ((*order_info)->lids)      ZOLTAN_FREE(&((*order_info)->lids));
  if ((*order_info)->sep_sizes) ZOLTAN_FREE(&((*order_info)->sep_sizes));
  if ((*order_info)->rank)      ZOLTAN_FREE(&((*order_info)->rank));
  if ((*order_info)->vtxdist)   ZOLTAN_FREE(&((*order_info)->vtxdist));

  ZOLTAN_FREE(order_info);
  order_info = NULL;

  /* ZOLTAN_TRACE_EXIT(zz, yo); */
  return (ierr);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
