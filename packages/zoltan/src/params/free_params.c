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


#include <stdio.h>
#include <stdlib.h>
#include "params_const.h"
#include "zoltan_mem.h"
#include "zoltan_types.h"
#include "zz_util_const.h"
#include "zz_const.h"


void Zoltan_Free_Params(
PARAM_LIST **params)				/* parameters structure */
{
/*
 * Free the list of new parameter values.
 */
    PARAM_LIST *ptr, *next;

    if (params == NULL) return;

    ptr = *params;
    while (ptr != NULL) {
	next = ptr->next;
	ZOLTAN_FREE(&(ptr->name));
	ZOLTAN_FREE(&(ptr->new_val));
	ZOLTAN_FREE(&ptr);
	ptr = next;
    }

    *params = NULL;
}

int Zoltan_Copy_Params(PARAM_LIST **to, PARAM_LIST const *from)
{
  PARAM_LIST *param;
  PARAM_LIST *prev;

  if (*to != NULL) {
    Zoltan_Free_Params(to);
  }

  prev = NULL;

  while (from) {
    
    param = (PARAM_LIST *) ZOLTAN_MALLOC(sizeof(PARAM_LIST));
    if (param == NULL) {
      Zoltan_Free_Params(to);
      return ZOLTAN_MEMERR;
    }

    param->name = Zoltan_Strdup(from->name);
    param->new_val = Zoltan_Strdup(from->new_val);
    param->index = from->index;
    param->next = NULL;

    if (prev){
      prev->next = param;
    }

    from = from->next;
    prev = param;

    if (*to == NULL){
      *to = param;
    }
  }

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
