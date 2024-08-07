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
 
size_t Zoltan_Serialize_Params_Size(struct Zoltan_Struct const *from) 
{
  /* Count the number of parameters */
  PARAM_LIST const *param = from->Params;
  int nParam = 0;
  while (param) {
    nParam++;
    param = param->next;
  }

  return sizeof(int)                          /* to store number of params */
       + nParam * (MAX_PARAM_STRING_LEN * 2); /* max per param: (name, value) */
}

int Zoltan_Serialize_Params(struct Zoltan_Struct const *from, char **buf)
{
  /* Serialize the parameters */
  char *bufptr = *buf;
  PARAM_LIST const *param = from->Params;
  size_t paramSize = Zoltan_Serialize_Params_Size(from);

  /* Pack number of parameters */
  int nParam = paramSize / (MAX_PARAM_STRING_LEN * 2);
  *((int *) bufptr) = nParam;
  bufptr += sizeof(int);

  /* Pack each parameter, using max string length bytes per string */
  while (param) {
    strcpy(bufptr, param->name);
    bufptr += MAX_PARAM_STRING_LEN;
    strcpy(bufptr, param->new_val);
    bufptr += MAX_PARAM_STRING_LEN;
    param = param->next;
  }
  *buf = bufptr;

  return ZOLTAN_OK;
}

int Zoltan_Deserialize_Params(struct Zoltan_Struct *to, char **buf)
{
  /* Serialize the parameters */
  char *bufptr = *buf;
  int i;

  /* Unpack number of parameters */
  int nParam = *((int *)bufptr);
  bufptr += sizeof(int);

  /* Unpack parameters' (name, value) pairs and set them */
  for (i = 0; i < nParam; i++) {
    char *pname = bufptr;
    char *pval = bufptr + MAX_PARAM_STRING_LEN;
    Zoltan_Set_Param(to, pname, pval);
    bufptr += 2 * MAX_PARAM_STRING_LEN;
  }
  *buf = bufptr;

  return ZOLTAN_OK;
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
