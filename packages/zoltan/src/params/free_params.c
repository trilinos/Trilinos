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

#include <stdio.h>
#include <stdlib.h>
#include "params_const.h"
#include "zoltan_mem.h"


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
