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
#include "params_const.h"


void Zoltan_Print_Params(
  PARAM_LIST *ptr)			/* pointer to list of parameters */
{
/*
 *  Function to print out list of set parameter values.
 */

    printf("Parameter Settings\n");
    while (ptr != NULL) {
       printf("%s = %s\n",ptr->name, ptr->new_val);
       ptr = ptr->next;
    }
    printf("\n");
}
