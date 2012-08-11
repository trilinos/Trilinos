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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zz_util_const.h"
#include "params_const.h"
#include "zoltan_mem.h"
#include "zoltan_util.h"
#include "zz_const.h"

int       Zoltan_Bind_Param(
PARAM_VARS *params,		/* parameter structure */
char *name,			/* parameter name */
void *var 			/* pointer to variable to be associated 
				   with the parameter name */
)
{
  return Zoltan_Bind_Param_Vec(params, name, var, 0);
}

int       Zoltan_Bind_Param_Vec(
PARAM_VARS *params,		/* parameter structure */
char *name,			/* parameter name */
void *var,			/* pointer to variable to be associated 
				   with the parameter name */
int dim				/* dimension of parameter vector */
)
{
/*
 *  Function to bind a parameter name to a variable.
 *  On output:
 *    ZOLTAN_OK indicates success.
 *    ZOLTAN_WARN indicates that parameter name was not found (misspelled?).
 *            No binding took place. A warning message is printed in this case.
 *    ZOLTAN_FATAL signals something more serious.
 */

    char     *yo = "Zoltan_Bind_Param";
    char      msg[256];
    char     *name2;		/* clean version of name */
    int       flag;		/* return value from function */
    PARAM_VARS *ptr;		/* pointer to a parameter */

    /* First convert to upper case & remove leading white space. */
    flag = Zoltan_Clean_String(name, &name2);
    if (flag)
	return (flag);

    /* Search through parameter array to find name2 */

    for (ptr = params; ptr->name != NULL; ptr++) {
	if (!strcmp(name2, ptr->name)) {	/* string match */
	    ptr->ptr = var;
	    ptr->length = dim;
            ZOLTAN_FREE(&name2);
	    return (ZOLTAN_OK);
	}
    }

    /* If we reach this point, the parameter name must be invalid */
    sprintf(msg, "Parameter name %s not found; it will"
                 "not be bound to any variable.", name2);
    ZOLTAN_PRINT_WARN(-1, yo, msg);
    ZOLTAN_FREE(&name2);
    return (ZOLTAN_WARN);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
