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
#include <ctype.h>
#include <string.h>
#include "params_const.h"
#include "zoltan_util.h"
#include "zz_const.h"


int       Zoltan_Check_Param(
const char *name,		/* name of parameter being reset */
const char *val,		/* new value for parameter */
PARAM_VARS * params,		/* structure describing parameters */
PARAM_UTYPE *result,		/* pointer to return value */
int *matched_index)		/* where in struct the match occurs */
{		
    char     *yo = "Zoltan_Check_Param";
    char      msg[256];
    int       i;		/* loop counter */
    int       status;		/* return code: */
    /* 0 => name found and value OK */
    /* 1 => name not found */
    /* 2 => name found, but value bad */

    status = 1;
    i = 0;
    while (params->name != NULL) {
	if (!strcmp(params->name, name)) {
	    status = 0;
	    break;
	}
	i++;
	params++;
    }

    if (status == 0) {		/* name found */
      *matched_index = i;
      if (!strcmp(val, "DEFAULT")){
         result->def = 1;
      }
      else {
	/* Figure out what type it is and read value. */
        result->def = 0;
	if (!strcmp(params->type, "INT") || !strcmp(params->type, "INTEGER")) {
	    /* First special case if True or False */
	    if (*val == 'T')
		(*result).ival = 1;
	    else if (*val == 'F')
		(*result).ival = 0;
	    else {
		/* Check that there's a digit here */
		for (i = strlen(val); i >= 0; i--)
		    if (isdigit((int)(val[i])))
			break;
		if (i < 0)
		    status = 2;
		else {
		    (*result).ival = atoi(val);
		}
	    }
	}

	else if ((!strcmp(params->type, "FLOAT")) ||
	         (!strcmp(params->type, "REAL")) ||
                 (!strcmp(params->type, "DOUBLE"))) {
	    /* Check that there's a digit here */
	    for (i = strlen(val); i >= 0; i--)
		if (isdigit((int)(val[i])))
		    break;
	    if (i < 0)
		status = 2;
	    else {
		(*result).fval = atof(val);
		(*result).dval = atof(val);
	    }
	}

	else if (!strcmp(params->type, "LONG")) {
	    /* First special case if True or False */
	    if (*val == 'T')
		(*result).lval = 1;
	    else if (*val == 'F')
		(*result).lval = 0;
	    else {
		/* Check that there's a digit here */
		for (i = strlen(val); i >= 0; i--)
		    if (isdigit((int)(val[i])))
			break;
		if (i < 0)
		    status = 2;
		else {
		    (*result).lval = atol(val);
		}
	    }
	}

	else if (!strcmp(params->type, "STRING")) {
	    strncpy((*result).sval, val, MAX_PARAM_STRING_LEN);
	}

	else if (!strcmp(params->type, "CHAR")) {
	    (*result).cval = *val;
	}

	else {
	    sprintf(msg, "Bad type for parameter `%s'", 
                    params->name);
            ZOLTAN_PRINT_WARN(-1, yo, msg);
	    status = 2;
	}
      }
    }
    else {			/* name not matched */
	*matched_index = -1;
    }

    return (status);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
