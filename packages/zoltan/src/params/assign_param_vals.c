/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_assignparamvalsc_id = "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "lb_const.h"
#include "lbi_const.h"
#include "params_const.h"


void      LB_Assign_Param_Vals(
LB_PARAM * change_list,		/* list of parameter values being changed */
PARAM_VARS * params)		/* structure describing parameters */
{	
    char     *name;		/* name of parameter being reset */
    char     *val;		/* new value for parameter */
    int       found;		/* is name found? */
    PARAM_VARS *param_ptr;       /* pointer to current param */

    while (change_list != NULL) {
        param_ptr = params;
	name = change_list->name;
	val = change_list->new_val;

	found = 0;
	while (param_ptr->name != NULL) {
	    if (!strcmp(param_ptr->name, name)) {
		found = 1;
		break;
	    }
	    param_ptr++;
	}

	if (found) {		/* name found */
	    /* Figure out what type it is and read value. */
	    if (!strcmp(param_ptr->type, "INT") || 
                !strcmp(param_ptr->type, "INTEGER")) {
		/* First special case if True or False */
		if (*val == 'T')
		    *((int *) param_ptr->ptr) = 1;
		else if (*val == 'F')
		    *((int *) param_ptr->ptr) = 0;
		else {
		    *((int *) param_ptr->ptr) = atoi(val);
		}
	    }

	    else if (!strcmp(param_ptr->type, "DOUBLE")) {
		*((double *) param_ptr->ptr) = atof(val);
	    }

	    else if (!strcmp(param_ptr->type, "LONG")) {
		/* First special case if True or False */
		if (*val == 'T')
		    *((long *) param_ptr->ptr) = 1;
		else if (*val == 'F')
		    *((long *) param_ptr->ptr) = 0;
		else {
		    *((long *) param_ptr->ptr) = atol(val);
		}
	    }

	    else if (!strcmp(param_ptr->type, "STRING")) {
		strncpy((char *) param_ptr->ptr, val, MAX_PARAM_STRING_LEN);
	    }

	    else if (!strcmp(param_ptr->type, "CHAR")) {
		*((char *) param_ptr->ptr) = *val;
	    }
	}

	change_list = change_list->next;
    }
}
