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
    int       i;		/* loop counter */

    while (change_list != NULL) {
	name = change_list->name;
	val = change_list->new_val;

	found = 0;
	while (params->name != NULL) {
	    if (!strcmp(params->name, name)) {
		found = 1;
		break;
	    }
	    params++;
	}

	if (found) {		/* name found */
	    /* Figure out what type it is and read value. */
	    if (!strcmp(params->type, "INT") || !strcmp(params->type, "INTEGER")) {
		/* First special case if True or False */
		if (*val == 'T')
		    *((int *) params->ptr) = 1;
		else if (*val == 'F')
		    *((int *) params->ptr) = 0;
		else {
		    *((int *) params->ptr) = atoi(val);
		}
	    }

	    else if (!strcmp(params->type, "DOUBLE")) {
		*((double *) params->ptr) = atof(val);
	    }

	    else if (!strcmp(params->type, "LONG")) {
		/* First special case if True or False */
		if (*val == 'T')
		    *((long *) params->ptr) = 1;
		else if (*val == 'F')
		    *((long *) params->ptr) = 0;
		else {
		    *((long *) params->ptr) = atol(val);
		}
	    }

	    else if (!strcmp(params->type, "STRING")) {
		strncpy((char *) params->ptr, val, MAX_PARAM_STRING_LEN);
	    }

	    else if (!strcmp(params->type, "CHAR")) {
		*((char *) params->ptr) = *val;
	    }
	}

	change_list = change_list->next;
    }
}
