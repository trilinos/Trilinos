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

#include "lb_const.h"
#include "lbi_const.h"
#include "params_const.h"


int      LB_Assign_Param_Vals(
LB_PARAM * change_list,		/* list of parameter values being changed */
PARAM_VARS * params,		/* structure describing parameters */
int debug_level,                /* level for output of debugging info */
int proc                        /* processor # (controls debug printing) */
)
{	
    char     *name;		/* name of parameter being reset */
    char     *val;		/* new value for parameter */
    int       found;		/* is name found? */
    int       ierr;		/* error code */
    PARAM_VARS *param_ptr;       /* pointer to current param */

    ierr = LB_OK;

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

          /* Check that param_ptr->ptr isn't NULL */
          if (param_ptr->ptr == NULL) {
             ierr = LB_WARN;
             if (debug_level > 0 && proc == 0) {
                fprintf(stderr, "Zoltan Warning: Parameter %s is not bound to any variable. "
                       "Parameter ignored.\n", param_ptr->name);
             }
          }
          else {

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
 
                if (debug_level > 0 && proc == 0) {
                    printf("Assign Parameter %s = %d\n", 
                            param_ptr->name, *((int *) param_ptr->ptr));
                }
	    }

	    else if (!strcmp(param_ptr->type, "DOUBLE")) {
		*((double *) param_ptr->ptr) = atof(val);
 
                if (debug_level > 0 && proc == 0) {
                    printf("Assign Parameter %s = %lf\n", 
                            param_ptr->name, *((double *) param_ptr->ptr));
                }
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
 
                if (debug_level > 0 && proc == 0) {
                    printf("Assign Parameter %s = %ld\n", 
                            param_ptr->name, *((long *) param_ptr->ptr));
                }
	    }

	    else if (!strcmp(param_ptr->type, "STRING")) {
		strncpy((char *) param_ptr->ptr, val, MAX_PARAM_STRING_LEN);
 
                if (debug_level > 0 && proc == 0) {
                    printf("Assign Parameter %s = %s\n", 
                            param_ptr->name, (char *) param_ptr->ptr);
                }
	    }

	    else if (!strcmp(param_ptr->type, "CHAR")) {
		*((char *) param_ptr->ptr) = *val;
 
                if (debug_level > 0 && proc == 0) {
                    printf("Assign Parameter %s = %c\n", 
                            param_ptr->name, *((char *) param_ptr->ptr));
                }
	    }
	}
      }

      change_list = change_list->next;
    }

    return ierr;
}
