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
#include <string.h>
#include "params_const.h"
#include "zoltan_types.h"
#include "zoltan_util.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Print_Assigned_Param_Vals(PARAM_VARS * );

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int      Zoltan_Assign_Param_Vals(
PARAM_LIST * change_list,		/* list of parameter values being changed */
PARAM_VARS * params,		/* structure describing parameters        */
int debug_level,                /* level for output of debugging info     */
int proc,                       /* processor # (controls debug printing)  */
int print_proc                  /* processor that should perform printing */
)
{	
    char     *yo = "Zoltan_Assign_Param_Vals";
    char      msg[256];
    char     *name;		/* name of parameter being reset */
    char     *val;		/* new value for parameter       */
    int       found;		/* is name found?                */
    int       ierr;		/* error code                    */
    PARAM_VARS *param_ptr;      /* pointer to current param      */

    ierr = ZOLTAN_OK;

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
             ierr = ZOLTAN_WARN;
             if (debug_level > 0 && proc == print_proc) {
                sprintf(msg, "Parameter %s is not bound "
                       "to any variable.  Parameter ignored.\n", 
                        param_ptr->name);
                ZOLTAN_PRINT_WARN(proc, yo, msg);
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
	    }

	    else if ((!strcmp(param_ptr->type, "FLOAT")) ||
	             (!strcmp(param_ptr->type, "REAL"))) {
		*((float *) param_ptr->ptr) = atof(val);
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
      }

      change_list = change_list->next;
    }

    if (debug_level > 0 && proc == print_proc)
        Zoltan_Print_Assigned_Param_Vals(params);

    return ierr;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Print_Assigned_Param_Vals(
PARAM_VARS * params 		/* structure describing parameters */
)
{
/* Prints the parameter values in PARAM_VARS *param.     */
PARAM_VARS *param_ptr;       /* pointer to current param */
param_ptr = params;

    while (param_ptr->name != NULL) {
      if (param_ptr->ptr != NULL) {
        if (!strcmp(param_ptr->type, "INT") || 
            !strcmp(param_ptr->type, "INTEGER")) {
 
            printf("ZOLTAN Parameter %s = %d\n", 
                    param_ptr->name, *((int *) param_ptr->ptr));
                
        }
        else if (!strcmp(param_ptr->type, "DOUBLE")) {
            printf("ZOLTAN Parameter %s = %f\n", 
                    param_ptr->name, *((double *) param_ptr->ptr));
        }
        else if (!strcmp(param_ptr->type, "LONG")) {
            printf("ZOLTAN Parameter %s = %ld\n", 
                    param_ptr->name, *((long *) param_ptr->ptr));
        }
        else if (!strcmp(param_ptr->type, "STRING")) {
            printf("ZOLTAN Parameter %s = %s\n", 
                    param_ptr->name, (char *) param_ptr->ptr);
        }
        else if (!strcmp(param_ptr->type, "CHAR")) {
            printf("ZOLTAN Parameter %s = %c\n", 
                    param_ptr->name, *((char *) param_ptr->ptr));
        }
      }
      param_ptr++;
    }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
