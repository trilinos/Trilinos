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

#include <stdio.h>
#include "lb_const.h"
#include "params_const.h"

/* 
 * Handle parameter changes for variables stored in LB object.
 * Currently, only example is Imbalance_Tol.
 */

int LB_Set_Key_Param(
LB *lb,                         /* load balance object */
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;			/* return code */
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */
    PARAM_VARS key_params[] = {
	{ "IMBALANCE_TOL", NULL, "DOUBLE" },
	{ "AUTO_MIGRATE", NULL, "INT" },
	{ NULL, NULL, NULL } };

    status = LB_Check_Param(name, val, key_params, &result, &index);

    if (status == 0) {
      if (index == 0) {		/* Imbalance_Tol */
	if (result.dval < 1.0) {
	    fprintf(stderr, "WARNING: Invalid Imbalance_Tol value (%g) "
		"being set to 1.0\n", result.dval);
	    result.dval = 1.0;
	}
	lb->Imbalance_Tol = result.dval;
	status = 3;		/* Don't add to Params field of LB */
      }
      else if (index == 1) {		/* Help_Migrate */
	lb->Migrate.Auto_Migrate = result.ival;
	status = 3;		/* Don't add to Params field of LB */
      }
   }

    return(status);
}
