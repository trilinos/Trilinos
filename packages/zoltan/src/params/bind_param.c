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

#include "lbi_const.h"
#include "lb_const.h"
#include "lb_util_const.h"
#include "params_const.h"

int       LB_Bind_Param(
PARAM_VARS *params,		/* parameter structure */
char *name,			/* parameter name */
void *var)			/* pointer to variable to be associated with the parameter name */
{
/*
 *  Function to bind a parameter name to a variable.
 *  On output:
 *    LB_OK indicates success.
 *    LB_WARN indicates that parameter name was not found (misspelled?).
 *            No binding took place. A warning message is printed in this case.
 *    LB_FATAL signals something more serious.
 */

    char     *name2;		/* clean version of name */
    int       flag;		/* return value from function */
    PARAM_VARS *ptr;		/* pointer to a parameter */

    /* First convert to upper case & remove leading white space. */
    flag = LB_clean_string(name, &name2);
    if (flag)
	return (flag);

    /* Search through parameter array to find name2 */

    for (ptr = params; ptr->name != NULL; ptr++) {
	if (!strcmp(name2, ptr->name)) {	/* string match */
	    ptr->ptr = var;
	    return (LB_OK);
	}
    }

    /* If we reach this point, the parameter name must be invalid */
    return (LB_WARN);
}
