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
#ifdef __STDC__
#include <string.h>
#else
#include <strings.h>
#endif  /* __STDC__ */
#include "lbi_const.h"
#include "lb_const.h"
#include "lb_util_const.h"
#include "params_const.h"
#include "all_allo_const.h"
#include "rcb_const.h"
#include "octupdate_const.h"
#include "parmetis_jostle_const.h"
#include "reftree_const.h"
#include "timer_const.h"
#include "ha_const.h"

static int add_param(LB *, char *, char *);

/* List of set_parameter functions to be called */
static LB_SET_PARAM_FN * Param_func[] = {
        LB_Set_Malloc_Param,
        LB_Set_RCB_Param,
        LB_Set_ParMetis_Param,
        LB_Set_Jostle_Param,
        LB_Set_Octpart_Param,
        LB_Set_Reftree_Param,
        LB_Set_Timer_Param,
     /* LB_Set_Machine_Param, */
        /*** Add your new parameter setting function here! ***/
        NULL /* Last entry must be NULL! */
};

int       LB_Set_Param(
LB *lb,				/* load balance structure */
char *name1,			/* parameter name */
char *val1)			/* value to set this parameter to */
{
/*
 *  Function to set a parameter value.
 *  On output:
 *    LB_OK indicates success.
 *    LB_WARN indicates that parameter was not set properly (misspelled?).
 *             A warning message is printed in this case.
 *    LB_FATAL signals something more serious.
 */

    char     *name, *val;	/* clean versions of name1, val1 */
    int       flag;		/* return value from function */
    int       status;		/* has character string been matched? */
    LB_SET_PARAM_FN **func;     /* pointer to parameter setting functions */

    /* Status flag is used as follows:
     *  0 -> matched and value passed all checks. 
     *  1 -> not matched. 
     *  2 -> matched, but failed checks. 
     *  3 -> matched and OK, but don't add it to params list.
     *  other -> more serious error (LB return codes) */

    /* First convert to upper case & remove leading white space. */
    flag = LB_clean_string(name1, &name);
    if (flag)
	return (flag);
    flag = LB_clean_string(val1, &val);
    if (flag) {
	LB_FREE(&name);
	return (flag);
    }

    /* Call the key parameter routine. This one is Zoltan-specific. */
    status = LB_Set_Key_Param(lb, name, val);

    /* Now call all the other parameter setting routines. */
    for (func = Param_func; *func != NULL; func++) {
        if (status == 1)
            status = (**func)(name, val);
    }

    /* All parameter setting routines have been called, now finish up. */

    if (status == 1)		/* Parameter name never found */
	fprintf(stderr, "Warning: parameter `%s' not found;"
                        " not reset to `%s'.\n", name, val);

    if (status == 0)		/* Parameter OK, add it to list */
	add_param(lb, name, val);
    else {
	LB_FREE(&name);
	LB_FREE(&val);
    }

    if (status == 0)
	flag = LB_OK;
    else if (status == 1 || status == 2)
	flag = LB_WARN;
    else
	flag = status;
    return (flag);
}


static int add_param(
LB *lb,				/* load balance structure */
char *name,			/* parameter name */
char *val)			/* value to set this parameter to */
{
/*
 * Parameter checked out OK.  Add it to linked list of param values.
 * Search through existing list to replace value if its there.
 * Otherwise, add it to the end of the list.
 */
    LB_PARAM *ptr, **pptr;	/* loops through parameter list */
    LB_PARAM *param;		/* parameter entry in list */


    pptr = &(lb->Params);
    while (*pptr != NULL) {
	ptr = *pptr;
	if (!strcmp(name, ptr->name)) {	/* string match */
	    LB_FREE(&(ptr->new_val));
	    ptr->new_val = val;
	    return (LB_OK);
	}
	pptr = &(ptr->next);
    }

    /* This is a new parameter, add it to list. */
    param = (LB_PARAM *) LB_MALLOC(sizeof(LB_PARAM));
    if (param == NULL) {
	LB_FREE(&name);
	LB_FREE(&val);
	return (LB_MEMERR);
    }
    *pptr = param;
    param->next = NULL;
    param->name = name;
    param->new_val = val;

    return (LB_OK);
}
