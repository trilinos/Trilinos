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
#include <ctype.h>
#include <string.h>
#include "lbi_const.h"
#include "lb_const.h"
#include "params_const.h"
#include "all_allo_const.h"
#include "rcb_const.h"
#include "octupdate_const.h"
#include "parmetis_jostle_const.h"
#include "timer_const.h"

static int add_param(LB *, char *, char *);
static int clean_string(char *, char **);


int       LB_Set_Param(
LB *lb,				/* load balance object */
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

    /* Status flag is used as follows:
     *  0 -> matched and value passed all checks. 
     *  1 -> not matched. 
     *  2 -> matched, but failed checks. 
     *  3 -> matched and OK, but don't add it to params list.
     *  other -> more serious error (LB return codes) */

    /* First convert to upper case & remove leading white space. */
    flag = clean_string(name1, &name);
    if (flag)
	return (flag);
    flag = clean_string(val1, &val);
    if (flag) {
	LB_FREE(&name);
	return (flag);
    }

    /* Now call all the parameter setting routines. */
    /* New parameter routines should be added here. */

    status = LB_Set_Key_Param(lb, name, val);

    if (status == 1)
        status = LB_Set_Malloc_Param(name, val);

    if (status == 1)
        status = LB_Set_RCB_Param(name, val);

    if (status == 1)
        status = LB_Set_ParMetis_Param(name, val);

    if (status == 1)
        status = LB_Set_Octpart_Param(name, val);

    if (status == 1)
        status = LB_Set_Timer_Param(name, val);
/*
    if (status == 1)
	status = LB_Set_SFC_Param(name, val);
*/

    /* All parameter setting routines have been called, now finish up. */

    if (status == 1)		/* Parameter name never found */
	printf("Warning: parameter `%s' not found; not reset to `%s'.\n",
	       name, val);

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


static int clean_string(
char *string1,			/* original string */
char **pstring2) 		/* cleaned string to return */
{
/* Remove leading & trailing white space and convert to upper case. */

    char     *string2;		/* cleaned up string */
    int       start, end;	/* indices bounding true string */
    int       length1;		/* length of string 1 */
    int       i;		/* loop counter */

    length1 = strlen(string1);
    start = 0;
    end = length1;
    while (start < length1 && isspace(string1[start]))
	start++;
    while (end > start && isspace(string1[end]))
	end--;

    string2 = (char *)
       LB_Malloc((end - start + 1) * sizeof(char), __FILE__, __LINE__);
    *pstring2 = string2;

    if (string2 == NULL)
	return (LB_MEMERR);

    for (i = start; i < end; i++) {
	*string2++ = toupper(string1[i]);
    }
    *string2 = '\0';

    return (LB_OK);
}


static int add_param(
LB *lb,				/* load balance object */
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
    param = (LB_PARAM *) LB_Malloc(sizeof(LB_PARAM), __FILE__, __LINE__);
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
