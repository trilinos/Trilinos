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
#include <string.h>
#include "key_params.h"
#include "zz_util_const.h"
#include "params_const.h"
#include "rcb_const.h"
#ifdef ZOLTAN_OCT
#include "octupdate_const.h"
#endif
#include "parmetis_jostle_const.h"
#include "reftree_const.h"
#include "timer_const.h"
#include "ha_const.h"
#include "rib_const.h"
#include "hsfc_const.h"
#include "all_allo_const.h"
#include "order_const.h"
#ifdef ZOLTAN_HG
#include "hg_const.h"
#include "phg_const.h"
#endif

static int add_param(ZZ *, char *, char *, int);
static int remove_param(ZZ *, char *, int);

/* List of set_parameter functions to be called */
static ZOLTAN_SET_PARAM_FN * Param_func[] = {
       Zoltan_Set_Malloc_Param,
       Zoltan_RCB_Set_Param,
       Zoltan_ParMetis_Set_Param,
       Zoltan_Jostle_Set_Param,
#ifdef ZOLTAN_OCT
       Zoltan_Oct_Set_Param,
#endif
       Zoltan_Reftree_Set_Param,
       Zoltan_RIB_Set_Param,
       Zoltan_HSFC_Set_Param,
       Zoltan_Order_Set_Param,
#ifdef ZOLTAN_HG
       Zoltan_HG_Set_Param,
       Zoltan_PHG_Set_Param,
#endif
       /* Zoltan_Set_Machine_Param, */
       /*** Add your new parameter setting function here! ***/
       NULL /* Last entry _must_ be NULL! */
};

int Zoltan_Set_Param(
ZZ *zz,				/* Zoltan structure */
char *name1,			/* parameter name */
char *val1)			/* value to set this parameter to */
{
    return Zoltan_Set_Param_Vec(zz, name1, val1, -1);
}

int Zoltan_Set_Param_Vec(
ZZ *zz,				/* Zoltan structure */
char *name1,			/* parameter name */
char *val1,			/* value to set this parameter to */
int index			/* index of vector parameter; -1 if scalar */
)
{
/*
 *  Function to set a parameter value.
 *  On output:
 *    ZOLTAN_OK indicates success.
 *    ZOLTAN_WARN indicates that parameter was not set properly (misspelled?).
 *             A warning message is printed in this case.
 *    ZOLTAN_FATAL signals something more serious.
 */

    char     *yo = "Zoltan_Set_Param_Vec";
    char      msg[256];
    char     *name, *val;	/* clean versions of name1, val1 */
    int       flag;		/* return value from function */
    int       status;		/* has character string been matched? */
    ZOLTAN_SET_PARAM_FN **func; /* pointer to parameter setting functions */

    /* Status flag is used as follows:
     *  0 -> matched and value passed all checks. 
     *  1 -> not matched. 
     *  2 -> matched, but failed checks. 
     *  3 -> matched and OK, but don't add it to params list.
     *  other -> more serious error (Zoltan return codes) */

    /* First convert to upper case & remove leading white space. */
    flag = Zoltan_Clean_String(name1, &name);
    if (flag)
	return (flag);
    flag = Zoltan_Clean_String(val1, &val);
    if (flag) {
	ZOLTAN_FREE(&name);
	return (flag);
    }

    /* Call the key parameter routine. This one is Zoltan-specific. */
    /* Currently, only Set_Key_Param allows vector parameters. */
    status = Zoltan_Set_Key_Param(zz, name, val, index);

    /* Now call all the other parameter setting routines. */
    /* Note that we don't pass on the index since vector parameters
       are currently only supported for key parameters. */
    for (func = Param_func; (status == 1) && (*func != NULL); func++) {
        status = (**func)(name, val);
    }

    /* All parameter setting routines have been called, now finish up. */

    if (status == 1) {		/* Parameter name never found */
	sprintf(msg, "Parameter `%s' not found; not reset to `%s'.\n", 
                name, val);
        ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
	ZOLTAN_FREE(&name);
    	ZOLTAN_FREE(&val);
    }
    else {
        if (!strcmp(val, "DEFAULT")){
	    remove_param(zz, name, index); /* Remove parameter from list */
            status = 0; 		/* "DEFAULT" is always valid */
    	    ZOLTAN_FREE(&name);
    	    ZOLTAN_FREE(&val);
        }
        else if (status == 0){		/* Parameter OK */
    	    add_param(zz, name, val, index); 	/* Add parameter to list */
        }
        else { 				/* Parameter not OK. Don't add.  */
    	    ZOLTAN_FREE(&name);             /* (It may be used to set values */
    	    ZOLTAN_FREE(&val);              /* directly in zz rather than in */
                                        /* the parameter list.)          */
        }
    }

    if (status == 0 || status == 3)
	flag = ZOLTAN_OK;
    else if (status == 1 || status == 2)
	flag = ZOLTAN_WARN;
    else
	flag = status;
    return (flag);
}


static int add_param(
ZZ *zz,				/* Zoltan structure */
char *name,			/* parameter name */
char *val,			/* value to set this parameter to */
int index			/* index of vector parameter; -1 if scalar */
)
{
/*
 * Parameter checked out OK.  Add it to linked list of param values.
 * Search through existing list to replace value if its there.
 * Otherwise, add it to the end of the list.
 */
    PARAM_LIST *ptr;             	/* loops through parameter list */
    PARAM_LIST *param;		/* parameter entry in list */

    /* printf("Debug: Adding parameter %s with value %s and index %d\n",
      name, val, index); */

    ptr = zz->Params;
    while (ptr != NULL) {
	if ((!strcmp(name, ptr->name)) && (index == ptr->index)){	
	    /* string and index match */
	    ZOLTAN_FREE(&(ptr->new_val));
	    ptr->new_val = val;
	    return (ZOLTAN_OK);
	}
	ptr = ptr->next;
    }

    /* This is a new parameter, add it to list. */
    param = (PARAM_LIST *) ZOLTAN_MALLOC(sizeof(PARAM_LIST));
    if (param == NULL) {
	ZOLTAN_FREE(&name);
	ZOLTAN_FREE(&val);
	return (ZOLTAN_MEMERR);
    }
    ptr = zz->Params;
    zz->Params = param;
    param->next = ptr;
    param->name = name;
    param->index = index;
    param->new_val = val;

    return (ZOLTAN_OK);
}

static int remove_param(
ZZ *zz,				/* Zoltan structure */
char *name,			/* parameter name */
int index 			/* index for vector param; -1 for scalars */
)
{
/*
 * Parameter checked out OK.  Remove it from linked list of param values.
 * If it is not in the list, do nothing.
 */
    PARAM_LIST *ptr, *oldptr;	/* loops through parameter list */

    oldptr = NULL;
    ptr = zz->Params;
    while (ptr != NULL) {
	if ((!strcmp(name, ptr->name)) && 
            ((index == ptr->index) || (index == -1))){
	    /* String and index match. (Index -1 matches anything.) */
            /* Remove parameter from list */
            if (oldptr == NULL)
               zz->Params = ptr->next;
            else
               oldptr->next = ptr->next;
            /* Free parameter */
            ZOLTAN_FREE(&(ptr->name));
            ZOLTAN_FREE(&(ptr->new_val));
            ZOLTAN_FREE(&ptr);
            /* Return OK */
	    return (ZOLTAN_OK);
	}
        oldptr = ptr;
	ptr = ptr->next;
    }

    /* Parameter was not in list */
    return (ZOLTAN_OK);
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
