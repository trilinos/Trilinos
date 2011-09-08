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
#include "reftree_const.h"
#include "ha_const.h"
#include "rib_const.h"
#include "simple_const.h"
#include "hsfc_const.h"
#include "all_allo_const.h"
#include "order_const.h"
#include "phg_const.h"
#include "graph_const.h"
#ifdef ZOLTAN_HIER
#include "hier.h"
#endif
#ifdef ZOLTAN_DRUM
#include "ha_drum.h"
#endif
#include "coloring_const.h"
#include "zz_const.h"

static int add_param(ZZ *, char **, char **, int);
static int remove_param(ZZ *, char *, int);

/* List of set_parameter functions to be called */
static ZOLTAN_SET_PARAM_FN * Param_func[] = {
       Zoltan_Set_Malloc_Param,
       Zoltan_RCB_Set_Param,
       Zoltan_Third_Set_Param,
#if defined(ZOLTAN_PARMETIS) || defined(ZOLTAN_METIS)
       Zoltan_ParMetis_Set_Param,
#endif 
#if defined(ZOLTAN_SCOTCH) || defined(ZOLTAN_PTSCOTCH)
       Zoltan_Scotch_Set_Param,
#endif
#ifdef ZOLTAN_OCT
       Zoltan_Oct_Set_Param,
#endif
       Zoltan_Reftree_Set_Param,
       Zoltan_RIB_Set_Param,
       Zoltan_HSFC_Set_Param,
       Zoltan_Order_Set_Param,
       Zoltan_PHG_Set_Param,
#ifdef ZOLTAN_HIER
       Zoltan_Hier_Set_Param,
#endif
#ifdef ZOLTAN_DRUM
       Zoltan_Drum_Set_Param,
#endif
       Zoltan_ZG_Set_Param,
       /* Zoltan_Set_Machine_Param, */
       Zoltan_Color_Set_Param,
       /*** Add your new parameter setting function here! ***/
       Zoltan_Graph_Package_Set_Param,
       Zoltan_Random_Set_Param,

       NULL /* Last entry _must_ be NULL! */
};

int Zoltan_Set_Param(
ZZ *zz,				/* Zoltan structure */
const char *name1,		/* parameter name */
const char *val1)		/* value to set this parameter to */
{
    return Zoltan_Set_Param_Vec(zz, name1, val1, -1);
}

int Zoltan_Set_Param_Vec(
ZZ *zz,				/* Zoltan structure */
const char *name1,		/* parameter name */
const char *val1,		/* value to set this parameter to */
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
    else if (status == 2) {
	sprintf(msg, "Invalid value `%s' for parameter `%s'; default "
                     "value will be used.\n", val, name);
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
    	    add_param(zz, &name, &val, index); 	/* Add parameter to list */
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
char **name,			/* parameter name */
char **val,			/* value to set this parameter to */
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
      *name, *val, index); */

    ptr = zz->Params;
    while (ptr != NULL) {
	if ((!strcmp(*name, ptr->name)) && (index == ptr->index)){	
	    /* string and index match */
	    ZOLTAN_FREE(name);
	    ZOLTAN_FREE(&(ptr->new_val));
	    ptr->new_val = *val;
	    return (ZOLTAN_OK);
	}
	ptr = ptr->next;
    }

    /* This is a new parameter, add it to list. */
    param = (PARAM_LIST *) ZOLTAN_MALLOC(sizeof(PARAM_LIST));
    if (param == NULL) {
	ZOLTAN_FREE(name);
	ZOLTAN_FREE(val);
	return (ZOLTAN_MEMERR);
    }
    ptr = zz->Params;
    zz->Params = param;
    param->next = ptr;
    param->name = *name;
    param->index = index;
    param->new_val = *val;

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

int      Zoltan_Filter_Params(
struct Zoltan_Struct *to_zz,   /* add to this ... */
struct Zoltan_Struct *from_zz, /* ... parameters of interest found here */
PARAM_VARS * params,		/* parameters of interest */
int debug_level,                /* level for output of debugging info     */
int proc,                       /* processor # (controls debug printing)  */
int print_proc                  /* processor that should perform printing */
)
{	
    char     *name;		/* name of parameter being reset */
    char     *val;		/* new value for parameter       */
    int       index;		/* index of parameter entry      */
    int       found;		/* is name found?                */
    int       ierr;		/* error code                    */
    PARAM_LIST *from_list;
    PARAM_VARS *param_ptr;      /* pointer to current param      */

    ierr = ZOLTAN_OK;

    from_list = from_zz->Params;

    while (from_list != NULL) {
        param_ptr = params;
	name = from_list->name;
	val = from_list->new_val;
	index = from_list->index;

        if (debug_level > 0 && proc == print_proc){
          printf("Incoming parameter list: %s = %s\n",name,val);
        }

	found = 0;
	while (param_ptr->name != NULL) {
	    if (!strcmp(param_ptr->name, name)) {
		found = 1;
		break;
	    }
	    param_ptr++;
	}

	if (found) {		/* name found */

          Zoltan_Set_Param_Vec(to_zz, name, val, index);
          if (debug_level > 0 && proc == print_proc){
            if (index >= 0)
              printf("Put %s[%d] = %s in outgoing parameter list\n",name,index,val);
            else
              printf("Put %s = %s in outgoing parameter list\n",name,val);
          }
        }
        from_list = from_list->next;
    }


    return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
