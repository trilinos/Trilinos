#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "lb_const.h"
#include "all_allo_const.h"


void LB_Free_Params(
LB *lb)				/* load balance object */
{
/*
 * Free the list of new parameter values.
 */
    LB_PARAM *ptr, *ptr2;	/* loops through parameter list */


    ptr = lb->Params;
    while (ptr != NULL) {
	ptr2 = ptr->next;
	LB_Free((void **) &(ptr->name));
	LB_Free((void **) &(ptr->new_val));
	LB_Free((void **) &ptr);
	ptr = ptr2;
    }

    lb->Params = NULL;
}
