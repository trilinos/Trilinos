#include <stdio.h>
#include "lb_const.h"


void      LB_Print_Params(
LB *lb)				/* load balance object */
{
/*
 *  Function to print out list of set parameter values.
 */

    LB_PARAM *ptr;

    ptr = lb->Params;
    printf("Parameter Settings\n");
    while (ptr != NULL) {
       printf("%s = %s\n",ptr->name, ptr->new_val);
       ptr = ptr->next;
    }
    printf("\n");
}
