#include "params_const.h"

int LB_Set_Machine_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    PARAM_VARS Mach_params[] = {
        { "USE_MACHINE_DESC", NULL, "INT" },
        { "MACHINE_DESC_FILE", NULL, "STRING" },
        { NULL, NULL, NULL } };

    status = LB_Check_Param(name, val, Mach_params, &result, &index);

    return(status);
}

