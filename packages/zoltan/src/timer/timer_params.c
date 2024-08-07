// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <string.h>

#include "timer.h"
#include "zoltan_timer.h"
#include "params_const.h"
#include "zoltan_util.h"
#include "zz_const.h"

/* Function to interpret and set timer parameters */

int Zoltan_Set_Timer_Param(
const char *name,                     /* input:  name of variable */
const char *val,                      /* input:  value of variable */
int *timer)                     /* output: timer type */
{
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    int status;
    PARAM_VARS Timer_params[] = {
        { "TIMER", NULL, "STRING", 0 },
        { NULL, NULL, NULL, 0 }
    };
    char *yo = "Zoltan_Set_Timer_Param";

    (*timer) = ZOLTAN_TIME_WALL;  /* default timer value */

    status = Zoltan_Check_Param(name, val, Timer_params, &result, &index);

    if (status == 0 && index == 0) {
        if (!strcmp(result.sval, "WALL"))
          (*timer) = ZOLTAN_TIME_WALL;
        else if (strcmp(result.sval, "CPU")==0) {
          (*timer) = ZOLTAN_TIME_CPU;
        }
        else if (strcmp(result.sval, "USER")==0){
#ifndef NO_TIMES
          (*timer) = ZOLTAN_TIME_USER;
#else
          ZOLTAN_PRINT_WARN(-1, yo, "User time not available;"
                          " CPU clock time will be used instead.");
          (*timer) = ZOLTAN_TIME_CPU;
#endif
        }
        else{
          char msg[256];
          sprintf(msg, "Unknown timer option %s.", result.sval);
          ZOLTAN_PRINT_WARN(-1, yo, msg);
          status = 2; /* Illegal parameter */
        }
    }

    return(status);
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
