/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


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
