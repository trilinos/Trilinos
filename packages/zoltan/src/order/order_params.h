
#ifndef __ORDER_PARAMS_H
#define __ORDER_PARAMS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "params_const.h"

/**********  parameters structure for ordering **********/
static PARAM_VARS Order_params[] = {
        { "ORDER_METHOD", NULL, "STRING", 0 },
        { "USE_ORDER_INFO", NULL, "INT", 0 },
        { NULL, NULL, NULL, 0 } };


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
