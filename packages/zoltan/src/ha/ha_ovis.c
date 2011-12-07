/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000-2012, Sandia National Laboratories.                    *
 *****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "zz_util_const.h"
#include "params_const.h"
#include "ha_ovis.h"

#ifdef ZOLTAN_OVIS

/***************************************************************/
/* Placeholder for OVIS interface.  Subject to massive change! */
/***************************************************************/


/*  Parameters structure for OVIS.  */
/*  Define parameter name and type here; must keep NULL line at end. */
static PARAM_VARS OVIS_params[] = {
                  { "OVIS_HELLO", NULL, "STRING", 0 },
                  { "OVIS_OUTPUT_LEVEL", NULL, "INT", 0 },
                  { "OVIS_MINVERSION", NULL, "DOUBLE", 0 },
                  { NULL, NULL, NULL, 0 } };


int Zoltan_OVIS_Setup(
  ZZ *zz        /* Zoltan structure; needed for parameters */
) 
{
  /* Allow OVIS parameters to be passed to Zoltan via Zoltan_Set_Param */
  /* Three example parameters below. */

  /* Declare variables for parameter values; initialize to default values */
  char   ovis_hello[MAX_PARAM_STRING_LEN];
  int    ovis_output_level = 1;
  double ovis_minversion = 0.1;
  strcpy(ovis_hello, "Howdy, Shorty!");

  /* Tell Zoltan to associate parameter names with the variables. */
  Zoltan_Bind_Param(OVIS_params, "OVIS_HELLO", 
                    ovis_hello);
  Zoltan_Bind_Param(OVIS_params, "OVIS_OUTPUT_LEVEL", 
                    (void *) &ovis_output_level);
  Zoltan_Bind_Param(OVIS_params, "OVIS_MINVERSION", 
                    (void *) &ovis_minversion);

  /* Tell Zoltan to look for parameters matching the names above */
  /*Zoltan_Assign_Param_Vals(zz->Params, OVIS_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);*/
  Zoltan_Assign_Param_Vals(zz->Params, OVIS_params, 7, zz->Proc,
                           zz->Debug_Proc);

  return ZOLTAN_OK;
}


/****************************************************************************/
/** Utility code; should not need to edit the code below to add parameters.**/
/****************************************************************************/

int Zoltan_OVIS_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */

    status = Zoltan_Check_Param(name, val, OVIS_params, &result, &index);

    return(status);
}

#endif /* ZOLTAN_OVIS */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
