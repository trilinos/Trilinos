/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "ha_const.h"
#include "params_const.h"

/* 
 * These routines build a description of the machine
 * defined by the processes in the MPI communicator
 * and the machine description file.
 */

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/****** Parameters structure for building machine description. *****/

static PARAM_VARS Mach_params[] = {
        { "USE_MACHINE_DESC", NULL, "INT" },
        { "MACHINE_DESC_FILE", NULL, "STRING" },
        { NULL, NULL, NULL } };

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

int Zoltan_Set_Machine_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
  int status;
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */

  status = Zoltan_Check_Param(name, val, Mach_params, &result, &index);

  return(status);
}

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

int Zoltan_Build_Machine_Desc(
   ZZ *zz              /* The Zoltan structure.                */
)
{
  char *yo = "Zoltan_Build_Machine_Desc";
  int ierr = ZOLTAN_OK;
  int use_mach_desc;
  char filename[256];

  Zoltan_Bind_Param(Mach_params, "USE_MACHINE_DESC", (void *) &use_mach_desc);
  Zoltan_Bind_Param(Mach_params, "MACHINE_DESC_FILE", (void *) filename);

  use_mach_desc = 0;
  strcpy(filename, MACHINE_DESC_FILE_DEFAULT);

  Zoltan_Assign_Param_Vals(zz->Params, Mach_params, zz->Debug_Level, zz->Proc,
                       zz->Debug_Proc);

  if (use_mach_desc > 0) {
    /* If zz->Machine_Desc already exists, don't rebuild it
     * unless USE_MACHINE_DESC has been set to 2. 
     */

    if ((zz->Machine_Desc == NULL) || (use_mach_desc==2)){
      /* Read machine description from file. 
       * Use Zoltan_Get_Processor_Name to extract the sub-machine
       * on which this Zoltan structure is running. 
       * Broadcast the machine structure to all procs.
       */
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Sorry, heterogeneous load-balancing "
                                  "is still under development!");
      ierr = ZOLTAN_WARN;
    }
  }
  else {
    zz->Machine_Desc = NULL;
  }

  return ierr;
}

