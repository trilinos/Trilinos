/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/* This routine builds a description of the machine
   defined by the processes in the MPI communicator
   and the machine description file.
 */

#include "lb_const.h"
#include "ha_const.h"
#include "params_const.h"

int LB_Build_Machine_Desc(
   LB *lb              /* The load-balancing structure.                */
)
{
  int ierr = LB_OK;
  int use_mach_desc;
  char *filename;

  PARAM_VARS Mach_params[] = {
      { "USE_MACHINE_DESC", NULL, "INT" },
      { "MACHINE_DESC_FILE", NULL, "STRING" },
      { NULL, NULL, NULL } };
  
  Mach_params[0].ptr = (void *) &use_mach_desc;
  Mach_params[1].ptr = (void *) &filename;

  use_mach_desc = 0;
  filename = MACHINE_DESC_FILE_DEFAULT;

  LB_Assign_Param_Vals(lb->Params, Mach_params);

  if (use_mach_desc) {
    /* If lb->Machine_Desc already exists, don't rebuild it. */
    /* NB: Can the processors change after LB_Create_Obj ?? */

    if (lb->Machine_Desc == NULL){
      /* Read machine description from file. 
       * Use LB_Get_Processor_Name to extract the sub-machine
       * that this LB structure is running on. 
       */
      printf("Sorry, heterogenous load-balancing is still under development!\n");
      ierr = LB_WARN;
    }
  }
  else {
    lb->Machine_Desc = NULL;
  }

  return ierr;
}

