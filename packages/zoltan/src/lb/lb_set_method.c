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
#include "lb_util_const.h"
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines to set the load-balancing method.
 *  These functions are all callable by the application.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Set_Method(LB *lb, char *method_name)
{
/*
 *  Function to set the load balancing method to be used.
 *  Input:
 *    lb                 --  The load balancing structure to which this method
 *                           applies.
 *    method_name        --  String specifying the desired method.
 *
 *  Output:
 *    lbf*               --  Appropriate fields set to designated values.
 */

  char *yo = "LB_Set_Method";
  char *method_upper;
  int error;

  /*
   *  Compare method_name string with standard strings for methods.
   *  If a match is found, set lb->Method and other pointers.
   *  But first free any left-over data from the previous method.
   */

  LB_Free_Structure(lb);

  /*
   *  Convert method_name to all upper case.
   *  Do not change the original string.
   */

  error = LB_clean_string(method_name, &method_upper);
  if (error) {
    fprintf(stderr, "%s Error %d returned from LB_clean_string; "
                    "No method set.\n", yo, error);
    LB_FREE(&method_upper);
    return error;
  }

  if (strcmp(method_upper, "RCB") == 0) {
    lb->Method = RCB;
    lb->LB_Fn = LB_rcb;
  }
  else if (strcmp(method_upper, "OCTPART") == 0) {
    lb->Method = OCTPART;
    lb->LB_Fn = LB_octpart;
  }
  else if (strcmp(method_upper, "PARMETIS") == 0) {
    lb->Method = PARMETIS;
    lb->LB_Fn = LB_ParMetis;
  }
  else if (strcmp(method_upper, "JOSTLE") == 0) {
    lb->Method = JOSTLE;
    lb->LB_Fn = LB_Jostle;
  }
  else if (strcmp(method_upper, "REFTREE") == 0) {
    lb->Method = REFTREE;
    lb->LB_Fn = LB_Reftree_Part;
  }
  else if (strcmp(method_upper, "IRB") == 0) {
    lb->Method = IRB;
    lb->LB_Fn = LB_irb;
  }
  else if (strcmp(method_upper, "NONE") == 0) {
    lb->Method = NONE;
    lb->LB_Fn = NULL;
  }

  /*
   *  SET OTHER METHODS HERE!!
   */

  else {  
    fprintf(stderr, "Error from %s:  Invalid LB method specified:  %s\n", 
            yo, method_name);
    return (LB_FATAL);
  }

  if (lb->Proc == lb->Debug_Proc && lb->Debug_Level >= LB_DEBUG_PARAMS) {
    printf("ZOLTAN Load balancing method = %d (%s)\n", lb->Method, method_name);
  }

  LB_FREE(&method_upper);

  return (LB_OK);
}
