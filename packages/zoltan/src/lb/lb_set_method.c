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

int Zoltan_LB_Set_Method(LB *lb, char *method_name)
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

  char *yo = "Zoltan_LB_Set_Method";
  char msg[256];
  char *method_upper;
  int error;

  /*
   *  Compare method_name string with standard strings for methods.
   *  If a match is found, set lb->Method and other pointers.
   *  But first free any left-over data from the previous method.
   */

  Zoltan_Free_Structure(lb);

  /*
   *  Convert method_name to all upper case.
   *  Do not change the original string.
   */

  error = Zoltan_Clean_String(method_name, &method_upper);
  if (error) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
      "Error returned from Zoltan_Clean_String; No method set.");
    ZOLTAN_FREE(&method_upper);
    return error;
  }

  if (strcmp(method_upper, "RCB") == 0) {
    lb->Method = RCB;
    lb->LB_Fn = Zoltan_RCB;
  }
  else if (strcmp(method_upper, "OCTPART") == 0) {
    lb->Method = OCTPART;
    lb->LB_Fn = Zoltan_Octpart;
  }
  else if (strcmp(method_upper, "PARMETIS") == 0) {
    lb->Method = PARMETIS;
    lb->LB_Fn = Zoltan_ParMetis;
  }
  else if (strcmp(method_upper, "JOSTLE") == 0) {
    lb->Method = JOSTLE;
    lb->LB_Fn = Zoltan_Jostle;
  }
  else if (strcmp(method_upper, "REFTREE") == 0) {
    lb->Method = REFTREE;
    lb->LB_Fn = Zoltan_Reftree_Part;
  }
  else if (strcmp(method_upper, "RIB") == 0) {
    lb->Method = RIB;
    lb->LB_Fn = Zoltan_RIB;
  }
  else if (strcmp(method_upper, "SFC") == 0) {
    lb->Method = SFC;
    lb->LB_Fn = Zoltan_SFC;
  }
  else if (strcmp(method_upper, "HSFC") == 0) {
    lb->Method = HSFC;
    lb->LB_Fn = Zoltan_HSFC;
  }
  else if (strcmp(method_upper, "NONE") == 0) {
    lb->Method = NONE;
    lb->LB_Fn = NULL;
  }

  /*
   *  SET OTHER METHODS HERE!!
   */

  else {  
    sprintf(msg, "Invalid LB method specified:  %s\n", method_name);
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, msg);
    ZOLTAN_FREE(&method_upper);
    return (ZOLTAN_FATAL);
  }

  if (lb->Proc == lb->Debug_Proc && lb->Debug_Level >= ZOLTAN_DEBUG_PARAMS) {
    printf("ZOLTAN Load balancing method = %d (%s)\n", lb->Method, method_name);
  }

  ZOLTAN_FREE(&method_upper);

  return (ZOLTAN_OK);
}
