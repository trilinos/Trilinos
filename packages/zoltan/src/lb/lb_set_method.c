/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "zz_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines to set the load-balancing method.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Set_LB_Method(ZZ *zz, char *method_name)
{
/*
 *  Function to set the load balancing method to be used.
 *  Input:
 *    zz                 --  The Zoltan structure to which this method
 *                           applies.
 *    method_name        --  String specifying the desired method.
 *
 *  Output:
 *    zz->LB.*           --  Appropriate fields set to designated values.
 */

  char *yo = "Zoltan_LB_Set_LB_Method";
  char msg[256];
  char *method_upper;
  int error = ZOLTAN_OK;

  /*
   *  Compare method_name string with standard strings for methods.
   *  If a match is found, set zz->LB.Method and other pointers.
   *  But first free any left-over data from the previous method.
   */

  if (zz->LB.Free_Structure != NULL)
    zz->LB.Free_Structure(zz);

  /*
   *  Convert method_name to all upper case.
   *  Do not change the original string.
   */

  error = Zoltan_Clean_String(method_name, &method_upper);
  if (error) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "Error returned from Zoltan_Clean_String; No method set.");
    goto End;
  }

  if (strcmp(method_upper, "SIMPLE") == 0) {
    zz->LB.Method = SIMPLE;
    zz->LB.LB_Fn = Zoltan_Simple;
    zz->LB.Free_Structure = NULL;
    zz->LB.Copy_Structure = NULL;
  }
  else if (strcmp(method_upper, "RANDOM") == 0) {
    zz->LB.Method = RANDOM;
    zz->LB.LB_Fn = Zoltan_Random;
    zz->LB.Free_Structure = NULL;
    zz->LB.Copy_Structure = NULL;
  }
  else if (strcmp(method_upper, "RCB") == 0) {
    zz->LB.Method = RCB;
    zz->LB.LB_Fn = Zoltan_RCB;
    zz->LB.Free_Structure = Zoltan_RCB_Free_Structure;
    zz->LB.Copy_Structure = Zoltan_RCB_Copy_Structure;
    zz->LB.Point_Assign = Zoltan_RB_Point_Assign;
    zz->LB.Box_Assign = Zoltan_RB_Box_Assign;
  }
  else if (strcmp(method_upper, "OCTPART") == 0) {
#ifdef ZOLTAN_OCT
    zz->LB.Method = OCTPART;
    zz->LB.LB_Fn = Zoltan_Octpart;
    zz->LB.Free_Structure = Zoltan_Oct_Free_Structure;
    zz->LB.Copy_Structure = NULL;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
#else
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "OCTPART method selected but not compiled into Zoltan; "
                       "Compile with ZOLTAN_OCT=1.");
    error = ZOLTAN_FATAL;
    goto End;
#endif
  }
  else if ((strcmp(method_upper, "GRAPH") == 0)
           || (strcmp(method_upper, "PARMETIS") == 0)) {
    zz->LB.Method = PARMETIS;
    zz->LB.LB_Fn = Zoltan_ParMetis;
    zz->LB.Free_Structure = NULL;
    zz->LB.Copy_Structure = NULL;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
  }
  else if (strcmp(method_upper, "JOSTLE") == 0) {
    zz->LB.Method = JOSTLE;
    zz->LB.LB_Fn = Zoltan_Jostle;
    zz->LB.Free_Structure = NULL;
    zz->LB.Copy_Structure = NULL;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
  }
  else if (strcmp(method_upper, "REFTREE") == 0) {
    zz->LB.Method = REFTREE;
    zz->LB.LB_Fn = Zoltan_Reftree_Part;
    zz->LB.Free_Structure = Zoltan_Reftree_Free_Structure;
    zz->LB.Copy_Structure = NULL;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
  }
  else if (strcmp(method_upper, "RIB") == 0) {
    zz->LB.Method = RIB;
    zz->LB.LB_Fn = Zoltan_RIB;
    zz->LB.Free_Structure = Zoltan_RIB_Free_Structure;
    zz->LB.Copy_Structure = Zoltan_RIB_Copy_Structure;
    zz->LB.Point_Assign = Zoltan_RB_Point_Assign;
    zz->LB.Box_Assign = Zoltan_RB_Box_Assign;
  }
  else if (strcmp(method_upper, "HSFC") == 0) {
    zz->LB.Method = HSFC;
    zz->LB.LB_Fn = Zoltan_HSFC;
    zz->LB.Free_Structure = Zoltan_HSFC_Free_Structure;
    zz->LB.Copy_Structure = Zoltan_HSFC_Copy_Structure;
    zz->LB.Point_Assign = Zoltan_HSFC_Point_Assign;
    zz->LB.Box_Assign = Zoltan_HSFC_Box_Assign;
  }
  else if ((strcmp(method_upper, "HYPERGRAPH") == 0) 
           || (strcmp(method_upper, "PHG") == 0)
           || (strcmp(method_upper, "PHG_REPART") == 0)
	   || (strcmp(method_upper, "PHG_REFINE") == 0)
	   || (strcmp(method_upper, "PHG_MULTILEVEL_REFINE") == 0)           
           || (strcmp(method_upper, "PATOH") == 0)
           || (strcmp(method_upper, "PARKWAY") == 0)){
    /* The hypergraph methods have a lot in common. We allow
       the user to either set the LB method to HYPERGRAPH and
       select a package via the parameter HYPERGRAPH_PACKAGE,
       or to set the LB method to be the package name (e.g. PHG, Patoh).
       Either way, LB.Method will eventually be set to the
       desired hypergraph package.
    */
#ifdef ZOLTAN_HG
    if (!strcmp(method_upper, "PATOH"))
      zz->LB.Method = PATOH;
    else if (!strcmp(method_upper, "PARKWAY"))
      zz->LB.Method = PARKWAY;
    else if (!strcmp(method_upper, "PHG_REPART"))
      zz->LB.Method = PHG_REPART;
    else if (!strcmp(method_upper, "PHG_REFINE"))
      zz->LB.Method = PHG_REFINE;
    else if (!strcmp(method_upper, "PHG_MULTILEVEL_REFINE"))
      zz->LB.Method = PHG_MULTILEVEL_REFINE;
    else /* HYPERGRAPH or PHG */
      zz->LB.Method = PHG;
    zz->LB.LB_Fn = Zoltan_PHG;
    zz->LB.Free_Structure = Zoltan_PHG_Free_Structure;
    zz->LB.Copy_Structure = NULL;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
#else
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Hypergraph method selected but not compiled into Zoltan; "
                       "Compile with ZOLTAN_HG=1.");
    error = ZOLTAN_FATAL;
    goto End;
#endif
  }
  else if (strcmp(method_upper, "HIER") == 0) {
#ifdef ZOLTAN_HIER
    zz->LB.Method = HIER;
    zz->LB.LB_Fn = Zoltan_Hier;
    zz->LB.Free_Structure = Zoltan_Hier_Free_Structure;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
#else
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "HIER method selected but not compiled into Zoltan; "
                       "Compile with ZOLTAN_HIER=1.");
    error = ZOLTAN_FATAL;
    goto End;
#endif
  }
  else if (strcmp(method_upper, "NONE") == 0) {
    zz->LB.Method = NONE;
    zz->LB.LB_Fn = NULL;
    zz->LB.Free_Structure = NULL;
    zz->LB.Copy_Structure = NULL;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
  }

  
  /*
   *  SET OTHER METHODS HERE!!
   */

  else {  
    sprintf(msg, "Invalid LB method specified:  %s\n", method_name);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    error = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->Proc == zz->Debug_Proc && zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS) {
    printf("ZOLTAN Load balancing method = %d (%s)\n", 
           zz->LB.Method, method_name);
  }

End:

  ZOLTAN_FREE(&method_upper);

  return (error);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
