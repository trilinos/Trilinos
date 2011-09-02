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

  if (strcmp(method_upper, "BLOCK") == 0) {
    zz->LB.Method = BLOCK;
    zz->LB.LB_Fn = Zoltan_Block;
    zz->LB.Free_Structure = NULL;
    zz->LB.Copy_Structure = NULL;
  }
  else if (strcmp(method_upper, "CYCLIC") == 0) {
    zz->LB.Method = CYCLIC;
    zz->LB.LB_Fn = Zoltan_Cyclic;
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
                     "Configure (via CMAKE) with -DZoltan_ENABLE_OCT:BOOL=ON.");
    error = ZOLTAN_FATAL;
    goto End;
#endif
  }
  else if (strcmp(method_upper, "GRAPH") == 0){
    zz->LB.Method = GRAPH;
    zz->LB.LB_Fn = Zoltan_Graph;
    /* Next two are useful only when using PHG */
    zz->LB.Free_Structure = Zoltan_PHG_Free_Structure;
    zz->LB.Copy_Structure = Zoltan_PHG_Copy_Structure;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
  }
  /* PARMETIS and JOSTLE are here for backward compatibility.
   * New way: LB_METHOD = GRAPH
   *          GRAPH_PACKAGE = PARMETIS or JOSTLE or PHG
   */
  else if (strcmp(method_upper, "PARMETIS") == 0){
#ifdef ZOLTAN_PARMETIS
    zz->LB.Method = GRAPH;
    zz->LB.LB_Fn = Zoltan_ParMetis;
    zz->LB.Free_Structure = NULL;
    zz->LB.Copy_Structure = NULL;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
#else
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "ParMETIS method selected but "
                       "ParMETIS not compiled into Zoltan; "
                       "Compile with --with-parmetis.");
    error = ZOLTAN_FATAL;
    goto End;
#endif
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
           || (strcmp(method_upper, "PHG") == 0)){

    /* HYPERGRAPH is a family of methods. */
    /* PHG is Zoltan's standard parallel hypergraph partitioner. */
    zz->LB.Method = HYPERGRAPH;
    zz->LB.LB_Fn = Zoltan_PHG;
    zz->LB.Free_Structure = Zoltan_PHG_Free_Structure;
    zz->LB.Copy_Structure = Zoltan_PHG_Copy_Structure;
    zz->LB.Point_Assign = NULL;
    zz->LB.Box_Assign = NULL;
  }
  else if (strcmp(method_upper, "HIER") == 0) {
#ifdef ZOLTAN_HIER
    zz->LB.Method = HIER;
    zz->LB.LB_Fn = Zoltan_Hier;
    zz->LB.Free_Structure = Zoltan_Hier_Free_Structure;
    zz->LB.Copy_Structure = Zoltan_Hier_Copy_Structure;
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
