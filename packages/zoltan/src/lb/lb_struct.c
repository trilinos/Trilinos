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
#include "all_allo_const.h"
#include "params_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines that create and destroy load-balancing 
 *  structures (struct LB_Struct).
 *  These functions are all callable by the application. 
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

LB *LB_Create(MPI_Comm communicator)
{
/*
 *  Function to create a load balancing structure.  May want more than one
 *  structure if using different decompositions with different techniques.
 *  This function allocates and initializes the structure.
 *  Output:
 *    LB *               --  Pointer to a LB structure.
 *
 */

char *yo = "LB_Create_Object";
LB *lb;

  /*
   * Allocate storage for the load-balancing structure.
   */

  lb = (LB *) LB_MALLOC(sizeof(LB));
  if (!lb) {
    fprintf(stderr, "Error from %s: Insufficient memory\n", yo);
    return NULL;
  }

  /*
   *  Set MPI values for lb:
   */

  if (communicator == MPI_COMM_NULL) {
    /*
     *  The processor is not in the communicator for the load-balancing
     *  structure.  Set lb->Communicator to MPI_COMM_NULL and give dummy 
     *  values to lb->Proc and lb->Num_Proc.
     */
    lb->Communicator = MPI_COMM_NULL;
    lb->Proc = -1;
    lb->Num_Proc = 0;
  }
  else {
    /*
     *  Set Communicator to the communicator passed in.
     */
    MPI_Comm_dup(communicator, &(lb->Communicator));
    MPI_Comm_size(lb->Communicator, &(lb->Num_Proc));
    MPI_Comm_rank(lb->Communicator, &(lb->Proc));
  }

  /*
   *  Set defaults for fields of lb:
   */

  lb->Method = RCB;    
  lb->LB_Fn = LB_rcb;
  lb->Debug_Level = LB_DEBUG_LEVEL_DEF;
  lb->Debug_Proc = LB_DEBUG_PROC_DEF;
  lb->Fortran = 0;
  lb->Timer = LB_TIMER_DEF;
  lb->Machine_Desc = NULL;
  lb->Params = NULL;
  lb->Imbalance_Tol = LB_IMBALANCE_TOL_DEF;
  lb->Deterministic = LB_DETERMINISTIC_DEF;
  lb->Obj_Weight_Dim = LB_OBJ_WEIGHT_DEF;
  lb->Comm_Weight_Dim = LB_COMM_WEIGHT_DEF;
  lb->Data_Structure = NULL;

  lb->Get_Num_Edges = NULL;
  lb->Get_Edge_List = NULL;
  lb->Get_Num_Geom = NULL;
  lb->Get_Geom = NULL;
  lb->Get_Num_Obj = NULL;
  lb->Get_Obj_List = NULL;
  lb->Get_First_Obj = NULL;
  lb->Get_Next_Obj = NULL;
  lb->Get_Num_Border_Obj = NULL;
  lb->Get_Border_Obj_List = NULL;
  lb->Get_First_Border_Obj = NULL;
  lb->Get_Next_Border_Obj = NULL;
  lb->Get_Num_Coarse_Obj = NULL;
  lb->Get_Coarse_Obj_List = NULL;
  lb->Get_First_Coarse_Obj = NULL;
  lb->Get_Next_Coarse_Obj = NULL;
  lb->Get_Num_Child = NULL;
  lb->Get_Child_List = NULL;
  lb->Get_Child_Weight = NULL;

  lb->Get_Num_Edges_Fort = NULL;
  lb->Get_Edge_List_Fort = NULL;
  lb->Get_Num_Geom_Fort = NULL;
  lb->Get_Geom_Fort = NULL;
  lb->Get_Num_Obj_Fort = NULL;
  lb->Get_Obj_List_Fort = NULL;
  lb->Get_First_Obj_Fort = NULL;
  lb->Get_Next_Obj_Fort = NULL;
  lb->Get_Num_Border_Obj_Fort = NULL;
  lb->Get_Border_Obj_List_Fort = NULL;
  lb->Get_First_Border_Obj_Fort = NULL;
  lb->Get_Next_Border_Obj_Fort = NULL;
  lb->Get_Num_Coarse_Obj_Fort = NULL;
  lb->Get_Coarse_Obj_List_Fort = NULL;
  lb->Get_First_Coarse_Obj_Fort = NULL;
  lb->Get_Next_Coarse_Obj_Fort = NULL;
  lb->Get_Num_Child_Fort = NULL;
  lb->Get_Child_List_Fort = NULL;
  lb->Get_Child_Weight_Fort = NULL;

  lb->Migrate.Auto_Migrate = LB_AUTO_MIGRATE_DEF;
  lb->Migrate.Pre_Migrate = NULL;
  lb->Migrate.Mid_Migrate = NULL;
  lb->Migrate.Post_Migrate = NULL;
  lb->Migrate.Pack_Obj = NULL;
  lb->Migrate.Unpack_Obj = NULL;
  lb->Migrate.Get_Obj_Size = NULL;
  
  lb->Migrate.Pre_Migrate_Fort = NULL;
  lb->Migrate.Mid_Migrate_Fort = NULL;
  lb->Migrate.Post_Migrate_Fort = NULL;
  lb->Migrate.Pack_Obj_Fort = NULL;
  lb->Migrate.Unpack_Obj_Fort = NULL;
  lb->Migrate.Get_Obj_Size_Fort = NULL;

  return(lb);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Destroy(LB **lb)
{
/*
 *  Function to free a load balancing structure.
 *  Input:
 *    LB ** lb           --  Pointer to a LB structure.
 *
 */

  LB_Free_Structure(*lb);

  LB_Free_Params(&((*lb)->Params));

  MPI_Comm_free(&((*lb)->Communicator));

  LB_FREE(lb);
}
