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
#include "lb_init_const.h"
#include "params_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines that create and destroy load-balancing 
 *  structures (struct Zoltan_Struct).
 *  These functions are all callable by the application. 
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Free_Structures(ZZ *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

ZZ *Zoltan_Create(MPI_Comm communicator)
{
/*
 *  Function to create a Zoltan structure.  May want more than one
 *  structure if using different decompositions with different techniques.
 *  This function allocates and initializes the structure.
 *  Output:
 *    ZZ *               --  Pointer to a Zoltan structure.
 *
 */

char *yo = "Zoltan_Create";
ZZ *zz;

  /*
   * Allocate storage for the Zoltan structure.
   */

  zz = (ZZ *) ZOLTAN_MALLOC(sizeof(ZZ));
  if (!zz) {
    int proc;
    MPI_Comm_rank(communicator, &proc);
    ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory to create structure.");
    return NULL;
  }

  /*
   *  Set MPI values for zz:
   */

  if (communicator == MPI_COMM_NULL) {
    /*
     *  The processor is not in the communicator for the load-balancing
     *  structure.  Set zz->Communicator to MPI_COMM_NULL and give dummy 
     *  values to zz->Proc and zz->Num_Proc.
     */
    zz->Communicator = MPI_COMM_NULL;
    zz->Proc = -1;
    zz->Num_Proc = 0;
  }
  else {
    /*
     *  Set Communicator to the communicator passed in.
     */
    MPI_Comm_dup(communicator, &(zz->Communicator));
    MPI_Comm_size(zz->Communicator, &(zz->Num_Proc));
    MPI_Comm_rank(zz->Communicator, &(zz->Proc));
  }

  /*
   *  Set defaults for fields of lb:
   */

  zz->Num_GID = ZOLTAN_NUM_ID_ENTRIES_DEF;
  zz->Num_LID = ZOLTAN_NUM_ID_ENTRIES_DEF;
  zz->Debug_Level = ZOLTAN_DEBUG_LEVEL_DEF;
  zz->Debug_Proc = ZOLTAN_DEBUG_PROC_DEF;
  zz->Fortran = 0;
  zz->Tflops_Special = ZOLTAN_TFLOPS_SPECIAL_DEF;
  zz->Timer = ZOLTAN_TIMER_DEF;
  zz->Machine_Desc = NULL;
  zz->Params = NULL;
  zz->Deterministic = ZOLTAN_DETERMINISTIC_DEF;
  zz->Obj_Weight_Dim = ZOLTAN_OBJ_WEIGHT_DEF;
  zz->Edge_Weight_Dim = ZOLTAN_EDGE_WEIGHT_DEF;

  zz->Get_Num_Edges = NULL;
  zz->Get_Edge_List = NULL;
  zz->Get_Num_Geom = NULL;
  zz->Get_Geom = NULL;
  zz->Get_Num_Obj = NULL;
  zz->Get_Obj_List = NULL;
  zz->Get_First_Obj = NULL;
  zz->Get_Next_Obj = NULL;
  zz->Get_Num_Border_Obj = NULL;
  zz->Get_Border_Obj_List = NULL;
  zz->Get_First_Border_Obj = NULL;
  zz->Get_Next_Border_Obj = NULL;
  zz->Get_Num_Coarse_Obj = NULL;
  zz->Get_Coarse_Obj_List = NULL;
  zz->Get_First_Coarse_Obj = NULL;
  zz->Get_Next_Coarse_Obj = NULL;
  zz->Get_Num_Child = NULL;
  zz->Get_Child_List = NULL;
  zz->Get_Child_Weight = NULL;

  zz->Get_Num_Edges_Fort = NULL;
  zz->Get_Edge_List_Fort = NULL;
  zz->Get_Num_Geom_Fort = NULL;
  zz->Get_Geom_Fort = NULL;
  zz->Get_Num_Obj_Fort = NULL;
  zz->Get_Obj_List_Fort = NULL;
  zz->Get_First_Obj_Fort = NULL;
  zz->Get_Next_Obj_Fort = NULL;
  zz->Get_Num_Border_Obj_Fort = NULL;
  zz->Get_Border_Obj_List_Fort = NULL;
  zz->Get_First_Border_Obj_Fort = NULL;
  zz->Get_Next_Border_Obj_Fort = NULL;
  zz->Get_Num_Coarse_Obj_Fort = NULL;
  zz->Get_Coarse_Obj_List_Fort = NULL;
  zz->Get_First_Coarse_Obj_Fort = NULL;
  zz->Get_Next_Coarse_Obj_Fort = NULL;
  zz->Get_Num_Child_Fort = NULL;
  zz->Get_Child_List_Fort = NULL;
  zz->Get_Child_Weight_Fort = NULL;

  zz->Pack_Obj = NULL;
  zz->Unpack_Obj = NULL;
  zz->Get_Obj_Size = NULL;
  zz->Pack_Obj_Multi = NULL;
  zz->Unpack_Obj_Multi = NULL;
  zz->Get_Obj_Size_Multi = NULL;
  
  zz->Pack_Obj_Fort = NULL;
  zz->Unpack_Obj_Fort = NULL;
  zz->Get_Obj_Size_Fort = NULL;

  Zoltan_LB_Init(&(zz->LB));
  Zoltan_Migrate_Init(&(zz->Migrate));

  return(zz);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void Zoltan_Destroy(ZZ **zz)
{
/*
 *  Function to free a Zoltan structure.
 *  Input:
 *    ZZ **zz           --  Pointer to a Zoltan structure.
 *
 */

  if (*zz != NULL) {

    Zoltan_Free_Structures(*zz);

    Zoltan_Free_Params(&((*zz)->Params));

    MPI_Comm_free(&((*zz)->Communicator));

    ZOLTAN_FREE(zz);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


static void Zoltan_Free_Structures(
  ZZ *zz)				/* Zoltan structure */
{
/*
 * Free any persistent memory associated with Zoltan modules.
 */

  /* Free load-balancing data */
  if (zz->LB.Free_Structure != NULL) 
    zz->LB.Free_Structure(zz);

 /* Add calls to additional module-specific free routines here.  */
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
