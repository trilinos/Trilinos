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
#include "ha_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines that create, copy, and destroy load-balancing 
 *  structures (struct Zoltan_Struct).
 *  These functions are all callable by the application. 
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Free_Structures(ZZ *);
static void Zoltan_Init(ZZ *);
static void Zoltan_Reset(ZZ *);
static void Zoltan_Free_Zoltan_Struct_Members(ZZ *);

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

  Zoltan_Init(zz);

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

  Zoltan_LB_Init(&(zz->LB), zz->Num_Proc);
  Zoltan_Migrate_Init(&(zz->Migrate));

  zz->ZTime = Zoltan_Timer_Create(ZOLTAN_TIMER_DEF);

  return(zz);
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

ZZ *Zoltan_Copy(ZZ *from)
{
  int fail=0;

  ZZ *to = Zoltan_Create(from->Communicator);

  fail = Zoltan_Copy_To(to, from);

  if (fail)
    {
    Zoltan_Destroy(&to);
    }

  return to; 
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

#define COPY_FIELD(f) to->f = from->f;

int Zoltan_Copy_To(ZZ *to, ZZ *from)
{
  /*
   * Copy one Zoltan_Struct to another.  "to" must be a valid 
   * Zoltan_Struct.
   */

  Zoltan_Reset(to);

  MPI_Comm_dup(from->Communicator, &(to->Communicator));

  COPY_FIELD(Proc)
  COPY_FIELD(Num_Proc)
  COPY_FIELD(Num_GID)
  COPY_FIELD(Num_LID)
  COPY_FIELD(Debug_Level)
  COPY_FIELD(Debug_Proc)
  COPY_FIELD(Fortran)
  COPY_FIELD(Tflops_Special)
  COPY_FIELD(Deterministic)
  COPY_FIELD(Obj_Weight_Dim)
  COPY_FIELD(Edge_Weight_Dim)
  COPY_FIELD(Timer)

  COPY_FIELD(Get_Partition_Multi)
  COPY_FIELD(Get_Partition_Multi_Fort)
  COPY_FIELD(Get_Partition_Multi_Data)
  COPY_FIELD(Get_Partition)
  COPY_FIELD(Get_Partition_Fort)
  COPY_FIELD(Get_Partition_Data)

  COPY_FIELD(Get_Num_Edges)
  COPY_FIELD(Get_Num_Edges_Fort)
  COPY_FIELD(Get_Num_Edges_Data)
  COPY_FIELD(Get_Num_Edges_Multi)
  COPY_FIELD(Get_Num_Edges_Multi_Fort)
  COPY_FIELD(Get_Num_Edges_Multi_Data)

  COPY_FIELD(Get_Edge_List)
  COPY_FIELD(Get_Edge_List_Fort)
  COPY_FIELD(Get_Edge_List_Data)
  COPY_FIELD(Get_Edge_List_Multi)
  COPY_FIELD(Get_Edge_List_Multi_Fort)
  COPY_FIELD(Get_Edge_List_Multi_Data)

  COPY_FIELD(Get_Num_Geom)
  COPY_FIELD(Get_Num_Geom_Fort)
  COPY_FIELD(Get_Num_Geom_Data)

  COPY_FIELD(Get_Geom_Multi)
  COPY_FIELD(Get_Geom_Multi_Fort)
  COPY_FIELD(Get_Geom_Multi_Data)

  COPY_FIELD(Get_Geom)
  COPY_FIELD(Get_Geom_Fort)
  COPY_FIELD(Get_Geom_Data)

  COPY_FIELD(Get_Num_Obj)
  COPY_FIELD(Get_Num_Obj_Fort)
  COPY_FIELD(Get_Num_Obj_Data)

  COPY_FIELD(Get_Obj_List)
  COPY_FIELD(Get_Obj_List_Fort)
  COPY_FIELD(Get_Obj_List_Data)
  COPY_FIELD(Get_First_Obj)
  COPY_FIELD(Get_First_Obj_Fort)
  COPY_FIELD(Get_First_Obj_Data)
  COPY_FIELD(Get_Next_Obj)
  COPY_FIELD(Get_Next_Obj_Fort)
  COPY_FIELD(Get_Next_Obj_Data)

  COPY_FIELD(Get_Num_Border_Obj)
  COPY_FIELD(Get_Num_Border_Obj_Fort)
  COPY_FIELD(Get_Num_Border_Obj_Data)
  COPY_FIELD(Get_Border_Obj_List)
  COPY_FIELD(Get_Border_Obj_List_Fort)
  COPY_FIELD(Get_Border_Obj_List_Data)
  COPY_FIELD(Get_First_Border_Obj)
  COPY_FIELD(Get_First_Border_Obj_Fort)
  COPY_FIELD(Get_First_Border_Obj_Data)
  COPY_FIELD(Get_Next_Border_Obj)
  COPY_FIELD(Get_Next_Border_Obj_Fort)
  COPY_FIELD(Get_Next_Border_Obj_Data)

  COPY_FIELD(Get_Num_Coarse_Obj)
  COPY_FIELD(Get_Num_Coarse_Obj_Fort)
  COPY_FIELD(Get_Num_Coarse_Obj_Data)

  COPY_FIELD(Get_Coarse_Obj_List)
  COPY_FIELD(Get_Coarse_Obj_List_Fort)
  COPY_FIELD(Get_Coarse_Obj_List_Data)
  COPY_FIELD(Get_First_Coarse_Obj)
  COPY_FIELD(Get_First_Coarse_Obj_Fort)
  COPY_FIELD(Get_First_Coarse_Obj_Data)
  COPY_FIELD(Get_Next_Coarse_Obj)
  COPY_FIELD(Get_Next_Coarse_Obj_Fort)
  COPY_FIELD(Get_Next_Coarse_Obj_Data)

  COPY_FIELD(Get_Num_Child)
  COPY_FIELD(Get_Num_Child_Fort)
  COPY_FIELD(Get_Num_Child_Data)

  COPY_FIELD(Get_Child_List)
  COPY_FIELD(Get_Child_List_Fort)
  COPY_FIELD(Get_Child_List_Data)

  COPY_FIELD(Get_Child_Weight)
  COPY_FIELD(Get_Child_Weight_Fort)
  COPY_FIELD(Get_Child_Weight_Data)

  COPY_FIELD(Get_Num_HG_Edges)
  COPY_FIELD(Get_Num_HG_Edges_Fort)
  COPY_FIELD(Get_Num_HG_Edges_Data)

  COPY_FIELD(Get_HG_Edge_List)
  COPY_FIELD(Get_HG_Edge_List_Fort)
  COPY_FIELD(Get_HG_Edge_List_Data)

  COPY_FIELD(Get_Num_HG_Pins)
  COPY_FIELD(Get_Num_HG_Pins_Fort)
  COPY_FIELD(Get_Num_HG_Pins_Data)

  COPY_FIELD(Get_Obj_Size)
  COPY_FIELD(Get_Obj_Size_Fort)
  COPY_FIELD(Get_Obj_Size_Data)
  COPY_FIELD(Get_Obj_Size_Multi)
  COPY_FIELD(Get_Obj_Size_Multi_Fort)
  COPY_FIELD(Get_Obj_Size_Multi_Data)

  COPY_FIELD(Pack_Obj)
  COPY_FIELD(Pack_Obj_Fort)
  COPY_FIELD(Pack_Obj_Data)
  COPY_FIELD(Pack_Obj_Multi)
  COPY_FIELD(Pack_Obj_Multi_Fort)
  COPY_FIELD(Pack_Obj_Multi_Data)

  COPY_FIELD(Unpack_Obj)
  COPY_FIELD(Unpack_Obj_Fort)
  COPY_FIELD(Unpack_Obj_Data)
  COPY_FIELD(Unpack_Obj_Multi)
  COPY_FIELD(Unpack_Obj_Multi_Fort)
  COPY_FIELD(Unpack_Obj_Multi_Data)

  COPY_FIELD(Get_Processor_Name)
  COPY_FIELD(Get_Processor_Name_Data)

  Zoltan_Copy_Machine_Desc(&(to->Machine_Desc), from->Machine_Desc);
  
  Zoltan_Copy_Params(&(to->Params), from->Params);

  Zoltan_Timer_Copy(&(to->ZTime), from->ZTime);

  Zoltan_LB_Copy_Struct(to, from);

  Zoltan_Migrate_Copy_Struct(&(to->Migrate), &(from->Migrate));

  return 0;
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

    Zoltan_Free_Zoltan_Struct_Members(*zz);

    MPI_Comm_free(&((*zz)->Communicator));

    ZOLTAN_FREE(zz);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void Zoltan_Free_Zoltan_Struct_Members(ZZ *zz)
{
  Zoltan_Free_Machine_Desc(&(zz->Machine_Desc));
  Zoltan_Free_Params(&(zz->Params));
  Zoltan_Timer_Destroy(&(zz->ZTime));
  Zoltan_Free_Structures(zz);  /* Algorithm-specific structures */
  Zoltan_LB_Free_Struct(&(zz->LB));
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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Init(ZZ* zz)
{
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

  zz->Get_Partition_Multi = NULL;
  zz->Get_Partition = NULL;
  zz->Get_Num_Edges_Multi = NULL;
  zz->Get_Num_Edges = NULL;
  zz->Get_Edge_List_Multi = NULL;
  zz->Get_Edge_List = NULL;
  zz->Get_Num_HG_Edges = NULL;
  zz->Get_HG_Edge_List = NULL;
  zz->Get_Num_Geom = NULL;
  zz->Get_Geom_Multi = NULL;
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

  zz->Get_Partition_Fort = NULL;
  zz->Get_Num_Edges_Fort = NULL;
  zz->Get_Edge_List_Fort = NULL;
  zz->Get_Num_Geom_Fort = NULL;
  zz->Get_Geom_Multi_Fort = NULL;
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

  zz->Get_Partition_Data = NULL;
  zz->Get_Num_Edges_Data = NULL;
  zz->Get_Edge_List_Data = NULL;
  zz->Get_Num_HG_Edges_Data = NULL;
  zz->Get_HG_Edge_List_Data = NULL;
  zz->Get_Num_Geom_Data = NULL;
  zz->Get_Geom_Data = NULL;
  zz->Get_Num_Obj_Data = NULL;
  zz->Get_Obj_List_Data = NULL;
  zz->Get_First_Obj_Data = NULL;
  zz->Get_Next_Obj_Data = NULL;
  zz->Get_Num_Border_Obj_Data = NULL;
  zz->Get_Border_Obj_List_Data = NULL;
  zz->Get_First_Border_Obj_Data = NULL;
  zz->Get_Next_Border_Obj_Data = NULL;
  zz->Get_Num_Coarse_Obj_Data = NULL;
  zz->Get_Coarse_Obj_List_Data = NULL;
  zz->Get_First_Coarse_Obj_Data = NULL;
  zz->Get_Next_Coarse_Obj_Data = NULL;
  zz->Get_Num_Child_Data = NULL;
  zz->Get_Child_List_Data = NULL;
  zz->Get_Child_Weight_Data = NULL;

  zz->Pack_Obj = NULL;
  zz->Unpack_Obj = NULL;
  zz->Get_Obj_Size = NULL;
  zz->Pack_Obj_Multi = NULL;
  zz->Unpack_Obj_Multi = NULL;
  zz->Get_Obj_Size_Multi = NULL;
  
  zz->Pack_Obj_Fort = NULL;
  zz->Unpack_Obj_Fort = NULL;
  zz->Get_Obj_Size_Fort = NULL;

  zz->Pack_Obj_Data = NULL;
  zz->Unpack_Obj_Data = NULL;
  zz->Get_Obj_Size_Data = NULL;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void Zoltan_Reset(ZZ *zz)
{
/*
 *  Function to reset fields of a Zoltan_Struct back to their
 *  initial state. 
 */
  Zoltan_Free_Zoltan_Struct_Members(zz);
  MPI_Comm_free(&(zz->Communicator));

  Zoltan_Init(zz);
  Zoltan_LB_Init(&(zz->LB), 0);
  Zoltan_Migrate_Init(&(zz->Migrate));

  zz->Communicator = MPI_COMM_NULL;
  zz->Proc = -1;
  zz->Num_Proc = 0;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
