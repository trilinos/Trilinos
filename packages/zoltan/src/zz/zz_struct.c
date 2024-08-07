// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "zz_rand.h"
#include "lb_init_const.h"
#include "params_const.h"
#include "ha_const.h"
#include "zz_util_const.h"

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
int proc;

  MPI_Comm_rank(communicator, &proc);

  /*
   * Allocate storage for the Zoltan structure.
   */

  zz = (ZZ *) ZOLTAN_MALLOC(sizeof(ZZ));
  if (!zz) {
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

ZZ *Zoltan_Copy(ZZ const *from)
{
  int fail=0;

  ZZ *to = Zoltan_Create(from->Communicator);

  fail = Zoltan_Copy_To(to, from);

  if (fail) {
    Zoltan_Destroy(&to);
  }

  return to; 
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_Copy_To(ZZ *to, ZZ const *from)
{
  /*
   * Copy one Zoltan_Struct to another.  "to" must be a valid 
   * Zoltan_Struct.
   */

  Zoltan_Free_Zoltan_Struct_Members(to);
  MPI_Comm_free(&(to->Communicator));

  *to = *from;

  MPI_Comm_dup(from->Communicator, &(to->Communicator));
  
  to->Params = NULL;
  Zoltan_Copy_Params(&(to->Params), from->Params);

  to->ZTime = Zoltan_Timer_Copy(from->ZTime);

  memset(&(to->LB), 0, sizeof(struct Zoltan_LB_Struct));
  Zoltan_LB_Copy_Struct(to, from);

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
  Zoltan_Free_Params(&(zz->Params));
  Zoltan_Timer_Destroy(&(zz->ZTime));
  Zoltan_Free_Structures(zz);  /* Algorithm-specific structures */
  Zoltan_LB_Free_Struct(&(zz->LB));
  Zoltan_Order_Free_Struct(&(zz->Order));
  Zoltan_TPL_Order_Free_Struct(&(zz->TPL_Order));
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
  zz->Seed = ZOLTAN_RAND_INIT;
  zz->Timer = ZOLTAN_TIMER_DEF;
  zz->Params = NULL;
  zz->Deterministic = ZOLTAN_DETERMINISTIC_DEF;
  zz->Obj_Weight_Dim = ZOLTAN_OBJ_WEIGHT_DEF;
  zz->Edge_Weight_Dim = ZOLTAN_EDGE_WEIGHT_DEF;

  zz->Get_Part_Multi = NULL;
  zz->Get_Part = NULL;
  zz->Get_Num_Edges_Multi = NULL;
  zz->Get_Num_Edges = NULL;
  zz->Get_Edge_List_Multi = NULL;
  zz->Get_Edge_List = NULL;
  zz->Get_HG_Size_CS = NULL;
  zz->Get_HG_CS = NULL;
  zz->Get_HG_Size_Edge_Wts = NULL;
  zz->Get_HG_Edge_Wts = NULL;
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
  zz->Get_Num_Fixed_Obj = NULL;
  zz->Get_Fixed_Obj_List = NULL;
  zz->Get_Part_Fort = NULL;
  zz->Get_Num_Edges_Fort = NULL;
  zz->Get_Edge_List_Fort = NULL;
  zz->Get_HG_Size_CS_Fort = NULL;
  zz->Get_HG_CS_Fort = NULL;
  zz->Get_HG_Size_Edge_Wts_Fort = NULL;
  zz->Get_HG_Edge_Wts_Fort = NULL;
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
  zz->Get_Num_Fixed_Obj_Fort = NULL;
  zz->Get_Fixed_Obj_List_Fort = NULL;

  zz->Get_Part_Data = NULL;
  zz->Get_Num_Edges_Data = NULL;
  zz->Get_Edge_List_Data = NULL;
  zz->Get_HG_Size_CS_Data = NULL;
  zz->Get_HG_CS_Data = NULL;
  zz->Get_HG_Size_Edge_Wts_Data = NULL;
  zz->Get_HG_Edge_Wts_Data = NULL;
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
  zz->Get_Num_Fixed_Obj_Data = NULL;
  zz->Get_Fixed_Obj_List_Data = NULL;

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

  zz->Get_Hier_Num_Levels = NULL;
  zz->Get_Hier_Part = NULL;
  zz->Get_Hier_Method = NULL;
  zz->Get_Hier_Num_Levels_Fort = NULL;
  zz->Get_Hier_Part_Fort = NULL;
  zz->Get_Hier_Method_Fort = NULL;
  zz->Get_Hier_Num_Levels_Data = NULL;
  zz->Get_Hier_Part_Data = NULL;
  zz->Get_Hier_Method_Data = NULL;

  zz->Order.needfree = 0;
  zz->TPL_Order.needfree = 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
size_t Zoltan_Serialize_Size(struct Zoltan_Struct const *zz) 
{
  /* Compute size of buffer needed to serialize Zoltan_Struct */
  size_t bufSize = 0;
  /* Copy 11 integers from Zoltan_Struct */
  bufSize += sizeof(int) * 11;
  bufSize += Zoltan_Serialize_Params_Size(zz);
  bufSize += Zoltan_LB_Serialize_Size(zz);
  return bufSize;
}

/*****************************************************************************/
int Zoltan_Serialize(ZZ const *zz, size_t bufSize, char *buf)
{
  int ierr = ZOLTAN_OK;
  char *bufptr = buf;

  /* Copy 11 integers; if add more, update Zoltan_Serialize_Size */
  int *intptr = (int *) bufptr;
  *intptr = zz->Num_GID; intptr++;
  *intptr = zz->Num_LID; intptr++;
  *intptr = zz->Debug_Level; intptr++;
  *intptr = zz->Debug_Proc; intptr++;
  *intptr = zz->Fortran; intptr++;
  *intptr = zz->Tflops_Special; intptr++;
  *intptr = *((int*) &(zz->Seed)); intptr++;
  *intptr = zz->Deterministic; intptr++;
  *intptr = zz->Obj_Weight_Dim; intptr++;
  *intptr = zz->Edge_Weight_Dim; intptr++;
  *intptr = zz->Timer; intptr++;
  bufptr = (char *) intptr;

  /* Need some params to set pointers on receiving side; advance bufptr */
  Zoltan_Serialize_Params(zz, &bufptr);

  /* Serialize LB information; advance bufptr */
  Zoltan_LB_Serialize(zz, &bufptr);

  if (bufptr - buf > bufSize) {
    ZOLTAN_PRINT_ERROR(zz->Proc, "Zoltan_Serialize", 
                       "Buffer provided is too small; Zoltan just overwrote "
                       "your memory.");
    ierr = ZOLTAN_MEMERR;
  }
  return ierr;
}

/*****************************************************************************/
int Zoltan_Deserialize(struct Zoltan_Struct *zz, size_t bufSize, char *buf)
{
  int ierr = ZOLTAN_OK;
  char *bufptr = buf;
  
  /* Copy 11 integers */
  int *intptr = (int *) bufptr;
  zz->Num_GID = *intptr; intptr++;
  zz->Num_LID = *intptr; intptr++;
  zz->Debug_Level = *intptr; intptr++;
  zz->Debug_Proc = *intptr; intptr++;
  zz->Fortran = *intptr; intptr++;
  zz->Tflops_Special = *intptr; intptr++;
  zz->Seed = *((unsigned int *)(intptr)); intptr++;
  zz->Deterministic = *intptr; intptr++;
  zz->Obj_Weight_Dim = *intptr; intptr++;
  zz->Edge_Weight_Dim = *intptr; intptr++;
  zz->Timer = *intptr; intptr++;
  bufptr = (char *) intptr;

  /* Consistent with Zoltan defaults, set default LB_METHOD to RCB; doing so
   * sets the various function pointers for RCB. 
   * Method and function pointers will be reset if parameter LB_METHOD has
   * been provided by user. */
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Deserialize_Params(zz, &bufptr);

  /* Deserialize LB information; advance bufptr */
  Zoltan_LB_Deserialize(zz, &bufptr);

  if (bufptr - buf > bufSize) {
    ZOLTAN_PRINT_ERROR(zz->Proc, "Zoltan_Deserialize", 
                       "Buffer provided is too small; Zoltan just copied "
                       "garbage memory.");
    ierr = ZOLTAN_MEMERR;
  }

  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
