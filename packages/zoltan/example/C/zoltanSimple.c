/*
** $Id$
**
** Basic example of using Zoltan to compute an RCB partitioning
*/

#include "zoltan.h"
#include "exzoltan.h"

#include <stdlib.h>
#include <stdio.h>

static int MyNumPts;
static int *Gids;
static float *Pts;

int main(int argc, char *argv[])
{
  int rc, me, nprocs;
  float ver;
  struct Zoltan_Struct *zz;
  int changes;
  int numGidEntries;
  int numLidEntries;
  int numImport;
  ZOLTAN_ID_PTR importGlobalGids;
  ZOLTAN_ID_PTR importLocalGids;
  int *importProcs;
  int *importToPart;
  int numExport;
  ZOLTAN_ID_PTR exportGlobalGids;
  ZOLTAN_ID_PTR exportLocalGids; 
  int *exportProcs;
  int *exportToPart;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /*
  ** Create a small test mesh and divide it across processors
  */
  
  exSetDivisions(32);    /* rectilinear mesh is div X div X div */
  
  MyNumPts = exInitializePoints(&Pts, &Gids, me, nprocs);
  
  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    free(Pts); free(Gids);
    exit(0);
  }

  /******************************************************************
  ** Partition it using RCB
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* RCB parameters */

  Zoltan_Set_Param(zz, "KEEP_CUTS", "1"); 
  Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"); 

  /* Call backs */

  Zoltan_Set_Num_Obj_Fn(zz, exGetNumberOfAssignedObjects, NULL);
  Zoltan_Set_Obj_List_Fn(zz, exGetObjectList, NULL);
  Zoltan_Set_Num_Geom_Fn(zz, exGetObjectSize, NULL);
  Zoltan_Set_Geom_Multi_Fn(zz, exGetObject, NULL);

  rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
    &numImport, &importGlobalGids, &importLocalGids, &importProcs, &importToPart,
    &numExport, &exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);

  rc = exGlobalSuccess(rc, nprocs, me, (me == 0));

  if (!rc){
    exPrintGlobalResult("Recursive Coordinate Bisection", nprocs, me,
                           MyNumPts, numImport, numExport, changes);
  }
  else{
    MPI_Finalize();
    free(Pts);
    free(Gids);
    Zoltan_Destroy(&zz);
    exit(0);
  }

  /*
  ** Clean up
  */

  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);

  free(Pts); free(Gids);

  /**********************
  ** all done ***********
  **********************/

  MPI_Finalize();

  return 0;
}
