/*
** $Id$
**
** Example of using Zoltan's RCB, RIB and HSFC partitioning algorithms,
**   and also using a Zoltan_Timer.
*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "zoltan.h"
#include "mpi.h"
#include "exzoltan.h"

/*
** Include these headers if we are using a Zoltan_Timer
*/
#include "zoltan_timer.h"
#include "zz_const.h"
/*
************************ 
*/

#define NTIMERS 3

static int timers[NTIMERS];
static char *events[NTIMERS]={
"RCB calculation",
"RIB calculation",
"HSFC calculation"
};
enum {rcb, rib, hsfc};

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

  if (getenv("WAIT_FOR_DEBUGGER"))
    sleep(60);

  /*
  ** Create a small test mesh and divide it across processors.
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

  timers[rcb] = Zoltan_Timer_Init(zz->ZTime, 1, events[rcb]);

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


  ZOLTAN_TIMER_START(zz->ZTime, timers[rcb], zz->Communicator);
      
  rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
    &numImport, &importGlobalGids, &importLocalGids, &importProcs, &importToPart,
    &numExport, &exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);

  ZOLTAN_TIMER_STOP(zz->ZTime, timers[rcb], zz->Communicator);

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

  MPI_Barrier(MPI_COMM_WORLD);

  /******************************************************************
  ** Partition it using RIB
  ******************************************************************/

  timers[rib] = Zoltan_Timer_Init(zz->ZTime, 1, events[rib]);

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "RIB");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* RIB parameters */

  Zoltan_Set_Param(zz, "KEEP_CUTS", "1"); 
  Zoltan_Set_Param(zz, "RIB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param(zz, "AVERAGE_CUTS", "1"); 

  /* Call backs */

  Zoltan_Set_Num_Obj_Fn(zz, exGetNumberOfAssignedObjects, NULL);
  Zoltan_Set_Obj_List_Fn(zz, exGetObjectList, NULL);
  Zoltan_Set_Num_Geom_Fn(zz, exGetObjectSize, NULL);
  Zoltan_Set_Geom_Multi_Fn(zz, exGetObject, NULL);

  ZOLTAN_TIMER_START(zz->ZTime, timers[rib], zz->Communicator);

  rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
    &numImport, &importGlobalGids, &importLocalGids, &importProcs, &importToPart,
    &numExport, &exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);

  ZOLTAN_TIMER_STOP(zz->ZTime, timers[rib], zz->Communicator);

  rc = exGlobalSuccess(rc, nprocs, me, (me == 0));

  if (!rc){
    exPrintGlobalResult("Recursive Inertial Bisection", nprocs, me,
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

  MPI_Barrier(MPI_COMM_WORLD);

  /******************************************************************
  ** Partition it using HSFC
  ******************************************************************/

  timers[hsfc] = Zoltan_Timer_Init(zz->ZTime, 1, events[hsfc]);

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "HSFC");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* HSFC parameters */

  Zoltan_Set_Param(zz, "KEEP_CUTS", "1"); 

  /* Call backs */

  Zoltan_Set_Num_Obj_Fn(zz, exGetNumberOfAssignedObjects, NULL);
  Zoltan_Set_Obj_List_Fn(zz, exGetObjectList, NULL);
  Zoltan_Set_Num_Geom_Fn(zz, exGetObjectSize, NULL);
  Zoltan_Set_Geom_Multi_Fn(zz, exGetObject, NULL);

  ZOLTAN_TIMER_START(zz->ZTime, timers[hsfc], zz->Communicator);

  rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
    &numImport, &importGlobalGids, &importLocalGids, &importProcs, &importToPart,
    &numExport, &exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);

  ZOLTAN_TIMER_STOP(zz->ZTime, timers[hsfc], zz->Communicator);

  rc = exGlobalSuccess(rc, nprocs, me, (me == 0));

  if (!rc){
    exPrintGlobalResult("Hilbert Space Filling Curve", nprocs, me,
                           MyNumPts, numImport, numExport, changes);
  }
  else{
    MPI_Finalize();
    free(Pts);
    free(Gids);
    Zoltan_Destroy(&zz);
    exit(0);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (me == 0){
    Zoltan_Timer_PrintAll(zz->ZTime, 0, stdout);
  }

  /*
  ** Clean up
  */

  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);

  MPI_Barrier(MPI_COMM_WORLD);

  free(Pts); free(Gids);

  /**********************
  ** all done ***********
  **********************/

  MPI_Finalize();

  return 0;
}
