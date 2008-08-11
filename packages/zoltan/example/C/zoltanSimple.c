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

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /*
  ** Initialize application data.  In this example,
  ** create a small test mesh and divide it across processors
  */
  
  exSetDivisions(32);    /* rectilinear mesh is div X div X div */
  
  MyNumPts = exInitializePoints(&Pts, &Gids, me, nprocs);
  
  /*  Initialize Zoltan */
  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    free(Pts); free(Gids);
    exit(0);
  }

  /******************************************************************
  ** Partition it using RCB
  ******************************************************************/
  /* Allocate and initialize memory for Zoltan structure */

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* Set general parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* Set RCB parameters */

  Zoltan_Set_Param(zz, "KEEP_CUTS", "1"); 
  Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"); 

  /* Register call-back query functions. */

  Zoltan_Set_Num_Obj_Fn(zz, exGetNumberOfAssignedObjects, NULL);
  Zoltan_Set_Obj_List_Fn(zz, exGetObjectList, NULL);
  Zoltan_Set_Num_Geom_Fn(zz, exGetObjectSize, NULL);
  Zoltan_Set_Geom_Multi_Fn(zz, exGetObjectCoords, NULL);

  /* Perform partitioning */
  rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
    &numImport, &importGlobalGids, &importLocalGids, &importProcs, &importToPart,
    &numExport, &exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);

  rc = exGlobalSuccess(rc, nprocs, me, (me == 0));

  /* Process partitioning results; 
  ** in this case, print information; 
  ** in a "real" application, migrate data here.
  */
  if (!rc){
    exPrintGlobalResult("Recursive Coordinate Bisection", nprocs, me,
                           MyNumPts, numImport, numExport, changes);
  }
  else{
    free(Pts);
    free(Gids);
    Zoltan_Destroy(&zz);
    MPI_Finalize();
    exit(0);
  }

  /*
  ** Clean up.
  */

  /* Free Zoltan memory allocated by Zoltan_LB_Partition. */
  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  /* Free Zoltan memory allocated by Zoltan_Create. */
  Zoltan_Destroy(&zz);

  /* Free Application memory */
  free(Pts); free(Gids);

  /**********************
  ** all done ***********
  **********************/

  MPI_Finalize();

  return 0;
}
