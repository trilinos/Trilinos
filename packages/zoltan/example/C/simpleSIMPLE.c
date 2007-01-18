#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "zoltan.h"

static int numObjects = 36;
static int numProcs, myRank, myNumObj;
static int myGlobalIDs[36];

static float objWeight(int globalID)
{
float w;
  if (globalID % numProcs == 0){
    w = 3;  /* simulate an initial imbalance */
  else
    w = (globalID % 3 + 1);

  return w;
}
static int get_number_of_objects(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return myNumObj;
}
static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
int i;
  *ierr = ZOLTAN_OK;
  for (i=0; i<myNumObj; i++){
    globalID[i] = myGlobalIDs[i];
    if (obj_wgts){
      obj_wgts[i] = objWeight(myGlobalIDs[i]);
    }
  }
}

int main(int argc, char *argv[])
{
  int rc, i;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids;
  ZOLTAN_ID_PTR exportGlobalGids, exportLocalGids; 
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  int ngids;
  int *gid_flags, *gid_list;
  float wgt;
  int j, nextIdx;

  /******************************************************************
  ** Initialize MPI and Zoltan
  ******************************************************************/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    printf("Error initializing Zoltan\n");
    MPI_Finalize();
    exit(0);
  }

  /******************************************************************
  ** Create a simple initial partitioning for this example
  ******************************************************************/

  for (i=0, myNumObj=0; i<numObjects; i++){
    if (i % numProcs == myRank){
      myGlobalIDs[myNumObj++] = i+1;
    }
  }

  /******************************************************************
  ** Create a Zoltan library structure for this instance of load
  ** balancing.  Set the parameters and query functions.
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "LB_METHOD", "SIMPLE");  /* Zoltan method: "SIMPLE" */
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); /* global ID is 1 integer */
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0"); /* no local IDs */
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* weights are 1 float */

  /* Query functions */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_object_list, NULL);

  /******************************************************************
  ** Call Zoltan to partition the objects.
  ******************************************************************/

  rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partitioning was changed, 0 otherwise */ 
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of vertices to be sent to me */
        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
        &importLocalGids,   /* Local IDs of vertices to be sent to me */
        &importProcs,    /* Process rank for source of each incoming vertex */
        &importToPart,   /* New partition for each incoming vertex */
        &numExport,      /* Number of vertices I must send to other processes*/
        &exportGlobalGids,  /* Global IDs of the vertices I must send */
        &exportLocalGids,   /* Local IDs of the vertices I must send */
        &exportProcs,    /* Process to which I send each of the vertices */
        &exportToPart);  /* Partition to which each vertex will belong */

  if (rc != ZOLTAN_OK){
    printf("Error in Zoltan library\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  /******************************************************************
  ** Visualize the new partitioning.
  ** Create a list of GIDs now assigned to my partition, let
  ** process zero display the partitioning.
  ******************************************************************/

  ngids = get_number_of_objects(NULL, &rc);
  gid_flags = (int *)calloc(sizeof(int) , numObjects);
  gid_list = (int *)malloc(sizeof(int) * ngids);
  get_object_list(NULL, 1, 1,
                  (ZOLTAN_ID_PTR)gid_list, NULL, 1, NULL, &rc);

  for (i=0; i <numProcs; i++){
    if (i == myRank){
      if (!i) printf("\nInitial Partitioning:\n========================\n");
      wgt = 0.0;
      printf("%d: ",i);
      for (j=0; j<ngids; j++){
        if (j && (j % 20 == 0)) printf("\n   ");
        printf("%d (%2.0f) ",gid_list[j],objWeight(gid_list[j]));
        wgt += objWeight(gid_list[j]);
      }
      printf("   weight: %f\n",wgt);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  for (i=0; i<ngids; i++){
    gid_flags[gid_list[i]-1] = 1;    /* my original vertices */
  }
  for (i=0; i<numImport; i++){ 
    gid_flags[importGlobalGids[i] - 1] = 1;  /* my imports */
  }
  for (i=0; i<numExport; i++){
    gid_flags[exportGlobalGids[i] - 1] = 0;  /* my exports */
  }
  nextIdx = 0;
  for (i=0; i<numObjects; i++){
    if (gid_flags[i]){   
      gid_flags[nextIdx] = i+1; /* my new GID list */ 
      nextIdx++;
    }
  }
  for (i=0; i <numProcs; i++){
    if (i == myRank){
      if (!i) printf("\nNew Partitioning:\n=========================\n");
      printf("%d: ",i);
      wgt=0.0;
      for (j=0; j<nextIdx; j++){
        if (j && (j % 20 == 0)) printf("\n   ");
        printf("%d (%2.0f) ",gid_flags[j],objWeight(gid_flags[j]));
        wgt += objWeight(gid_flags[j]);
      }
      printf("   weight: %f\n",wgt);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  if (gid_flags) free(gid_flags);
  if (gid_list) free(gid_list);

  /******************************************************************
  ** Free the arrays allocated by Zoltan_LB_Partition, and free
  ** the storage allocated for the Zoltan structure.
  ******************************************************************/

  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);

  MPI_Finalize();

  return 0;
}
