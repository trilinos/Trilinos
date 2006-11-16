/***************************************************************
** The most basic example of using Zoltan. 
** The objects have only a global ID and a weight. 
** The partitioning method is SIMPLE.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "zoltan.h"
#include "simpleGraph.h"
#include "simpleQueries.h"

int main(int argc, char *argv[])
{
  int rc, i, ngids, nextIdx;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  float *wgt_list;
  int *gid_flags, *gid_list, *lid_list;

  /******************************************************************
  ** Initialize MPI and Zoltan
  ******************************************************************/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    MPI_Finalize();
    exit(0);
  }

  /******************************************************************
  ** Create a Zoltan library structure for this instance of load
  ** balancing.  Set the parameters and query functions that will
  ** govern the library's calculation.  See the Zoltan User's
  ** Guide for the definition of these and many other parameters.
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "LB_METHOD", "SIMPLE");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* Query functions - defined in simpleQueries.h, to return
   * information about objects defined in simpleGraph.h      */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_object_list, NULL);

  /******************************************************************
  ** Zoltan can now partition the vertices in the simple mesh.
  ** In this simple example, we assume the number of partitions is
  ** equal to the number of processes.  Process rank 0 will own
  ** partition 0, process rank 1 will own partition 1, and so on.
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
    printf("sorry...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  /******************************************************************
  ** In a real application, you would rebalance the problem now by
  ** sending the objects to their new partitions.
  ******************************************************************/

  /******************************************************************
  ** Visualize the new partitioning.
  ** Create a list of GIDs now assigned to my partition, let
  ** process zero display the partitioning.
  ******************************************************************/

  ngids = get_number_of_objects(NULL, &rc);
  gid_flags = (int *)calloc(sizeof(int) , simpleNumVertices);
  gid_list = (int *)malloc(sizeof(int) * ngids);
  lid_list = (int *)malloc(sizeof(int) * ngids);
  wgt_list = (float *)malloc(sizeof(float) * simpleNumVertices);
  get_object_list(NULL, 1, 1,
                  (ZOLTAN_ID_PTR)gid_list, (ZOLTAN_ID_PTR)lid_list,
                  1, wgt_list, &rc);

  draw_partitions("initial distribution", ngids, gid_list, 1, wgt_list, 0);

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
  for (i=0; i<simpleNumVertices; i++){
    if (gid_flags[i]){
      gid_flags[nextIdx] = i+1; /* my new GID list */
      wgt_list[nextIdx] = simpleNumEdges[i];
      nextIdx++;
    }
  }
  draw_partitions("new partitioning", nextIdx, gid_flags, 1, wgt_list, 0);

  if (gid_flags) free(gid_flags);
  if (gid_list) free(gid_list);
  if (lid_list) free(lid_list);
  if (wgt_list) free(wgt_list);

  /******************************************************************
  ** Free the arrays allocated by Zoltan_LB_Partition, and free
  ** the storage allocated for the Zoltan structure.
  ******************************************************************/

  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);

  /**********************
  ** all done ***********
  **********************/

  MPI_Finalize();

  return 0;
}
