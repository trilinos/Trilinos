/***************************************************************
** Basic example of using Zoltan to compute an RCB partitioning
** of a very simple mesh or graph.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "zoltan.h"
#include "simpleGraph.h"

static int myRank;
static int numProcs;

/* Definition of query functions to be registered with the Zoltan library */

int get_number_of_objects(void *data, int *ierr)
{
int i, numobj;

  /* Return the number of vertices I own.  Initially, we divide
   * the vertices in round robin fashion based on process rank.
   * The vertices are the objects that Zoltan will balance
   * across partitions.
   */

  numobj = 0;

  for (i=0; i<numvertices; i++){
    if (i % numProcs == myRank) numobj++;
  }

  *ierr = ZOLTAN_OK;

  return numobj;
}
void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
int numNeighbors;
int i, j, next;

  /* Sanity check - By setting parameters, I previously told the
   * Zoltan library how long my global and local IDs are, and
   * whether I am weighting the vertices, and if so, how many
   * floats I am using to weight each vertex. */

  if ( (sizeGID != 1) ||  /* My global IDs are 1 integer */
       (sizeLID != 1) ||  /* My local IDs are 1 integer */
       (wgt_dim != 1)){   /* My vertex weights are 1 float */
    *ierr = ZOLTAN_FATAL;
    return;
  }

  /* Return the list of vertices I own, and their weights.  I 
   * will assign a weight to each vertex equal to the number
   * of neighbors it has in the graph.
   */

  for (i=0, next=0; i<numvertices; i++){
    if (i % numProcs == myRank){
      globalID[next] = i+1;   /* application wide global ID */
      localID[next] = next;   /* process specific local ID */

      numNeighbors = 0;
      for (j=0; j<5; j++){
        if (edges[i][j] == -1) break;
        numNeighbors++;
      }
      obj_wgts[next] = (float)numNeighbors;
      next++;
    }
  }

  *ierr = ZOLTAN_OK;

  return;
}
int get_num_geometry(void *data, int *ierr)
{
  /* Return the number of values needed to express the geometry of a vertex */

  *ierr = ZOLTAN_OK;
  return 2;
}
void get_geometry_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr)
{
int i;
int row, col;

  /* Sanity check */

  if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 2)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  /* Return the coordinates of the requested vertices */

  for (i=0;  i < num_obj ; i++){
    row = (globalID[i] - 1) / 5;
    col = (globalID[i] - 1) % 5;

    geom_vec[2*i] = (double)col;
    geom_vec[2*i + 1] = (double)row;
  }

  *ierr = ZOLTAN_OK;

  return;
}

int main(int argc, char *argv[])
{
  int rc;
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

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* RCB parameters */

  Zoltan_Set_Param(zz, "KEEP_CUTS", "1"); 
  Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"); 

  /* Query functions */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_object_list, NULL);
  Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, NULL);
  Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, NULL);

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
