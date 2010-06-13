/**************************************************************
* Stress test that can create a very large hypergraph to test
* the large memory problems.
*
* Argument is number of vertices in hypergraph.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <signal.h>
#include "zoltan.h"

static int myRank, numProcs, numPins, nborCount;
static int *vertex_part = NULL;

/* saves memory to use unit weights */
#define UNIT_WEIGHTS yes

/* Graph (Zoltan will create a hypergraph from it) */

#define NUM_GLOBAL_VERTICES     2500000000
#define VERTEX_WEIGHT_DIMENSION 1
#define EDGE_WEIGHT_DIMENSION 1

static ZOLTAN_GNO_TYPE numGlobalVertices;
static int vwgt_dim, ewgt_dim;
static ZOLTAN_ID_TYPE *vertexGIDs=NULL;
#ifndef UNIT_WEIGHTS
static float *vwgts=NULL;
#endif
static int edgeCutCost;

static int numMyVertices;

static int create_a_graph();
extern void Zoltan_write_linux_meminfo(int, char *);

#define proc_vertex_gid(proc, lid) (vertexGIDs[proc] + lid)

void meminfo_signal_handler(int sig)
{
  char msg[128];

  sprintf(msg,"(%d) Received signal %d\n",myRank,sig);

  // Signal handler for Linux that helps us to understand
  // whether failure was due to insufficient memory.

  signal(SIGINT, SIG_IGN);
  signal(SIGTERM, SIG_IGN);
  signal(SIGABRT, SIG_IGN);
  signal(SIGSEGV, SIG_IGN);
  signal(SIGFPE, SIG_IGN);

  Zoltan_write_linux_meminfo(1, msg);

  exit(sig);
}


/* Zoltan hypergraph query functions. */

static int get_number_of_vertices(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return numMyVertices;
}

static void get_vertex_list(void *data, int sizeGID, int sizeLID,
                  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
  int i;
  *ierr = ZOLTAN_OK;

  for (i=0; i < numMyVertices; i++){
    globalID[i] = proc_vertex_gid(myRank, i);
    localID[i] = i;
#ifndef UNIT_WEIGHTS
    obj_wgts[i] = vwgts[i];
#else
    obj_wgts[i] = 1.0;
#else
#endif
  }
}

static void get_local_hypergraph_size(void *data, 
                  int *num_lists, int *num_pins, int *format, int *ierr)
{
  /* Each vertex represents a "hyperedge".  
     It and its neighbors are the vertices in the hyperedge. */

  *ierr = ZOLTAN_OK;
  *num_lists = numMyVertices;
  *num_pins = numPins;
  *format = ZOLTAN_COMPRESSED_EDGE;
}                    

static void get_local_hypergraph(void *data, 
         int sizeGID, int num_edges, int num_pins, int format,
         ZOLTAN_ID_PTR edgeGID, int *index, ZOLTAN_ID_PTR pinGID, int *ierr)
{
  int i, next;
  ZOLTAN_ID_TYPE left, right;
  *ierr = ZOLTAN_OK;

  for (i=0,next=0; i < num_edges; i++){
    edgeGID[i] = proc_vertex_gid(myRank, i);
    index[i] = i * nborCount;

    if (i==0){
      left = proc_vertex_gid(myRank, num_edges - 1);
    }
    else{
      left = edgeGID[i] - 1;
    }
    if (i==num_edges-1){
      right = proc_vertex_gid(myRank, 0);
    }
    else{
      right = edgeGID[i] + 1;
    }

    pinGID[next++] = edgeGID[i];

    pinGID[next++] = left;

    pinGID[next++] = right;

    if (myRank > 0){
      pinGID[next++] = edgeGID[i] - num_edges;
    }
    if (myRank < numProcs-1){
      pinGID[next++] = edgeGID[i] + num_edges;
    }
  }
}

static void get_edge_weights_size(void *data, int *num_edges, int *ierr)
{
  *ierr = ZOLTAN_OK;
  *num_edges = numMyVertices;
}

static void get_edge_weights(void *data, int sizeGID, int sizeLID,
           int num_edges, int edge_weight_dim, ZOLTAN_ID_PTR edgeGID,
           ZOLTAN_ID_PTR edgeLID, float *edgeWeight, int *ierr)
{
  int i;

  for (i=0; i < numMyVertices; i++){
    edgeGID[i] = proc_vertex_gid(myRank, i);
    edgeWeight[i] = (float)edgeCutCost;
  }
}

/* Zoltan partition query functions. */

static void get_partition_list(void *data, int sizeGID, int sizeLID, int num_obj,
        ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int*parts, int *ierr)
{
  int i;
  *ierr = ZOLTAN_OK;

  if (vertex_part == NULL){
    /* All vertices are in my partition, we have not repartitioned yet */
    for (i=0; i < num_obj; i++){
      parts[i] = myRank;
    }
  }
  else{
    /* We have repartitioned the vertices */
    for (i=0; i < num_obj; i++){
      parts[i] = vertex_part[localID[i]];
    }
  }
  return;
}


int main(int argc, char *argv[])
{
  int i, rc;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  signal(SIGSEGV, meminfo_signal_handler);
  signal(SIGINT, meminfo_signal_handler);
  signal(SIGTERM, meminfo_signal_handler);
  signal(SIGABRT, meminfo_signal_handler);
  signal(SIGFPE, meminfo_signal_handler);

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
  ** Problem size
  ******************************************************************/

  numGlobalVertices = NUM_GLOBAL_VERTICES;
/*
  vwgt_dim = VERTEX_WEIGHT_DIMENSION;
  ewgt_dim = EDGE_WEIGHT_DIMENSION;
*/
  vwgt_dim = 1;
  ewgt_dim = 1;

  if (argc > 1){
    sscanf(argv[1], "%zd", &numGlobalVertices);
/*
    if (argc > 2){
      vwgt_dim = atoi(argv[2]);
      if (argc > 3){
        ewgt_dim = atoi(argv[3]);
      }
    }
*/
  }

  create_a_graph();

  /******************************************************************
  ** Create a Zoltan library structure for this instance of load
  ** balancing.  Set the parameters and query functions that will
  ** govern the library's calculation.  See the Zoltan User's
  ** Guide for the definition of these and many other parameters.
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "5");
  Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
  Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* number of weights per vertex */
  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");/* number of weights per hyperedge */

  /* PHG parameters  - see the Zoltan User's Guide for many more
   *   (The "REPARTITION" approach asks Zoltan to create a partitioning that is
   *    better but is not too far from the current partitioning, rather than partitioning 
   *    from scratch.  It may be faster but of lower quality that LB_APPROACH=PARTITION.)

  Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
  */

  /* Application defined query functions */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, NULL);
  Zoltan_Set_HG_Size_CS_Fn(zz, get_local_hypergraph_size, NULL);
  Zoltan_Set_HG_CS_Fn(zz, get_local_hypergraph, NULL);
  Zoltan_Set_HG_Size_Edge_Wts_Fn(zz, get_edge_weights_size, NULL);
  Zoltan_Set_HG_Edge_Wts_Fn(zz, get_edge_weights, NULL);

  Zoltan_Set_Part_Multi_Fn(zz, get_partition_list, NULL);

  /******************************************************************
  ** Zoltan can now partition the vertices of hypergraph.
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
  ** Check the balance of the partitions before running zoltan.
  ** The query function get_partition_list() will give the
  ** partitions of the vertices before we called Zoltan.
  ******************************************************************/

  if (myRank == 0){
    printf("\nBALANCE before running Zoltan\n");
  }

  rc = Zoltan_LB_Eval_HG(zz, 1, NULL);

  if (rc != ZOLTAN_OK){
    printf("sorry first LB_Eval_HG...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  /******************************************************************
  ** Print out the balance of the new partitions.
  ******************************************************************/

  vertex_part = (int *)malloc(sizeof(int) * numMyVertices);

  if (!vertex_part){
    printf("sorry memory error...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  for (i=0; i < numMyVertices; i++){
    vertex_part[i] = myRank;
  }

  if (numExport > 0){
    for (i=0; i < numExport; i++){
      vertex_part[exportLocalGids[i]] = exportToPart[i];
    }
  }

  if (myRank == 0){
    printf("\nBALANCE after running Zoltan\n");
  }

  rc = Zoltan_LB_Eval_HG(zz, 1, NULL);

  if (rc != ZOLTAN_OK){
    printf("sorry second LB_Eval_HG...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }


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

  if (vertex_part) free(vertex_part);
#ifndef UNIT_WEIGHTS
  if (vwgts) free(vwgts);
#endif

  return 0;
}

/* Create a simple graph.  The vertices of the graph are the objects for Zoltan to partition.
 * The graph itself can be used to generate a hypergraph in this way: A vertex and all of its
 * neighbors represent a hypergraph.
 */
static int create_a_graph()
{
  int i, nvtxs, num4, nborCount, mid;
  int nvtx;
#ifndef UNIT_WEIGHTS
  int random_weights = (vwgt_dim > 0);
  int heavyPart = (myRank % 3 == 0);
#endif

  nvtxs = (int)(numGlobalVertices / numProcs);

  if (nvtxs > 4){ 
    num4 = nvtxs / 4;
    nvtxs = num4 * 4;
  }
  else{
    num4 = 1;
    nvtxs = 4;
  }

  numGlobalVertices = (ZOLTAN_GNO_TYPE)nvtxs * numProcs;
  numMyVertices = nvtxs;

  if (numProcs == 1){
    nborCount = 3;
  }
  else if ((myRank == 0) || (myRank == (numProcs-1))){
    nborCount = 4;
  }
  else{
    nborCount = 5;
  }

  numPins = numMyVertices * nborCount;

  vertexGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * (numProcs + 1));
  vertexGIDs[0] = 0;

  for (i=1; i <= numProcs; i++){
    vertexGIDs[i] = vertexGIDs[i-1] + nvtxs;
  }

  if (myRank == 0){
    printf("Hypergraph will have %zd hyperedges, %d on each process\n", numGlobalVertices, nvtxs);
  }

  if (vwgt_dim == 0) 
    vwgt_dim = 1;

#ifndef UNIT_WEIGHTS
  vwgts = (float *)calloc( vwgt_dim * nvtxs, sizeof(float));
  if (!vwgts) return 1;
  
  srand(0);

  for (i = 0; i < nvtxs; i++)  {
    if (!random_weights){ /* Unit weights if no weights were requested. */
      vwgts[i] = 1.0;
    }
    else{
      vwgts[i*vwgt_dim] = ((float) rand())/RAND_MAX;
      if (heavyPart){
        vwgts[i*vwgt_dim] += .2;
      }
    }
  }
#endif

  /* Make edge cut costs higher in the "center" */

  mid = (int)(numProcs / 2.0);
  edgeCutCost = ((myRank < mid) ? myRank + 1 : numProcs - myRank);

  return 0;
}
