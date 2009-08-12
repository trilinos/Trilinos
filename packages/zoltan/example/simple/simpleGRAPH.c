/**************************************************************
*  Basic example of using Zoltan to partition a graph.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

/* Name of file containing graph to be partitioned */

static char *fname="graph.txt";

/* Structure to hold graph data */

typedef struct{
  int numGlobalVertices; /* total vertices in graph */
  int numMyVertices;     /* total in my partition */
  int numAllNbors;  /* total number of neighbors of my vertices */
  int *vertexGID;  /* global ID of each of my vertices */
  int *nborIndex;  /* nborIndex[i] is location of start of neighbors for vertex i */
  int *nborGID;    /* nborGIDs[nborIndex[i]] is first neighbor of vertex i */
  int *nborProc;   /* process owning each nbor in nborGID */
} GRAPH_DATA;

/* Application defined query functions */

static int get_number_of_vertices(void *data, int *ierr);
static void get_vertex_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr);
static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nborGID, int *nborProc,
        int wgt_dim, float *ewgts, int *ierr);

/* Functions to read graph in from file, distribute it, view it, handle errors */

static int get_next_line(FILE *fp, char *buf, int bufsize);
static void input_file_error(int numProcs, int tag, int startProc);
static void showSimpleGraphPartitions(int myProc, int numIDs, int *GIDs, int *parts);
static void read_input_file(int myRank, int numProcs, char *fname, GRAPH_DATA *myData);

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
  FILE *fp;
  GRAPH_DATA myGraph;

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
  ** Read graph from input file and distribute it 
  ******************************************************************/

  fp = fopen(fname, "r");
  if (!fp){
    if (myRank == 0) fprintf(stderr,"ERROR: Can not open %s\n",fname);
    MPI_Finalize();
    exit(1);
  }
  fclose(fp);

  read_input_file(myRank, numProcs, fname, &myGraph);

  /******************************************************************
  ** Create a Zoltan library structure for this instance of load
  ** balancing.  Set the parameters and query functions that will
  ** govern the library's calculation.  See the Zoltan User's
  ** Guide for the definition of these and many other parameters.
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
  Zoltan_Set_Param(zz, "FINAL_OUTPUT", "0");  /* save all info needed for stats */

#ifdef USE_EDGE_WEIGHTS
  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");
#else
  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");
#endif

  /* Graph parameters */

  Zoltan_Set_Param(zz, "PHG_EDGE_WEIGHT_OPERATION", "SUM");
  Zoltan_Set_Param(zz, "CHECK_GRAPH", "2"); 
  Zoltan_Set_Param(zz, "PHG_FROM_GRAPH_METHOD", "neighbors");
  Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

  /* Query functions - defined in simpleQueries.h */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, NULL);
  Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, NULL);
  Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, NULL);

  /* get balance and cut info before partitioning */

  rc = Zoltan_LB_Eval_Graph(zz, 1, &evalInfo);

  /******************************************************************
  ** Zoltan can now partition the simple graph.
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
  ** sending the objects to their new partitions.  Your query 
  ** functions would need to reflect the new partitioning.
  ******************************************************************/

  /******************************************************************
  ** Visualize the partitions just created.
  ** Create a list of GIDs now assigned to my partition, let
  ** process zero display the partitioning.
  ******************************************************************/
  

  ngids = get_number_of_vertices(NULL, &rc);
  gid_flags = (int *)calloc(sizeof(int) , simpleNumVertices);
  gid_list = (int *)malloc(sizeof(int) * ngids);
  lid_list = (int *)malloc(sizeof(int) * ngids);
  wgt_list = (float *)malloc(sizeof(float) * simpleNumVertices);
  get_vertex_list(NULL, 1, 1,
                  (ZOLTAN_ID_PTR)gid_list, (ZOLTAN_ID_PTR)lid_list,
                  1, wgt_list, &rc);

  draw_partitions("initial distribution", ngids, gid_list, 1, wgt_list, 1);

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
  draw_partitions("new partitioning", nextIdx, gid_flags, 1, wgt_list, 1);

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

  if (myGraph.numMyVertices > 0){
    free(myGraph.vertexGIDs);
    free(myGraph.nborIndex);
    free(myGraph.nborGIDs);
  }

  return 0;
}

/* Application defined query functions */

static int get_number_of_vertices(void *data, int *ierr)
{
  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;
  return graph->numMyVertices;
}

static void get_vertex_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
int i;

  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  /* In this example, return the IDs of our vertices, but no weights.
   * Zoltan will assume equally weighted vertices.
   */

  for (i=0; i<graph->numMyVertices; i++){
    globalID[i] = graph->vertexGID[i];
    localID[i] = i;
  }
}

static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr)
{
int i, idx;

  GRAPH_DATA *graph = (GRAPH_DATA *)data;

  if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->numMyVertices)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0;  i < num_obj ; i++){
    idx = localID[i];
    numEdges[i] = num = graph->nborIndex[idx+1] - graph->nborIndex[idx];
  }

  *ierr = ZOLTAN_OK;
  return;
}

static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nborGID, int *nborProc,
        int wgt_dim, float *ewgts, int *ierr)
{
int i, j, idx, num;
int *nextNbor, *nextProc;

  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  if ( (sizeGID != 1) || (sizeLID != 1) || 
       (num_obj != graph->numMyVertices)||
       (wgt_dim != 0)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  nextNbor = (int *)nborGID;
  nextProc = nborProc;

  for (i=0; i < num_obj; i++){

    /*
     * In this example, we are not setting edge weights.  Zoltan will
     * set each edge to weight 1.0.
     */

    idx = localID[i];
    num = graph->nborIndex[idx+1] - graph->nborIndex[idx];
    if (num != num_edges[i]){
      *ierr = ZOLTAN_FATAL;
      return;
    }

    idx = graph->nborIndex[idx];

    for (j=0; j < num; j++){

      *nextNbor++ = graph->nborGID[idx];
      *nextProc++ = graph->nborProc[idx];
      idx++;
    }
  }
  return;
}

/* Function to find next line of information in input file */
 
static int get_next_line(FILE *fp, char *buf, int bufsize)
{
int i, cval, len;
char *c;

  while (1){

    c = fgets(buf, bufsize, fp);

    if (c == NULL)
      return 0;  /* end of file */

    len = strlen(c);

    for (i=0, c=buf; i < len; i++, c++){
      cval = (int)*c; 
      if (isspace(cval) == 0) break;
    }
    if (i == len) continue;   /* blank line */
    if (*c == '#') continue;  /* comment */

    if (c != buf){
      strcpy(buf, c);
    }
    break;
  }

  return strlen(buf);  /* number of characters */
}

/* Proc 0 notifies others of error and exits */

static void input_file_error(int numProcs, int tag, int startProc)
{
int i, val;

  val = -1;

  fprintf(stderr,"ERROR in input file.\n");

  for (i=startProc; i < numProcs; i++){
    /* these procs have posted receive for "tag" */
    MPI_Send(&val, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
  }
  for (i=1; i < startProc; i++){
    /* these procs are done */
    MPI_Send(&val, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  exit(1);
}

/* Draw the partition assignments of the objects */

void showSimpleMeshPartitions(int myProc, int numIDs, int *GIDs, int *parts)
{
int partAssign[25], allPartAssign[25];
int i, j, part;

  memset(partAssign, 0, sizeof(int) * 25);

  for (i=0; i < numIDs; i++){
    partAssign[GIDs[i]-1] = parts[i];
  }

  MPI_Reduce(partAssign, allPartAssign, 25, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myProc == 0){

    for (i=20; i >= 0; i-=5){
      for (j=0; j < 5; j++){
        part = allPartAssign[i + j];
        if (j < 4)
          printf("%d-----",part);
        else
          printf("%d\n",part);
      }
      if (i > 0)
        printf("|     |     |     |     |\n");
    }
    printf("\n");
  }
}

typedef struct{
  int numGlobalVertices; /* total vertices in graph */
  int numMyVertices;     /* total in my partition */
  int numAllNbors;  /* total number of neighbors of my vertices */
  int *vertexGID;  /* global ID of each of my vertices */
  int *nborIndex;  /* nborIndex[i] is location of start of neighbors for vertex i */
  int *nborGID;    /* nborGIDs[nborIndex[i]] is first neighbor of vertex i */
  int *nborProc;   /* process owning each nbor in nborGID */
} GRAPH_DATA;

void read_input_file(int myRank, int numProcs, char *fname, GRAPH_DATA *graph)
{
char *buf;
int bufsize = 512;
int num, nobj, remaining, ack=0;
int i, j;
int *gids, *procs, *nbors, *idx;
FILE *fp;
MPI_Status status;
int ack_tag = 5, count_tag = 10, id_tag = 15;

  if (myRank == 0){

    buf = (char *)malloc(sizeof(char) * bufsize);
    fp = fopen(fname, "r");

    num = get_next_line(fp, buf, bufsize);
    if (num == 0) input_file_error(numProcs, count_tag, 1);
    num = sscanf(buf, "%d", &graph->numGlobalVertices);
    if (num != 1) input_file_error(numProcs, count_tag, 1);

    if (numProcs > 1){
      nobj = graph->numGlobalObjects / 2;
      remaining = graph->numGlobalObjects - nobj;
    }
    else{
      nobj = graph->numGlobalObjects;
      remaining = 0;
    }

    graph->myGlobalIDs = (int *)malloc(sizeof(int) * nobj);
    graph->numMyObjects = nobj;

    for (i=0; i < nobj; i++){

      num = get_next_line(fp, buf, bufsize);
      if (num == 0) input_file_error(numProcs, count_tag, 1);
      num = sscanf(buf, "%d", graph->myGlobalIDs + i);
      if (num != 1) input_file_error(numProcs, count_tag, 1);
  
    }

    gids = (int *)malloc(sizeof(int) * (nobj + 1));

    for (i=1; i < numProcs; i++){
    
      if (remaining > 1){
        nobj = remaining / 2;
        remaining -= nobj;
      }
      else if (remaining == 1){
        nobj = 1;
        remaining = 0;
      }
      else{
        nobj = 0;
      }

      if ((i == numProcs - 1) && (remaining > 0))
        nobj += remaining;

      if (nobj > 0){
        for (j=0; j < nobj; j++){
          num = get_next_line(fp, buf, bufsize);
          if (num == 0) input_file_error(numProcs, count_tag, i);
          num = sscanf(buf, "%d", gids + j);
          if (num != 1) input_file_error(numProcs, count_tag, i);
        }
      }

      MPI_Send(&nobj, 1, MPI_INT, i, count_tag, MPI_COMM_WORLD);
      MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

      if (nobj > 0)
        MPI_Send(gids, nobj, MPI_INT, i, id_tag, MPI_COMM_WORLD);
      
    }

    free(gids);
    fclose(fp);
    free(buf);

    /* signal all procs it is OK to go on */
    ack = 0;
    for (i=1; i < numProcs; i++){
      MPI_Send(&ack, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
  else{

    MPI_Recv(&graph->numMyObjects, 1, MPI_INT, 0, count_tag, MPI_COMM_WORLD, &status);
    ack = 0;
    if (graph->numMyObjects > 0){
      graph->myGlobalIDs = (int *)malloc(sizeof(int) * graph->numMyObjects);
      MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);
      MPI_Recv(graph->myGlobalIDs, graph->numMyObjects, MPI_INT, 0, 
               id_tag, MPI_COMM_WORLD, &status);
    }
    else if (graph->numMyObjects == 0){
      MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);
    }
    else{
      MPI_Finalize();
      exit(1);
    }

    MPI_Recv(&ack, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    if (ack < 0){
      MPI_Finalize();
      exit(1);
    }
  }
}
