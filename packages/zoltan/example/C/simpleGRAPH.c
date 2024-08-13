// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/**************************************************************
*  Basic example of using Zoltan to partition a graph.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

/* Name of file containing graph to be partitioned */

static char *global_fname="graph.txt";

/* Structure to hold graph data 
   ZOLTAN_ID_TYPE is defined when Zoltan is compiled.  It's size can
   be obtained at runtime by a library call.  (See zoltan_types.h).
*/

typedef struct{
  int numMyVertices; /* total vertices in in my partition */
  int numAllNbors;   /* total number of neighbors of my vertices */
  ZOLTAN_ID_TYPE *vertexGID;    /* global ID of each of my vertices */
  int *nborIndex;    /* nborIndex[i] is location of start of neighbors for vertex i */
  ZOLTAN_ID_TYPE *nborGID;      /* nborGIDs[nborIndex[i]] is first neighbor of vertex i */
  int *nborProc;     /* process owning each nbor in nborGID */
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
static int get_line_ints(char *buf, int bufsize, int *vals);
static void input_file_error(int numProcs, int tag, int startProc);
static void showGraphPartitions(int myProc, int numIDs, ZOLTAN_ID_TYPE *GIDs, int *parts, int nparts);
static void read_input_file(int myRank, int numProcs, char *fname, GRAPH_DATA *myData);
static unsigned int simple_hash(unsigned int *key, unsigned int n);


int main(int argc, char *argv[])
{
  int i, rc;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  int myRank, numProcs;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  int *parts = NULL;
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

  fp = fopen(global_fname, "r");
  if (!fp){
    if (myRank == 0) fprintf(stderr,"ERROR: Can not open %s\n",global_fname);
    MPI_Finalize();
    exit(1);
  }
  fclose(fp);

  read_input_file(myRank, numProcs, global_fname, &myGraph);

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
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* Graph parameters */

  Zoltan_Set_Param(zz, "CHECK_GRAPH", "2"); 
  Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

  /* Query functions - defined in simpleQueries.h */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &myGraph);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &myGraph);
  Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &myGraph);
  Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &myGraph);

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
  ** Visualize the graph partitioning before and after calling Zoltan.
  ******************************************************************/

  parts = (int *)malloc(sizeof(int) * myGraph.numMyVertices);

  for (i=0; i < myGraph.numMyVertices; i++){
    parts[i] = myRank;
  }

  if (myRank== 0){
    printf("\nGraph partition before calling Zoltan\n");
  }

  showGraphPartitions(myRank, myGraph.numMyVertices, myGraph.vertexGID, parts, numProcs);

  for (i=0; i < numExport; i++){
    parts[exportLocalGids[i]] = exportToPart[i];
  }

  if (myRank == 0){
    printf("Graph partition after calling Zoltan\n");
  }

  showGraphPartitions(myRank, myGraph.numMyVertices, myGraph.vertexGID, parts, numProcs);

  if (parts) free(parts);

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
    free(myGraph.vertexGID);
    free(myGraph.nborIndex);
    if (myGraph.numAllNbors > 0){
      free(myGraph.nborGID);
      free(myGraph.nborProc);
    }
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
    numEdges[i] = graph->nborIndex[idx+1] - graph->nborIndex[idx];
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
int i, j, from, to;
int *nextProc;
ZOLTAN_ID_TYPE *nextNbor;

  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  if ( (sizeGID != 1) || (sizeLID != 1) || 
       (num_obj != graph->numMyVertices)||
       (wgt_dim != 0)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  nextNbor = nborGID;
  nextProc = nborProc;

  for (i=0; i < num_obj; i++){

    /*
     * In this example, we are not setting edge weights.  Zoltan will
     * set each edge to weight 1.0.
     */

    to = graph->nborIndex[localID[i]+1];
    from = graph->nborIndex[localID[i]];
    if ((to - from) != num_edges[i]){
      *ierr = ZOLTAN_FATAL;
      return;
    }

    for (j=from; j < to; j++){

      *nextNbor++ = graph->nborGID[j];
      *nextProc++ = graph->nborProc[j];
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

/* Function to return the list of non-negative integers in a line */

static int get_line_ints(char *buf, int bufsize, int *vals)
{
char *c = buf;
int count=0;

  while (1){
    if ( (c-buf) >= bufsize) break;

    while (!(isdigit(*c))){
      if ((c - buf) >= bufsize) break;
      c++;
    }
  
    if ( (c-buf) >= bufsize) break;
  
    vals[count++] = atoi(c);
  
    while (isdigit(*c)){
      if ((c - buf) >= bufsize) break;
      c++;
    }
  
    if ( (c-buf) >= bufsize) break;
  }

  return count;
}


/* Proc 0 notifies others of error and exits */

static void input_file_error(int numProcs, int tag, int startProc)
{
int i, val[2];

  val[0] = -1;   /* error flag */

  fprintf(stderr,"ERROR in input file.\n");

  for (i=startProc; i < numProcs; i++){
    /* these procs have posted a receive for "tag" expecting counts */
    MPI_Send(val, 2, MPI_INT, i, tag, MPI_COMM_WORLD);
  }
  for (i=1; i < startProc; i++){
    /* these procs are done and waiting for ok-to-go */
    MPI_Send(val, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  exit(1);
}

/* Draw the partition assignments of the objects */

static void showGraphPartitions(int myProc, int numIDs, ZOLTAN_ID_TYPE *GIDs, int *parts, int nparts)
{
int partAssign[25], allPartAssign[25];
int i, j, part, cuts, prevPart=-1;
float imbal, localImbal, sum;
int *partCount;


  memset(partAssign, 0, sizeof(int) * 25);

  for (i=0; i < numIDs; i++){
    partAssign[GIDs[i]-1] = parts[i];
  }

  MPI_Reduce(partAssign, allPartAssign, 25, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myProc == 0){
    partCount = (int *)calloc(sizeof(int), nparts);

    cuts = 0;

    for (i=20; i >= 0; i-=5){
      for (j=0; j < 5; j++){
        part = allPartAssign[i + j];
        partCount[part]++;
        if (j > 0){
          if (part == prevPart){
            printf("-----%d",part);
          }
          else{
            printf("--x--%d",part);
            cuts++;
            prevPart = part;
          }
        }
        else{
          printf("%d",part);
          prevPart = part;
        }
      }
      printf("\n");
      if (i > 0){
        for (j=0; j < 5; j++){
          if (allPartAssign[i+j] != allPartAssign[i+j-5]){
            printf("x     ");
            cuts++;
          }
          else{
            printf("|     ");
          }
        }
        printf("\n");
      }
    }
    printf("\n");

    for (sum=0, i=0; i < nparts; i++){
      sum += partCount[i];
    }
    imbal = 0;
    for (i=0; i < nparts; i++){
      /* An imbalance measure.  1.0 is perfect balance, larger is worse */
      localImbal = (nparts * partCount[i]) / sum;
      if (localImbal > imbal) imbal = localImbal;
    }

    printf("Object imbalance (1.0 perfect, larger numbers are worse): %f\n",imbal);
    printf("Total number of edge cuts: %d\n\n",cuts);

    if (nparts) free(partCount);
  }

}

/*
 * Read the graph in the input file and distribute the vertices.
 */

void read_input_file(int myRank, int numProcs, char *fname, GRAPH_DATA *graph)
{
char buf[512];
int bufsize;
int numGlobalVertices, numGlobalNeighbors;
int num, nnbors, ack=0;
int vGID;
int i, j, procID;
int vals[128], send_count[2];
int *idx;
unsigned int id;
FILE *fp;
MPI_Status status;
int ack_tag = 5, count_tag = 10, id_tag = 15;
GRAPH_DATA *send_graph;

  if (myRank == 0){

    bufsize = 512;

    fp = fopen(fname, "r");

    /* Get the number of vertices */

    num = get_next_line(fp, buf, bufsize);
    if (num == 0) input_file_error(numProcs, count_tag, 1);
    num = sscanf(buf, "%d", &numGlobalVertices);
    if (num != 1) input_file_error(numProcs, count_tag, 1);

    /* Get the number of vertex neighbors  */

    num = get_next_line(fp, buf, bufsize);
    if (num == 0) input_file_error(numProcs, count_tag, 1);
    num = sscanf(buf, "%d", &numGlobalNeighbors);
    if (num != 1) input_file_error(numProcs, count_tag, 1);

    /* Allocate arrays to read in entire graph */

    graph->vertexGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numGlobalVertices);
    graph->nborIndex = (int *)malloc(sizeof(int) * (numGlobalVertices + 1));
    graph->nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numGlobalNeighbors);
    graph->nborProc = (int *)malloc(sizeof(int) * numGlobalNeighbors);

    graph->nborIndex[0] = 0;

    for (i=0; i < numGlobalVertices; i++){

      num = get_next_line(fp, buf, bufsize);
      if (num == 0) input_file_error(numProcs, count_tag, 1);

      num = get_line_ints(buf, strlen(buf), vals);

      if (num < 2) input_file_error(numProcs, count_tag, 1);

      vGID = vals[0];
      nnbors = vals[1];

      if (num < (nnbors + 2)) input_file_error(numProcs, count_tag, 1);

      graph->vertexGID[i] = (ZOLTAN_ID_TYPE)vGID;

      for (j=0; j < nnbors; j++){
        graph->nborGID[graph->nborIndex[i] + j] = (ZOLTAN_ID_TYPE)vals[2 + j];
      }

      graph->nborIndex[i+1] = graph->nborIndex[i] + nnbors;
    }

    fclose(fp);

    /* Assign each vertex to a process using a hash function */

    for (i=0; i <numGlobalNeighbors; i++){
      id = (unsigned int)graph->nborGID[i];
      graph->nborProc[i] = simple_hash(&id, numProcs);
    } 

    /* Create a sub graph for each process */

    send_graph = (GRAPH_DATA *)calloc(sizeof(GRAPH_DATA) , numProcs);

    for (i=0; i < numGlobalVertices; i++){
      id = (unsigned int)graph->vertexGID[i];
      procID = simple_hash(&id, numProcs);
      send_graph[procID].numMyVertices++;
    }

    for (i=0; i < numProcs; i++){
      num = send_graph[i].numMyVertices;
      send_graph[i].vertexGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * num);
      send_graph[i].nborIndex = (int *)calloc(sizeof(int) , (num + 1));
    }

    idx = (int *)calloc(sizeof(int), numProcs);

    for (i=0; i < numGlobalVertices; i++){

      id = (unsigned int)graph->vertexGID[i];
      nnbors = graph->nborIndex[i+1] - graph->nborIndex[i];
      procID = simple_hash(&id, numProcs);

      j = idx[procID];
      send_graph[procID].vertexGID[j] = (ZOLTAN_ID_TYPE)id;
      send_graph[procID].nborIndex[j+1] = send_graph[procID].nborIndex[j] + nnbors;

      idx[procID] = j+1;
    }

    for (i=0; i < numProcs; i++){

      num = send_graph[i].nborIndex[send_graph[i].numMyVertices];

      send_graph[i].nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * num);
      send_graph[i].nborProc= (int *)malloc(sizeof(int) * num);

      send_graph[i].numAllNbors = num;
    }

    memset(idx, 0, sizeof(int) * numProcs);

    for (i=0; i < numGlobalVertices; i++){

      id = (unsigned int)graph->vertexGID[i];
      nnbors = graph->nborIndex[i+1] - graph->nborIndex[i];
      procID = simple_hash(&id, numProcs);
      j = idx[procID];

      if (nnbors > 0){
        memcpy(send_graph[procID].nborGID + j, graph->nborGID + graph->nborIndex[i],
               nnbors * sizeof(ZOLTAN_ID_TYPE));
  
        memcpy(send_graph[procID].nborProc + j, graph->nborProc + graph->nborIndex[i],
               nnbors * sizeof(int));
  
        idx[procID] = j + nnbors;
      }
    }

    free(idx);

    /* Process zero sub-graph */

    free(graph->vertexGID);
    free(graph->nborIndex);
    free(graph->nborGID);
    free(graph->nborProc);

    *graph = send_graph[0];

    /* Send other processes their subgraph */

    for (i=1; i < numProcs; i++){
      send_count[0] = send_graph[i].numMyVertices;
      send_count[1] = send_graph[i].numAllNbors;

      MPI_Send(send_count, 2, MPI_INT, i, count_tag, MPI_COMM_WORLD);
      MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

      if (send_count[0] > 0){

        MPI_Send(send_graph[i].vertexGID, send_count[0], ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
        free(send_graph[i].vertexGID);

        MPI_Send(send_graph[i].nborIndex, send_count[0] + 1, MPI_INT, i, id_tag + 1, MPI_COMM_WORLD);
        free(send_graph[i].nborIndex);

        if (send_count[1] > 0){
          MPI_Send(send_graph[i].nborGID, send_count[1], ZOLTAN_ID_MPI_TYPE, i, id_tag + 2, MPI_COMM_WORLD);
          free(send_graph[i].nborGID);

          MPI_Send(send_graph[i].nborProc, send_count[1], MPI_INT, i, id_tag + 3, MPI_COMM_WORLD);
          free(send_graph[i].nborProc);
        }
      }
    }

    free(send_graph);

    /* signal all procs it is OK to go on */
    ack = 0;
    for (i=1; i < numProcs; i++){
      MPI_Send(&ack, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
  else{

    MPI_Recv(send_count, 2, MPI_INT, 0, count_tag, MPI_COMM_WORLD, &status);

    if (send_count[0] < 0){
      MPI_Finalize();
      exit(1);
    }

    ack = 0;

    graph->numMyVertices = send_count[0];
    graph->numAllNbors  = send_count[1];

    if (send_count[0] > 0){

      graph->vertexGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
      graph->nborIndex = (int *)malloc(sizeof(int) * (send_count[0] + 1));

      if (send_count[1] > 0){
        graph->nborGID  = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
        graph->nborProc = (int *)malloc(sizeof(int) * send_count[1]);
      }
    }

    MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);

    if (send_count[0] > 0){
      MPI_Recv(graph->vertexGID,send_count[0],ZOLTAN_ID_MPI_TYPE, 0, id_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(graph->nborIndex,send_count[0] + 1, MPI_INT, 0, id_tag + 1, MPI_COMM_WORLD, &status);

      if (send_count[1] > 0){
        MPI_Recv(graph->nborGID,send_count[1], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 2, MPI_COMM_WORLD, &status);
        MPI_Recv(graph->nborProc,send_count[1], MPI_INT, 0, id_tag + 3, MPI_COMM_WORLD, &status);
      }
    }

    /* ok to go on? */

    MPI_Recv(&ack, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    if (ack < 0){
      MPI_Finalize();
      exit(1);
    }
  }
}

unsigned int simple_hash(unsigned int *key, unsigned int n)
{
  unsigned int h, rest, *p, bytes, num_bytes;
  char *byteptr;

  num_bytes = (unsigned int) sizeof(int);

  /* First hash the int-sized portions of the key */
  h = 0;
  for (p = (unsigned int *)key, bytes=num_bytes;
       bytes >= (unsigned int) sizeof(int);
       bytes-=sizeof(int), p++){
    h = (h*2654435761U) ^ (*p);
  }

  /* Then take care of the remaining bytes, if any */
  rest = 0;
  for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++){
    rest = (rest<<8) | (*byteptr);
  }

  /* Merge the two parts */
  if (rest)
    h = (h*2654435761U) ^ rest;

  /* Return h mod n */
  return (h%n);
}

