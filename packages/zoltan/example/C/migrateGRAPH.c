// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/**************************************************************
*  An expansion of simpleGRAPH.c.
*
*  We use Zoltan's distributed directory utility to create a
*  global directory of object (graph vertex) locations.  We use
*  Zoltan's migration capabilities to move the graph vertices.
*
*  In this example, part numbers and process ranks are the same.
*  In general, the number of parts need not equal the number
*  of processes.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

/***************************************************************************/
/* Name of file containing graph to be partitioned */

static char *global_fname="graph.txt";

/***************************************************************************/
/* Structure to hold graph data */
typedef struct{
  int numMyVertices;         /* total vertices in my part */
  ZOLTAN_ID_TYPE *vertexGID; /* global ID of each of my vertices */
  int *nborIndex;            /* index into nbor info;
                                last element is total nbors */
  ZOLTAN_ID_TYPE *nborGID;   /* nborGIDs[nborIndex[i]] is first nbor of vtx i */
  int *nborPart;             /* process owning each nbor in nborGID */
  int vertex_capacity;       /* size of vertexGID buffer */
  int nbor_capacity;         /* size of nborGID & nborPart buffers */
  struct Zoltan_DD_Struct *dd;
} GRAPH_DATA;

/***************************************************************************/
/* Application-defined query functions used in partitioning ****************/
/***************************************************************************/

static int get_number_of_vertices(void *data, int *ierr)
{
  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;
  return graph->numMyVertices;
}

/***************************************************************************/
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
    localID[i] = (ZOLTAN_ID_TYPE)i;
  }
}

/***************************************************************************/
static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr)
{
int i;
ZOLTAN_ID_TYPE idx;

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

/***************************************************************************/
static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nborGID, int *nborPart,
        int wgt_dim, float *ewgts, int *ierr)
{
int i, j, from, to;
ZOLTAN_ID_TYPE *nextNbor;
int *nextProc;

  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  if ( (sizeGID != 1) || (sizeLID != 1) ||
       (num_obj != graph->numMyVertices)||
       (wgt_dim != 0)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  nextNbor = nborGID;
  nextProc = nborPart;

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
      *nextProc++ = graph->nborPart[j];
    }
  }
  return;
}

/***************************************************************************/
/* Application-defined query functions used in migrating *******************
 *
 * Migration strategy:
 *   The data sent for each vertex is the list 
 *   of the vertex neighbor global IDs.
 *   At the pack query function, copy the vertex
 *   data to the Zoltan-supplied buffer.
 *   At the mid-migrate query function, rewrite 
 *   the local graph data omitting the
 *   exported vertices.
 *   At the unpack query function, copy the 
 *   new vertices to the local graph data.
 */
/***************************************************************************/

static void get_message_sizes(void *data, int gidSize, int lidSize, int num_ids,
                              ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                              int *sizes, int *ierr)
{
  GRAPH_DATA *graph;
  int i, len;

  graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  for (i=0; i < num_ids; i++){
    len = graph->nborIndex[localID[i]+1] - graph->nborIndex[localID[i]];
    sizes[i] = (len+1) * sizeof(ZOLTAN_ID_TYPE);
  }
}

/***************************************************************************/
static void pack_object_messages(void *data, int gidSize, int lidSize,
                                 int num_ids, ZOLTAN_ID_PTR globalIDs,
                                 ZOLTAN_ID_PTR localIDs, int *dests, int *sizes,
                                 int *idx, char *buf, int *ierr)
{
  int i, j, num_nbors;
  ZOLTAN_ID_TYPE *nbors=NULL, *ibuf=NULL;
  ZOLTAN_ID_TYPE lid;
  GRAPH_DATA *graph;
  *ierr = ZOLTAN_OK;
  graph = (GRAPH_DATA *)data;

  /* For each exported vertex, write its neighbor 
   * global IDs to the supplied buffer 
   */

  for (i=0; i < num_ids; i++){
    lid = localIDs[i];
    nbors = graph->nborGID + graph->nborIndex[lid];
    num_nbors = graph->nborIndex[lid+1] - graph->nborIndex[lid];

    ibuf = (ZOLTAN_ID_TYPE *)(buf + idx[i]);
    ibuf[0] = num_nbors;

    for (j=1; j <= num_nbors; j++){
      ibuf[j] = *nbors++;
    }
  }
}

/***************************************************************************/
static void mid_migrate(void *data, int gidSize, int lidSize,
                        int numImport, ZOLTAN_ID_PTR importGlobalID,
                        ZOLTAN_ID_PTR importLocalID,
                        int *importProc, int *importPart,
                        int numExport, ZOLTAN_ID_PTR exportGlobalID,
                        ZOLTAN_ID_PTR exportLocalID,
                        int *exportProc, int *exportPart, int *ierr)
{
  GRAPH_DATA *graph;
  int i, len, next_vertex, next_nbor;
  int *exports;

  *ierr = ZOLTAN_OK;
  graph = (GRAPH_DATA *)data;

  /* The exported vertices have been packed.  Remove them from local graph. */

  exports = (int *)calloc(sizeof(int) , graph->numMyVertices);
  for (i=0; i <numExport; i++){
    exports[exportLocalID[i]] = 1;
  }

  next_vertex = 0;
  next_nbor = 0;

  graph->nborIndex[0] = 0;

  for (i=0; i < graph->numMyVertices; i++){
    if (exports[i] == 0){
      len = graph->nborIndex[i+1] - graph->nborIndex[i];

      if (i > next_vertex){
        graph->vertexGID[next_vertex] = graph->vertexGID[i];
        if (len > 0){
          memcpy(graph->nborGID + next_nbor,
                 graph->nborGID + graph->nborIndex[i],
                 sizeof(ZOLTAN_ID_TYPE) * len);
          memcpy(graph->nborPart + next_nbor,
                 graph->nborPart + graph->nborIndex[i], sizeof(int) * len);
        }
        graph->nborIndex[next_vertex+1] = graph->nborIndex[next_vertex] + len;
      }
      next_nbor += len;
      next_vertex++;
    }
  }

  free(exports);
  graph->numMyVertices = next_vertex;
}

/***************************************************************************/
static void unpack_object_messages(void *data, int gidSize, int num_ids,
                     ZOLTAN_ID_PTR globalIDs,
                     int *size, int *idx, char *buf, int *ierr)
{
  int i, len, num_nbors, num_vertex, next_vertex, next_nbor;
  ZOLTAN_ID_TYPE *ibuf=NULL;
  GRAPH_DATA *graph;
  *ierr = ZOLTAN_OK;
  graph = (GRAPH_DATA *)data;

  /* Add incoming vertices to local graph */

  next_vertex = graph->numMyVertices;
  next_nbor = graph->nborIndex[next_vertex];

  num_nbors = 0;
  for (i=0; i < num_ids; i++){
    num_nbors += (size[i]);
  }
  num_nbors /= sizeof(ZOLTAN_ID_TYPE); /* number of incoming neighbors */
  num_nbors -= num_ids;

  num_nbors += next_nbor;              /* plus number of existing neighbors */

  num_vertex = next_vertex + num_ids;

  if (num_vertex > graph->vertex_capacity){
    graph->vertexGID = (ZOLTAN_ID_TYPE *)
                        realloc(graph->vertexGID,
                                sizeof(ZOLTAN_ID_TYPE) * num_vertex);
    graph->nborIndex = (int *)realloc(graph->nborIndex,
                                      sizeof(int) * (num_vertex+1));
    graph->vertex_capacity = num_vertex;
  }

  if (num_nbors > graph->nbor_capacity){
    graph->nborGID = (ZOLTAN_ID_TYPE *)
                      realloc(graph->nborGID,
                              sizeof(ZOLTAN_ID_TYPE) * num_nbors);
    graph->nborPart = (int *)realloc(graph->nborPart, sizeof(int) * num_nbors);
    graph->nbor_capacity = num_nbors;
  }

  for (i=0; i < num_ids; i++){
    graph->vertexGID[next_vertex] = globalIDs[i];
    ibuf = (ZOLTAN_ID_TYPE *)(buf + idx[i]);
    len = ibuf[0];

    if (len > 0)
      memcpy(graph->nborGID + next_nbor, &ibuf[1],
             len * sizeof(ZOLTAN_ID_TYPE));
    graph->nborIndex[next_vertex+1] = graph->nborIndex[next_vertex] + len;
    next_vertex++;
    next_nbor += len;
  }

  graph->numMyVertices += num_ids;
}

/****************************************************************************/
/* Functions to read graph from file, distribute it, view it, handle errors */
/****************************************************************************/

static int get_next_line(FILE *, char *, int);
static int get_line_ints(char *, int, int *);
static void input_file_error(int, int, int);
static void printGraph(int, GRAPH_DATA *);
static void showGraphParts(int, struct Zoltan_DD_Struct *);
static void read_input_file(int, int, char *, GRAPH_DATA *);

/****************************************************************************/
/****************************************************************************/

int main(int argc, char *argv[])
{
  int i, rc;
  int myRank, numProcs;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids;
  ZOLTAN_ID_PTR exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  int *parts=NULL;
  ZOLTAN_ID_PTR lids=NULL;
  FILE *fp;
  struct Zoltan_DD_Struct *dd;
  GRAPH_DATA myGraph;
  int gid_length = 1;   /* our global IDs consist of 1 integer */
  int lid_length = 1;   /* our local IDs consist of 1 integer */

  /******************************************************************
  ** Initialize MPI and Zoltan
  ******************************************************************/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    printf("Failed to initialize Zoltan -- sorry...\n");
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
  ** Create a distributed data directory which maps vertex
  ** global IDs to their current part number.  We'll use this
  ** after migrating vertices, to update the parts in which
  ** our vertices neighbors are.
  **
  ** Our local IDs (array "lids") are of type ZOLTAN_ID_TYPE because
  ** we are using Zoltan's distributed data directory.  It assumes
  ** that a global ID is a sequence of "gid_length" ZOLTAN_ID_TYPEs.
  ** It assumes that a local ID is a sequence of "lid_length"
  ** ZOLTAN_ID_TYPEs.
  ******************************************************************/

  rc = Zoltan_DD_Create(&dd, MPI_COMM_WORLD,
                           gid_length,    /* length of a global ID */
                           lid_length,    /* length of a local ID */
                           0,             /* length of user data  */
                           myGraph.numMyVertices,  /* hash table size */
                           0);                     /* debug level */

  parts = (int *) malloc(myGraph.numMyVertices * sizeof(int));
  lids = (ZOLTAN_ID_TYPE *) 
          malloc(myGraph.numMyVertices * sizeof(ZOLTAN_ID_TYPE));

  for (i=0; i < myGraph.numMyVertices; i++){
    parts[i] = myRank;   /* part number of this vertex */
    lids[i] = (ZOLTAN_ID_TYPE)i;  /* local ID on my process for this vertex */
  }

  rc = Zoltan_DD_Update(dd, myGraph.vertexGID, lids, NULL, parts,
                        myGraph.numMyVertices);


  myGraph.dd = dd;

  /******************************************************************
  ** Create a Zoltan structure for this instance of load
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
  Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all,
                                                              1-remove none */

  /* Query functions, defined in this source file */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &myGraph);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &myGraph);
  Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &myGraph);
  Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &myGraph);

  Zoltan_Set_Obj_Size_Multi_Fn(zz, get_message_sizes,&myGraph);
  Zoltan_Set_Pack_Obj_Multi_Fn(zz, pack_object_messages,&myGraph);
  Zoltan_Set_Unpack_Obj_Multi_Fn(zz, unpack_object_messages,&myGraph);
  Zoltan_Set_Mid_Migrate_PP_Fn(zz, mid_migrate,&myGraph);

  /******************************************************************
  ** Visualize the graph partition before calling Zoltan.
  ******************************************************************/

  if (myRank== 0){
    printf("\nGraph partition before calling Zoltan\n");
  }

  printGraph(myRank, &myGraph);
  showGraphParts(myRank, myGraph.dd);

  /******************************************************************
  ** Zoltan can now partition the simple graph.
  ** In this simple example, we assume the number of parts is
  ** equal to the number of processes.  Process rank 0 will own
  ** part 0, process rank 1 will own part 1, and so on.
  ******************************************************************/

  rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partition was changed, 0 otherwise */
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of vertices to be sent to me */
        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
        &importLocalGids,   /* Local IDs of vertices to be sent to me */
        &importProcs,    /* Process rank for source of each incoming vertex */
        &importToPart,   /* New part for each incoming vertex */
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
  ** Update the data directory with the new part numbers
  ******************************************************************/

  for (i=0; i < numExport; i++){
    parts[exportLocalGids[i]] = exportToPart[i];
  }

  rc = Zoltan_DD_Update(dd, myGraph.vertexGID, lids, NULL, parts,
                        myGraph.numMyVertices);

  /******************************************************************
  ** Migrate vertices to new parts
  ******************************************************************/

  rc = Zoltan_Migrate(zz,
                      numImport, importGlobalGids, importLocalGids,
                      importProcs, importToPart,
                      numExport, exportGlobalGids, exportLocalGids,
                      exportProcs, exportToPart);


  /******************************************************************
  ** Use the data dictionary to find neighbors' parts
  ******************************************************************/

  rc = Zoltan_DD_Find(dd,
             (ZOLTAN_ID_PTR)(myGraph.nborGID), NULL, NULL,
              myGraph.nborPart, myGraph.nborIndex[myGraph.numMyVertices], NULL);

  /******************************************************************
  ** Visualize the graph partition after calling Zoltan.
  ******************************************************************/

  if (myRank == 0){
    printf("Graph partition after calling Zoltan\n");
  }
  printGraph(myRank, &myGraph);
  showGraphParts(myRank, myGraph.dd);

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

  if (myGraph.vertex_capacity > 0){
    free(myGraph.vertexGID);
    free(myGraph.nborIndex);
    if (myGraph.nbor_capacity > 0){
      free(myGraph.nborGID);
      free(myGraph.nborPart);
    }
  }
  Zoltan_DD_Destroy(&dd);

  if (parts) free(parts);
  if (lids) free(lids);

  MPI_Finalize();
  return 0;
}

/***************************************************************************/
/* Find next line of information in input file */
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


/* Return the list of non-negative integers in a line */
static int get_line_ints(char *buf, int bufsize, int *vals)
{
char *c = buf;
int count=0;

  while (1){
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
static void input_file_error(int nproc , int tag, int startProc)
{
int i, val[2];

  val[0] = -1;   /* error flag */

  fprintf(stderr,"ERROR in input file.\n");

  for (i=startProc; i < nproc; i++){
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


/* Print the entries in the graph data structure */
static void printGraph(int me, GRAPH_DATA *graph)
{
int i, j;
  printf("%d PG nVtx %d vtxCap %d nborCap %d\n",
         me, graph->numMyVertices,
         graph->vertex_capacity, graph->nbor_capacity);

  for (i = 0; i < graph->numMyVertices; i++) {
    printf("%d PG %d with %d: ",
           me, graph->vertexGID[i], graph->nborIndex[i+1]-graph->nborIndex[i]);
    for (j = graph->nborIndex[i]; j < graph->nborIndex[i+1]; j++)
      printf("(%d,%d) ", graph->nborGID[j], graph->nborPart[j]);
    printf("\n");
    fflush(stdout);
  }
}

/* Draw part assignments of the objects.  We know globals IDs are [1,25] */
static void showGraphParts(int me, struct Zoltan_DD_Struct *dd)
{
int i, j, part, cuts, prevPart=-1, nparts=0;
int *partCount;
float imbal, localImbal, sum;
ZOLTAN_ID_TYPE gid[25];
int parts[25];

  for (i=0; i < 25; i++){
    gid[i] = i+1;
  }

  Zoltan_DD_Find(dd, gid, NULL, NULL, parts, 25, NULL);

  if (me > 0){
    return;
  }

  for (i=0; i < 25; i++){
    if (parts[i] > nparts)
      nparts = parts[i];
  }

  nparts++;

  partCount = (int *)calloc(sizeof(int), nparts);

  cuts = 0;

  for (i=20; i >= 0; i-=5){
    for (j=0; j < 5; j++){
      part = parts[i + j];
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
        if (parts[i+j] != parts[i+j-5]){
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

  printf("Imbalance (1.0 perfect, larger numbers are worse): %f\n",imbal);
  printf("Total number of edge cuts: %d\n\n",cuts);

  if (nparts) free(partCount);
}


/* Read the graph in the input file and distribute the vertices.  */
void read_input_file(int me , int nproc, char *fname, GRAPH_DATA *graph)
{
char buf[512];
int bufsize;
int numGlobalVertices, numGlobalNeighbors;
int num, nnbors, ack=0;
int vGID;
int i, j, procID;
int vals[128], send_count[2];
int *idx;
ZOLTAN_ID_TYPE id;
FILE *fp;
MPI_Status status;
int ack_tag = 5, count_tag = 10, id_tag = 15;
GRAPH_DATA *send_graph;

  if (me == 0){

    bufsize = 512;

    fp = fopen(fname, "r");

    /* Get the number of vertices */

    num = get_next_line(fp, buf, bufsize);
    if (num == 0) input_file_error(nproc, count_tag, 1);
    num = sscanf(buf, "%d", &numGlobalVertices);
    if (num != 1) input_file_error(nproc, count_tag, 1);

    /* Get the number of vertex neighbors  */

    num = get_next_line(fp, buf, bufsize);
    if (num == 0) input_file_error(nproc, count_tag, 1);
    num = sscanf(buf, "%d", &numGlobalNeighbors);
    if (num != 1) input_file_error(nproc, count_tag, 1);

    /* Allocate arrays to read in entire graph */

    graph->vertexGID = (ZOLTAN_ID_TYPE *)
                        malloc(sizeof(ZOLTAN_ID_TYPE) * numGlobalVertices);
    graph->nborIndex = (int *)malloc(sizeof(int) * (numGlobalVertices + 1));
    graph->nborGID = (ZOLTAN_ID_TYPE *)
                      malloc(sizeof(ZOLTAN_ID_TYPE) * numGlobalNeighbors);
    graph->nborPart = (int *)malloc(sizeof(int) * numGlobalNeighbors);

    graph->vertex_capacity = numGlobalVertices;
    graph->nbor_capacity = numGlobalNeighbors;

    graph->nborIndex[0] = 0;

    for (i=0; i < numGlobalVertices; i++){

      num = get_next_line(fp, buf, bufsize);
      if (num == 0) input_file_error(nproc, count_tag, 1);

      num = get_line_ints(buf, bufsize, vals);

      if (num < 2) input_file_error(nproc, count_tag, 1);

      vGID = vals[0];
      nnbors = vals[1];

      if (num < (nnbors + 2)) input_file_error(nproc, count_tag, 1);

      graph->vertexGID[i] = (ZOLTAN_ID_TYPE)vGID;

      for (j=0; j < nnbors; j++){
        graph->nborGID[graph->nborIndex[i] + j] = (ZOLTAN_ID_TYPE)vals[2 + j];
      }

      graph->nborIndex[i+1] = graph->nborIndex[i] + nnbors;
    }

    fclose(fp);

    /* Assign each vertex to a process using a hash function */

    for (i=0; i <numGlobalNeighbors; i++){
      id = graph->nborGID[i];
      graph->nborPart[i] = id % nproc;  /* Cyclic distribution */
    }

    /* Create a sub graph for each process */

    send_graph = (GRAPH_DATA *)calloc(sizeof(GRAPH_DATA) , nproc);

    for (i=0; i < numGlobalVertices; i++){
      id = graph->vertexGID[i];
      procID = id % nproc;  /* Cyclic distribution */
      send_graph[procID].numMyVertices++;
    }

    for (i=0; i < nproc; i++){
      num = send_graph[i].numMyVertices;
      send_graph[i].vertexGID = (ZOLTAN_ID_TYPE *)
                                 malloc(sizeof(ZOLTAN_ID_TYPE) * num);
      send_graph[i].nborIndex = (int *)calloc(sizeof(int) , (num + 1));
    }

    idx = (int *)calloc(sizeof(int), nproc);

    for (i=0; i < numGlobalVertices; i++){

      id = graph->vertexGID[i];
      nnbors = graph->nborIndex[i+1] - graph->nborIndex[i];
      procID = id % nproc;  /* Cyclic distribution */

      j = idx[procID];
      send_graph[procID].vertexGID[j] = id;
      send_graph[procID].nborIndex[j+1] =send_graph[procID].nborIndex[j]+nnbors;
      idx[procID] = j+1;
    }

    for (i=0; i < nproc; i++){

      num = send_graph[i].nborIndex[send_graph[i].numMyVertices];
      send_graph[i].nborGID = (ZOLTAN_ID_TYPE *)
                               malloc(sizeof(ZOLTAN_ID_TYPE) * num);
      send_graph[i].nborPart= (int *)malloc(sizeof(int) * num);
    }

    memset(idx, 0, sizeof(int) * nproc);

    for (i=0; i < numGlobalVertices; i++){
      id = graph->vertexGID[i];
      nnbors = graph->nborIndex[i+1] - graph->nborIndex[i];
      procID = id % nproc;  /* Cyclic distribution */
      j = idx[procID];

      if (nnbors > 0){
        memcpy(send_graph[procID].nborGID + j,
               graph->nborGID + graph->nborIndex[i],
               nnbors * sizeof(ZOLTAN_ID_TYPE));
        memcpy(send_graph[procID].nborPart + j,
               graph->nborPart + graph->nborIndex[i],
               nnbors * sizeof(int));
        idx[procID] = j + nnbors;
      }
    }

    free(idx);

    /* Process zero sub-graph */

    free(graph->vertexGID);
    free(graph->nborIndex);
    free(graph->nborGID);
    free(graph->nborPart);

    *graph = send_graph[0];

    /* Send other processes their subgraph */

    for (i=1; i < nproc; i++){
      send_count[0] = send_graph[i].numMyVertices;
      send_count[1] = send_graph[i].nborIndex[send_graph[i].numMyVertices];

      MPI_Send(send_count, 2, MPI_INT, i, count_tag, MPI_COMM_WORLD);
      MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

      if (send_count[0] > 0){

        MPI_Send(send_graph[i].vertexGID, send_count[0], ZOLTAN_ID_MPI_TYPE,
                 i, id_tag, MPI_COMM_WORLD);
        free(send_graph[i].vertexGID);
        MPI_Send(send_graph[i].nborIndex, send_count[0] + 1, MPI_INT,
                 i, id_tag + 1, MPI_COMM_WORLD);
        free(send_graph[i].nborIndex);

        if (send_count[1] > 0){
          MPI_Send(send_graph[i].nborGID, send_count[1], ZOLTAN_ID_MPI_TYPE,
                   i, id_tag + 2, MPI_COMM_WORLD);
          free(send_graph[i].nborGID);
          MPI_Send(send_graph[i].nborPart, send_count[1], MPI_INT,
                   i, id_tag + 3, MPI_COMM_WORLD);
          free(send_graph[i].nborPart);
        }
      }
    }

    free(send_graph);

    /* signal all procs it is OK to go on */
    ack = 0;
    for (i=1; i < nproc; i++){
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
    if (send_count[0] > 0){
      graph->vertexGID = (ZOLTAN_ID_TYPE *)
                          malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
      graph->nborIndex = (int *)malloc(sizeof(int) * (send_count[0] + 1));

      if (send_count[1] > 0){
        graph->nborGID  = (ZOLTAN_ID_TYPE *)
                           malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
        graph->nborPart = (int *)malloc(sizeof(int) * send_count[1]);
      }
    }

    MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);

    if (send_count[0] > 0){
      MPI_Recv(graph->vertexGID,send_count[0], ZOLTAN_ID_MPI_TYPE, 0,
               id_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(graph->nborIndex,send_count[0] + 1, MPI_INT, 0,
               id_tag + 1, MPI_COMM_WORLD, &status);

      if (send_count[1] > 0){
        MPI_Recv(graph->nborGID,send_count[1], ZOLTAN_ID_MPI_TYPE, 0,
                 id_tag + 2, MPI_COMM_WORLD, &status);
        MPI_Recv(graph->nborPart,send_count[1], MPI_INT, 0, 
                 id_tag + 3, MPI_COMM_WORLD, &status);
      }
    }
    /* ok to go on? */
    MPI_Recv(&ack, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    if (ack < 0){
      MPI_Finalize();
      exit(1);
    }
    graph->vertex_capacity = send_count[0];
    graph->nbor_capacity = send_count[1];
  }
}
