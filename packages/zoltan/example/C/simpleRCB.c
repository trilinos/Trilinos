// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/***************************************************************
** Basic example of using Zoltan to compute an RCB partitioning
** of a very simple mesh or graph.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"


/* Name of file containing the mesh to be partitioned */

static char *global_fname="mesh.txt";

/* Structure to hold mesh data */

typedef struct{
  int numGlobalPoints;
  int numMyPoints;
  ZOLTAN_ID_PTR myGlobalIDs;
  float *x;
  float *y;
} MESH_DATA;

/* Application defined query functions */

static int get_number_of_objects(void *data, int *ierr);
static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
static int get_num_geometry(void *data, int *ierr);
static void get_geometry_list(void *data, int sizeGID, int sizeLID,
             int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr);

/* read in and display input mesh, handle errors */

static int get_next_line(FILE *fp, char *buf, int bufsize);
static void input_file_error(int numProcs, int tag, int startProc);
void read_input_objects(int myRank, int numProcs, char *fname, MESH_DATA *myData);
void showSimpleMeshPartitions(int myProc, int numIDs, ZOLTAN_ID_PTR IDs, int *parts);

int main(int argc, char *argv[])
{
  int rc, i, myRank, numProcs;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  int *parts;
  FILE *fp;
  MESH_DATA myMesh;

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
  ** Read geometry from input file and distribute it unevenly
  ******************************************************************/

  fp = fopen(global_fname, "r");
  if (!fp){
    if (myRank == 0) fprintf(stderr,"ERROR: Can not open %s\n",global_fname);
    MPI_Finalize();
    exit(1);
  }
  fclose(fp);

  read_input_objects(myRank, numProcs, global_fname, &myMesh);

  /******************************************************************
  ** Create a Zoltan library structure for this instance of load
  ** balancing.  Set the parameters and query functions that will
  ** govern the library's calculation.  See the Zoltan User's
  ** Guide for the definition of these and many other parameters.
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "1");
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* RCB parameters */

  Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"); 
  /*Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "0"); */

  /* Query functions, to provide geometry to Zoltan */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &myMesh);
  Zoltan_Set_Obj_List_Fn(zz, get_object_list, &myMesh);
  Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &myMesh);
  Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, &myMesh);

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
  ** Visualize the mesh partitioning before and after calling Zoltan.
  ******************************************************************/

  parts = (int *)malloc(sizeof(int) * myMesh.numMyPoints);

  for (i=0; i < myMesh.numMyPoints; i++){
    parts[i] = myRank;
  }

  if (myRank== 0){
    printf("\nMesh partition assignments before calling Zoltan\n");
  }

  showSimpleMeshPartitions(myRank, myMesh.numMyPoints, myMesh.myGlobalIDs, parts);

  for (i=0; i < numExport; i++){
    parts[exportLocalGids[i]] = exportToPart[i];
  }

  if (myRank == 0){
    printf("Mesh partition assignments after calling Zoltan\n");
  }

  showSimpleMeshPartitions(myRank, myMesh.numMyPoints, myMesh.myGlobalIDs, parts);

  free(parts);

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

  if (myMesh.numMyPoints > 0){
    free(myMesh.myGlobalIDs);
    free(myMesh.x);
    free(myMesh.y);
  }

  return 0;
}

/* Application defined query functions */

static int get_number_of_objects(void *data, int *ierr)
{
  MESH_DATA *mesh= (MESH_DATA *)data;
  *ierr = ZOLTAN_OK;
  return mesh->numMyPoints;
}

static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
int i;
  MESH_DATA *mesh= (MESH_DATA *)data;
  *ierr = ZOLTAN_OK;

  /* In this example, return the IDs of our objects, but no weights.
   * Zoltan will assume equally weighted objects.
   */

  for (i=0; i<mesh->numMyPoints; i++){
    globalID[i] = mesh->myGlobalIDs[i];
    localID[i] = i;
  }
}

static int get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 2;
}

static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr)
{
int i;

  MESH_DATA *mesh= (MESH_DATA *)data;

  if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 2)){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  *ierr = ZOLTAN_OK;

  for (i=0;  i < num_obj ; i++){
    geom_vec[2*i] = (double)mesh->x[i];
    geom_vec[2*i + 1] = (double)mesh->y[i];
  }

  return;
}

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


/* Proc 0 reads the points in the input file and divides them across processes */

void read_input_objects(int myRank, int numProcs, char *fname, MESH_DATA *myMesh)
{
char *buf;
int bufsize = 512;
int num, nobj, remaining, ack=0;
int i, j;
ZOLTAN_ID_PTR gids;
float *xcoord, *ycoord;
FILE *fp;
MPI_Status status;
int ack_tag = 5, count_tag = 10, id_tag = 15;
int x_tag = 20, y_tag = 25;

  if (myRank == 0){

    buf = (char *)malloc(sizeof(char) * bufsize);
    fp = fopen(fname, "r");

    num = get_next_line(fp, buf, bufsize);
    if (num == 0) input_file_error(numProcs, count_tag, 1);
    num = sscanf(buf, "%d", &myMesh->numGlobalPoints);
    if (num != 1) input_file_error(numProcs, count_tag, 1);

    if (numProcs > 1){
      nobj = myMesh->numGlobalPoints / 2;
      remaining = myMesh->numGlobalPoints - nobj;
    }
    else{
      nobj = myMesh->numGlobalPoints;
      remaining = 0;
    }

    myMesh->myGlobalIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nobj);
    myMesh->x = (float *)malloc(sizeof(float) * nobj);
    myMesh->y = (float *)malloc(sizeof(float) * nobj);
    myMesh->numMyPoints= nobj;

    for (i=0; i < nobj; i++){

      num = get_next_line(fp, buf, bufsize);
      if (num == 0) input_file_error(numProcs, count_tag, 1);

      num = sscanf(buf, ZOLTAN_ID_SPEC "%f %f", myMesh->myGlobalIDs + i,
                                                myMesh->x + i, myMesh->y + i);

      if (num != 3) input_file_error(numProcs, count_tag, 1);

    }

    gids = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * (nobj + 1));
    xcoord = (float *)malloc(sizeof(float) * (nobj + 1));
    ycoord = (float *)malloc(sizeof(float) * (nobj + 1));

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
          num = sscanf(buf, ZOLTAN_ID_SPEC "%f %f", gids+j, xcoord+j, ycoord+j);

          if (num != 3) input_file_error(numProcs, count_tag, i);
        }
      }

      MPI_Send(&nobj, 1, MPI_INT, i, count_tag, MPI_COMM_WORLD);
      MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

      if (nobj > 0){
        MPI_Send(gids, nobj, ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
        MPI_Send(xcoord, nobj, MPI_FLOAT, i, x_tag, MPI_COMM_WORLD);
        MPI_Send(ycoord, nobj, MPI_FLOAT, i, y_tag, MPI_COMM_WORLD);
      }
    }

    free(gids);
    free(xcoord);
    free(ycoord);
    fclose(fp);
    free(buf);

    /* signal all procs it is OK to go on */
    ack = 0;
    for (i=1; i < numProcs; i++){
      MPI_Send(&ack, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
  else{

    MPI_Recv(&myMesh->numMyPoints, 1, MPI_INT, 0, count_tag, MPI_COMM_WORLD, &status);
    ack = 0;
    if (myMesh->numMyPoints > 0){
      myMesh->myGlobalIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * myMesh->numMyPoints);
      myMesh->x = (float *)malloc(sizeof(float) * myMesh->numMyPoints);
      myMesh->y = (float *)malloc(sizeof(float) * myMesh->numMyPoints);
      MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);
      MPI_Recv(myMesh->myGlobalIDs, myMesh->numMyPoints, ZOLTAN_ID_MPI_TYPE, 0,
               id_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(myMesh->x, myMesh->numMyPoints, MPI_FLOAT, 0,
               x_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(myMesh->y, myMesh->numMyPoints, MPI_FLOAT, 0,
               y_tag, MPI_COMM_WORLD, &status);
    }
    else if (myMesh->numMyPoints == 0){
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
void showSimpleMeshPartitions(int myProc, int numIDs, ZOLTAN_ID_PTR GIDs, int *parts)
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

