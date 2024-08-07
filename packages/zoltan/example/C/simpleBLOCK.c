// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/* Basic example of using Zoltan to compute a quick partitioning
** of a set of objects.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "zoltan.h"

/* Name of file containing objects to be partitioned */

static char *global_fname="objects.txt";

/* Structure to hold object data */

typedef struct{
  int numGlobalObjects;
  int numMyObjects;
  ZOLTAN_ID_TYPE *myGlobalIDs;
} OBJECT_DATA;

static int get_number_of_objects(void *data, int *ierr);

static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
 
static int get_next_line(FILE *fp, char *buf, int bufsize);

static void input_file_error(int numProcs, int tag, int startProc);

static void showSimpleMeshPartitions(int myProc, int numIDs, ZOLTAN_ID_TYPE *GIDs, int *parts);

static void read_input_objects(int myRank, int numProcs, char *fname, OBJECT_DATA *myData);

int main(int argc, char *argv[])
{
  int rc, i;
  int myRank, numProcs;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids;
  ZOLTAN_ID_PTR exportGlobalGids, exportLocalGids; 
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  int *parts = NULL;

  FILE *fp;
  OBJECT_DATA myData;

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
  ** Read objects from input file and distribute them unevenly
  ******************************************************************/

  fp = fopen(global_fname, "r");
  if (!fp){
    if (myRank == 0) fprintf(stderr,"ERROR: Can not open %s\n",global_fname);
    MPI_Finalize();
    exit(1);
  }
  fclose(fp);

  read_input_objects(myRank, numProcs, global_fname, &myData);

  /******************************************************************
  ** Create a Zoltan library structure for this instance of load
  ** balancing.  Set the parameters and query functions.
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "LB_METHOD", "BLOCK");  /* Zoltan method: "BLOCK" */
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); /* global ID is 1 integer */
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1"); /* local ID is 1 integer */
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); /* we omit object weights */

  /* Query functions */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &myData);
  Zoltan_Set_Obj_List_Fn(zz, get_object_list, &myData);

  /******************************************************************
  ** Call Zoltan to partition the objects.
  ******************************************************************/

  rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partitioning was changed, 0 otherwise */ 
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of objects to be sent to me */
        &importGlobalGids,  /* Global IDs of objects to be sent to me */
        &importLocalGids,   /* Local IDs of objects to be sent to me */
        &importProcs,    /* Process rank for source of each incoming object */
        &importToPart,   /* New partition for each incoming object */
        &numExport,      /* Number of objects I must send to other processes*/
        &exportGlobalGids,  /* Global IDs of the objects I must send */
        &exportLocalGids,   /* Local IDs of the objects I must send */
        &exportProcs,    /* Process to which I send each of the objects */
        &exportToPart);  /* Partition to which each object will belong */

  if (rc != ZOLTAN_OK){
    printf("Error in Zoltan library\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  /******************************************************************
  ** Visualize the object partitioning before and after calling Zoltan.
  **
  ** In this example, partition number equals process rank.
  ******************************************************************/

  parts = (int *)malloc(sizeof(int) * myData.numMyObjects);

  for (i=0; i < myData.numMyObjects; i++){
    parts[i] = myRank;
  }

  if (myRank== 0){
    printf("\nObject partition assignments before calling Zoltan\n");
  }

  showSimpleMeshPartitions(myRank, myData.numMyObjects, myData.myGlobalIDs, parts);

  for (i=0; i < numExport; i++){
    parts[exportLocalGids[i]] = exportToPart[i];
  }

  if (myRank == 0){
    printf("Object partition assignments after calling Zoltan\n");
  }

  showSimpleMeshPartitions(myRank, myData.numMyObjects, myData.myGlobalIDs, parts);

  /******************************************************************
  ** Free the arrays allocated by Zoltan_LB_Partition, and free
  ** the storage allocated for the Zoltan structure.
  ******************************************************************/

  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);
  if (myData.numMyObjects) {
    free(parts);
    free(myData.myGlobalIDs);
  }

  MPI_Finalize();

  return 0;
}
/* Application defined query functions */

static int get_number_of_objects(void *data, int *ierr)
{
  OBJECT_DATA *objects = (OBJECT_DATA *)data;
  *ierr = ZOLTAN_OK;
  return objects->numMyObjects;
}

static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
int i;
  OBJECT_DATA *objects = (OBJECT_DATA *)data;
  *ierr = ZOLTAN_OK;

  /* In this example, return the IDs of our objects, but no weights.
   * Zoltan will assume equally weighted objects.
   */

  for (i=0; i<objects->numMyObjects; i++){
    globalID[i] = objects->myGlobalIDs[i];
    localID[i] = i;
  }
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

void showSimpleMeshPartitions(int myProc, int numIDs, ZOLTAN_ID_TYPE *GIDs, int *parts)
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

/* Proc 0 reads the objects in the input file and divides them across processes */

void read_input_objects(int myRank, int numProcs, char *fname, OBJECT_DATA *myData)
{
char *buf;
int bufsize = 512;
int num, nobj, remainingObj, ack=0;
int i, j;
ZOLTAN_ID_TYPE *gids;
FILE *fp;
MPI_Status status;
int obj_ack_tag = 5, obj_count_tag = 10, obj_id_tag = 15;

  if (myRank == 0){

    buf = (char *)malloc(sizeof(char) * bufsize);
    fp = fopen(fname, "r");

    num = get_next_line(fp, buf, bufsize);
    if (num == 0) input_file_error(numProcs, obj_count_tag, 1);
    num = sscanf(buf, "%d", &myData->numGlobalObjects);
    if (num != 1) input_file_error(numProcs, obj_count_tag, 1);

    if (numProcs > 1){
      nobj = myData->numGlobalObjects / 2;
      remainingObj = myData->numGlobalObjects - nobj;
    }
    else{
      nobj = myData->numGlobalObjects;
      remainingObj = 0;
    }

    myData->myGlobalIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nobj);
    myData->numMyObjects = nobj;

    for (i=0; i < nobj; i++){

      num = get_next_line(fp, buf, bufsize);
      if (num == 0) input_file_error(numProcs, obj_count_tag, 1);
      num = sscanf(buf, ZOLTAN_ID_SPEC , myData->myGlobalIDs + i);
      if (num != 1) input_file_error(numProcs, obj_count_tag, 1);
  
    }

    gids = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * (nobj + 1));

    for (i=1; i < numProcs; i++){
    
      if (remainingObj > 1){
        nobj = remainingObj / 2;
        remainingObj -= nobj;
      }
      else if (remainingObj == 1){
        nobj = 1;
        remainingObj = 0;
      }
      else{
        nobj = 0;
      }

      if ((i == numProcs - 1) && (remainingObj > 0))
        nobj += remainingObj;

      if (nobj > 0){
        for (j=0; j < nobj; j++){
          num = get_next_line(fp, buf, bufsize);
          if (num == 0) input_file_error(numProcs, obj_count_tag, i);
          num = sscanf(buf, ZOLTAN_ID_SPEC, gids + j);
          if (num != 1) input_file_error(numProcs, obj_count_tag, i);
        }
      }

      MPI_Send(&nobj, 1, MPI_INT, i, obj_count_tag, MPI_COMM_WORLD);
      MPI_Recv(&ack, 1, MPI_INT, i, obj_ack_tag, MPI_COMM_WORLD, &status);

      if (nobj > 0)
        MPI_Send(gids, nobj, ZOLTAN_ID_MPI_TYPE, i, obj_id_tag, MPI_COMM_WORLD);
      
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

    MPI_Recv(&myData->numMyObjects, 1, MPI_INT, 0, obj_count_tag, MPI_COMM_WORLD, &status);
    ack = 0;
    if (myData->numMyObjects > 0){
      myData->myGlobalIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * myData->numMyObjects);
      MPI_Send(&ack, 1, MPI_INT, 0, obj_ack_tag, MPI_COMM_WORLD);
      MPI_Recv(myData->myGlobalIDs, myData->numMyObjects, ZOLTAN_ID_MPI_TYPE, 0, 
               obj_id_tag, MPI_COMM_WORLD, &status);
    }
    else if (myData->numMyObjects == 0){
      MPI_Send(&ack, 1, MPI_INT, 0, obj_ack_tag, MPI_COMM_WORLD);
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
