// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/***************************************************************
** Example of using Zoltan to compute an RIB partitioning
** of a possibly very large collection of vertices and weights.
** Really a stress test for Zoltan developers.
**
** usage  stressTestRIB [global number of vertices] [vertex weight dim] [vertex dim]
**
** Compile with a math library.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include "zz_const.h"


static int myRank, numProcs;
static double mbytes=0;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif


/* Mesh data */

#define NUM_GLOBAL_VERTICES     2500000000
#define VERTEX_WEIGHT_DIMENSION 1
#define VERTEX_DIMENSION        3

static int numLocalVertices=0;
static float *v_x=NULL;
static float *v_y=NULL;
static float *v_z=NULL;
static float *vertex_weight=NULL;
static ZOLTAN_ID_TYPE *vertex_gid=NULL;
int *vertex_part=NULL;
static ZOLTAN_ID_TYPE first_gid;

extern void Zoltan_write_linux_meminfo(int, char *, int);

static int create_vertices(ZOLTAN_GNO_TYPE gnvtxs, int ndim, int vwgt_dim, int nprocs, int rank);
static int vertexDim, vertexWeightDim;

/* Application defined query functions */

static int get_number_of_objects(void *data, int *ierr);
static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
static int get_num_geometry(void *data, int *ierr);
static void get_geometry_list(void *data, int sizeGID, int sizeLID,
             int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr);
static void get_partition_list(void *data, int sizeGID, int sizeLID, int num_obj,
        ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int*parts, int *ierr);

void meminfo_signal_handler(int sig)
{
  char msg[128];

  sprintf(msg,"(%d) Received signal %d\n",myRank,sig);

  /* Signal handler for Linux that helps us to understand */
  /* whether failure was due to insufficient memory. */

  signal(SIGINT, SIG_IGN);
  signal(SIGTERM, SIG_IGN);
  signal(SIGABRT, SIG_IGN);
  signal(SIGSEGV, SIG_IGN);
  signal(SIGFPE, SIG_IGN);

  Zoltan_write_linux_meminfo(1, msg, 0);

  exit(sig);
}


int main(int argc, char *argv[])
{
  int rc, i; 
  ZOLTAN_GNO_TYPE numGlobalVertices;
  float ver;
  char dimstring[16];
  double min, max, avg, local;

  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

#ifdef HOST_LINUX
  signal(SIGSEGV, meminfo_signal_handler);
  signal(SIGINT, meminfo_signal_handler);
  signal(SIGTERM, meminfo_signal_handler);
  signal(SIGABRT, meminfo_signal_handler);
  signal(SIGFPE, meminfo_signal_handler);
#endif

  /******************************************************************
  ** Problem size
  ******************************************************************/

  numGlobalVertices = NUM_GLOBAL_VERTICES;
  vertexWeightDim = VERTEX_WEIGHT_DIMENSION;
  vertexDim = VERTEX_DIMENSION;

  if (argc > 1){
    sscanf(argv[1], "%zd", &numGlobalVertices);
    if (argc > 2){
      vertexWeightDim = atoi(argv[2]);
      if (argc > 3){
        vertexDim = atoi(argv[3]);
      }
    }
  }

  sprintf(dimstring,"%d",vertexWeightDim);

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
    exit(1);
  }

  Zoltan_Memory_Debug(2);

  /******************************************************************
  ** Create vertices
  ******************************************************************/

  rc = create_vertices(numGlobalVertices, vertexDim, vertexWeightDim, numProcs, myRank);

  if (rc){
    fprintf(stderr,"Process rank %d: insufficient memory\n",myRank);
    MPI_Finalize();
    exit(1);
  }

  first_gid = vertex_gid[myRank];

  /******************************************************************
  ** Create a Zoltan library structure for this instance of load
  ** balancing.  Set the parameters and query functions that will
  ** govern the library's calculation.  See the Zoltan User's
  ** Guide for the definition of these and many other parameters.
  ******************************************************************/

  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "RIB");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", dimstring);
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* RIB parameters */

  Zoltan_Set_Param(zz, "RIB_OUTPUT_LEVEL", "0");

  /* Query functions, to provide geometry to Zoltan */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_object_list, NULL);
  Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, NULL);
  Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, NULL);
  Zoltan_Set_Part_Multi_Fn(zz, get_partition_list, NULL);

  /******************************************************************
  ** Zoltan can now partition the vertices in the simple mesh.
  ** In this simple example, we assume the number of partitions is
  ** equal to the number of processes.  Process rank 0 will own
  ** partition 0, process rank 1 will own partition 1, and so on.
  ******************************************************************/

  if (myRank == 0){
    printf("Run Zoltan\n");
  }

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
    if (myRank == 0)printf("sorry...\n");
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

  rc = Zoltan_LB_Eval_Balance(zz, 1, NULL);

  if (rc != ZOLTAN_OK){
    printf("sorry first LB_Eval_Balance...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  /******************************************************************
  ** Print out the balance of the new partitions.
  ******************************************************************/
 
  vertex_part = (int *)malloc(sizeof(int) * numLocalVertices);

  if (!vertex_part){
    printf("sorry memory error...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  for (i=0; i < numLocalVertices; i++){
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

  rc = Zoltan_LB_Eval_Balance(zz, 1, NULL);

  if (rc != ZOLTAN_OK){
    printf("sorry second LB_Eval_Balance...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  /******************************************************************
  ** Free the arrays allocated by Zoltan_LB_Partition, and free
  ** the storage allocated for the Zoltan structure.
  ******************************************************************/

  if (myRank == 0){
    printf("Free structures\n");
  }

  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);

  if (vertex_part) free(vertex_part);
  if (v_x) free(v_x);
  if (v_y) free(v_y);
  if (v_z) free(v_z);
  if (vertex_weight) free(vertex_weight);
  if (vertex_gid) free(vertex_gid);

  /**********************
  ** all done ***********
  **********************/

  local= (double)Zoltan_Memory_Usage(ZOLTAN_MEM_STAT_MAXIMUM)/(1024.0*1024);
  MPI_Reduce(&local, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  avg /= (double)numProcs;
  MPI_Reduce(&local, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  if (myRank == 0){
    printf("Total MBytes in use by test while Zoltan is running: %12.3lf\n",
             mbytes/(1024.0*1024));
    printf("Min/Avg/Max of maximum MBytes in use by Zoltan:    %12.3lf / %12.3lf / %12.3lf\n",
             min, avg, max);
  }

  MPI_Finalize();

  return 0;
}

/* Application defined query functions */

static int get_number_of_objects(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return numLocalVertices;
}

static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
int i,j;
float *src, *dest;

  *ierr = ZOLTAN_OK;

  src = vertex_weight;
  dest = obj_wgts;

  for (i=0; i < numLocalVertices; i++){
    globalID[i] = first_gid + i;
    localID[i] = i;
    for (j=0; j < vertexWeightDim; j++){
      *dest++ = *src++;
    }
  }
}

static int get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return vertexDim;
}

static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr)
{
int i;
double *dest = geom_vec;

  *ierr = ZOLTAN_OK;

  for (i=0;  i < num_obj ; i++){
    *dest++ = (double)v_x[localID[i]];
    if (num_dim > 1){
      *dest++ = (double)v_y[localID[i]];
      if (num_dim > 2){
        *dest++ = (double)v_z[localID[i]];
      }
    }
  }

  return;
}

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


static void initialize_vertex_global_id_info(int numMyGIDs, int numProc)
{
  int i;
  vertex_gid = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * (numProc + 1));
  vertex_gid[0] = 0;

  for (i=1; i <= numProc; i++){
    vertex_gid[i] = vertex_gid[i-1] + numMyGIDs;
  }

  mbytes += sizeof(ZOLTAN_ID_TYPE) * (numProc + 1);
}


static int create_vertices(ZOLTAN_GNO_TYPE gnvtxs, int ndim, int vwgt_dim, int nprocs, int rank)
{
  int    nvtxs, num4, i, j;
  double theta, delta, radius, m, length, step;
  int heavyProc = ((rank % 3 == 0));

  /* for simplicity coerce number of vertices on a process to a multiple of 4 */

  nvtxs = (int)(gnvtxs / nprocs);

  if (nvtxs > 4){
    num4 = nvtxs / 4;
    nvtxs = num4 * 4;
  }
  else{
    num4 = 1;
    nvtxs = 4;
  }

  gnvtxs = (ZOLTAN_GNO_TYPE)nvtxs;
  gnvtxs *= nprocs;
  numLocalVertices = nvtxs;

  if (rank == 0){
    printf("Graph will have %zd vertices, %d on each process\n", gnvtxs, nvtxs);
  }

  /* Each process has the same number of vertices.  Let's determine their global IDs */

  initialize_vertex_global_id_info(nvtxs, nprocs);

  /* Calculate vertex coordinates */
  v_x = (float *) malloc(nvtxs * sizeof(float));
  mbytes += nvtxs * sizeof(float);
  if (ndim > 1){
     v_y = (float *) malloc(nvtxs * sizeof(float));
     mbytes += nvtxs * sizeof(float);
     if (ndim > 2){
        v_z = (float *) malloc(nvtxs * sizeof(float));
        mbytes += nvtxs * sizeof(float);
     }
  }

  vertex_weight = (float *) malloc(vwgt_dim*nvtxs * sizeof(float));
  mbytes += vwgt_dim * nvtxs * sizeof(float);

  if (ndim == 1){
    /* a line */

    step = 1.0 / 500.0;
    length = (double)nvtxs * step;
    v_x[0] = length * (float)rank;

    for (i=1; i < nvtxs; i++){
      v_x[i] = v_x[i+1] + step;
    }
  }
  else if (ndim == 2){
    /* a circle */
    radius = (double)nvtxs/500.0;
    theta = (2 * M_PI ) / (double)nprocs;
    delta = theta / (double)nvtxs;
    m = (theta * rank);

    for (i=0; i < nvtxs; i++, m += delta){
      v_x[i] = radius * cos(m);
      v_y[i] = radius * sin(m);
    }
  }
  else if (ndim == 3){
    /* a cylinder */

    radius = (double)nvtxs/500.0;
    delta = M_PI_2 / (double)(num4 + 1);
    theta = delta;
    i = 0;

    while (theta < M_PI_2){
      /* points along first quadrant of a circle in the plane z=rank */
      v_x[i] = radius * cos(theta);
      v_y[i] = radius * sin(theta);
      v_z[i] = (float)rank;
      theta += delta;
      i++;
    }

    for (i=0; i < num4; i++){
      /* second quadrant */
      v_x[num4+i] = -v_x[num4 - i - 1];
      v_y[num4+i] = v_y[num4 - i - 1];
      v_z[num4+i] = (float)rank;

      /* third quadrant */
      v_x[2*num4+i] = -v_x[i];
      v_y[2*num4+i] = -v_y[i];
      v_z[2*num4+i] = (float)rank;

      /* third quadrant */
      v_x[3*num4+i] = v_x[num4 - i - 1];
      v_y[3*num4+i] = -v_y[num4 - i - 1];
      v_z[3*num4+i] = (float)rank;
    }
  }

  srand(0);
  for (i = 0; i < nvtxs; i++)  {
    if (vwgt_dim == 0) /* Unit weights if no weights were requested. */
      vertex_weight[i] = 1.0;
    else
      if (heavyProc){
        for (j = 0; j < vwgt_dim; j++)  {
            vertex_weight[i*vwgt_dim+j] = .2 + ((float) rand())/RAND_MAX;
        }
      }
      else{
        for (j = 0; j < vwgt_dim; j++)  {
            vertex_weight[i*vwgt_dim+j] = ((float) rand())/RAND_MAX;
        }
      }
  }

  return 0;
}

