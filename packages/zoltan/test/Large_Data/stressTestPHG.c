/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
/**************************************************************
* Stress test that can create a very large hypergraph to test
* the large memory problems.
*
* When compiled to use 64-bit global IDs, the option
*   --use_high_order_bits
* asks that the 32 high order bits be used as part of the
* global ID.
*
* This test uses hypergraph query functions unless the option
*    --use_graph_queries
* is chosen.
*
* The approximate number of vertices in the graph is given with
* the option
*   --size=
*
* The vertices have unit weights unless the option
*   --use_varying_weights
* is chosen.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <signal.h>
#include <getopt.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include "zz_const.h"

#ifndef M_PI
#define M_PI 3.14159
#endif

extern long long atoll(const char *);

static int myRank, numProcs, numPins, nborCount;
static int *vertex_part = NULL;

static double mbytes=0;

#define NUM_GLOBAL_VERTICES     2500
#define VERTEX_WEIGHT_DIMENSION 1
#define EDGE_WEIGHT_DIMENSION 1

static int64_t gid_base = 0x000000000;
static int64_t high_order_bit = 0x100000000;
static long long numGlobalVertices;
static int vwgt_dim=1;
static int unit_weights=1;

static ZOLTAN_ID_TYPE *vertexGIDs=NULL;
static float *vwgts=NULL;
static double *x_coord;
static double *y_coord;
static int lid_quarter[5];
static int edgeCutCost;
static int numMyVertices;

static int create_a_graph();
extern void Zoltan_write_linux_meminfo(int, char *, int);

#define proc_vertex_gid(proc, lid) (vertexGIDs[proc] + lid)

static void usage()
{
  printf( "\nUsage: --do_rcb\n");
  printf( "\n       --use_rcb_coarsening\n");
  printf( "\n       --hybrid_reduction_levels={num levels}\n");
  printf( "\n       --hybrid_reduction_factor={float factor}\n");
  printf( "\n       --use_high_order_bits\n");
  printf( "\n       --use_graph_queries\n");
  printf( "\n       --size={global number of vertices}\n");
  printf( "\n       --use_varying_weights\n");

  printf( "\nDefault is to do PHG without rcb coarsening.\n");
  printf( "\nDefault is to not use the 32 high order bits of the 64 bit global ID\n");
  printf( "\nDefault is to use hypergraph queries, not graph queries.\n");
  printf( "\nDefault global number of vertices is %d\n",NUM_GLOBAL_VERTICES);
  printf( "\nDefault is unit vertex weights.\n");
}

static int get_partition_quality(struct Zoltan_Struct *zz, float *cutn, float *cutl, float *imbalance)
{
ZOLTAN_HG_EVAL result;
int rc;

  rc = Zoltan_LB_Eval_HG(zz, 0, &result);

  if (rc != ZOLTAN_OK){
    fprintf(stderr,"%d Failure in Zoltan_LB_Eval_HG\n",myRank);
    return 1;
  }

  *cutn = result.cutn[EVAL_GLOBAL_SUM];
  *cutl = result.cutl[EVAL_GLOBAL_SUM];
  *imbalance = result.imbalance;

  return 0;
}

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


/* Zoltan query functions. */

static int get_number_of_vertices(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return numMyVertices;
}

static int get_vertex_coordinate(
  int lid, double *x, double *y, double *z)
{
  *z = (double)myRank;
  if (lid < lid_quarter[1]){
   *x = x_coord[lid];  
   *y = y_coord[lid];  
  }
  else if (lid < lid_quarter[2]){
   *x = x_coord[lid_quarter[2] -1 - lid];  
   *y = -y_coord[lid_quarter[2] -1 - lid];  
  }
  else if (lid < lid_quarter[3]){
   *x = -x_coord[lid - lid_quarter[2]];
   *y = -y_coord[lid - lid_quarter[2]];
  }
  else{
   *x = -x_coord[lid_quarter[4] -1 - lid];  
   *y = y_coord[lid_quarter[4] -1 - lid];  
  }
  return 0;
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

    if (!unit_weights)
      obj_wgts[i] = vwgts[i];
    else
      obj_wgts[i] = 1.0;
  }
}

static int get_coordinate_dimension(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}

static void get_geometry(void *data, int num_gid, int num_lid,
  int num_obj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int num_dim,
  double *geom_vec, int *ierr)
{
  int i;
  double *c;
  if (num_dim != 3){
    *ierr = ZOLTAN_FATAL;
    return;
  }

  c = geom_vec;
  for (i=0; i < num_obj; i++){
    get_vertex_coordinate(lids[i], c, c+1, c+2);
    c += 3;
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

    /* "left" neighbor is on my process */
    if (i==0){
      left = proc_vertex_gid(myRank, num_edges - 1);
    }
    else{
      left = edgeGID[i] - 1;
    }
    pinGID[next++] = left;

    /* "right" neighbor is on my process */
    if (i==num_edges-1){
      right = proc_vertex_gid(myRank, 0);
    }
    else{
      right = edgeGID[i] + 1;
    }
    pinGID[next++] = right;

    /* vertex generating the hyperedge is on my process */
    pinGID[next++] = edgeGID[i];

    if (myRank > 0){ /* pin belongs to process just "before" me */
      pinGID[next++] = proc_vertex_gid(myRank-1, i);
    }
    if (myRank < numProcs-1){ /* pin belongs to process just "after" me */
      pinGID[next++] = proc_vertex_gid(myRank+1, i);
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

static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr)
{
int i;

  *ierr = ZOLTAN_OK;

  for (i=0; i < num_obj; i++){
    /* Every vertex has nborCount neighbors */
    numEdges[i] = nborCount - 1;
  }
}
static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nborGID, int *nborProc,
        int wgt_dim, float *ewgts, int *ierr)
{
int i;
ZOLTAN_ID_TYPE lid, before, after, left=0, right=0;
float wgt;

  *ierr = ZOLTAN_OK;

  wgt = (float)edgeCutCost / (float)(nborCount-1);

  for (i=0; i < num_obj; i++){

    lid = localID[i];

    if (lid==0){
      before = proc_vertex_gid(myRank, num_obj-1);
    }
    else{
      before = globalID[i] - 1;
    }

    if (lid==num_obj-1){
      after = proc_vertex_gid(myRank, 0);
    }
    else{
      after = globalID[i] + 1;
    }

    if (myRank > 0){
      left = proc_vertex_gid(myRank-1, lid);
    }

    if (myRank < numProcs-1){
      right = proc_vertex_gid(myRank+1, lid);
    }

    *nborGID++ = before;  *nborProc++ = myRank;  *ewgts++ = wgt;
    *nborGID++ = after;   *nborProc++ = myRank;  *ewgts++ = wgt;
    if (myRank > 0){
      *nborGID++ = left;   *nborProc++ = myRank-1;  *ewgts++ = wgt;
    }
    if (myRank < numProcs-1){
      *nborGID++ = right;  *nborProc++ = myRank+1;  *ewgts++ = wgt;
    }
  }
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


int main(int argc, char *argv[])
{
  int i, rc, status;
  float ver;
  struct option opts[20];
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  double localMBytes, min, max, avg;
  float cutn[2], cutl[2], imbalance[2];
  struct timeval t1, t2;
  char factorBuf[64], levelBuf[64];
  time_t startusecs, endusecs, diff;

  int use_hg = 1;
  int use_graph = 0;
  int use_rcb_with_phg = 0;
  int do_rcb =0;
  float hybrid_reduction_factor =.1;
  int hybrid_reduction_levels  =2;

  cutn[0] = cutn[1] = cutl[0] = cutl[1] = imbalance[0] = imbalance[1] = 0.0;

#ifdef LINUX_HOST
  signal(SIGSEGV, meminfo_signal_handler);
  signal(SIGINT, meminfo_signal_handler);
  signal(SIGTERM, meminfo_signal_handler);
  signal(SIGABRT, meminfo_signal_handler);
  signal(SIGFPE, meminfo_signal_handler);
#endif

  /******************************************************************
  ** Initialize MPI and Zoltan
  ******************************************************************/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    if (myRank == 0) printf("sorry...\n");
    MPI_Finalize();
    return 1;
  }



  /******************************************************************
  ** Arguments
  ******************************************************************/

  opts[0].name = "do_rcb";
  opts[0].has_arg = 0;
  opts[0].flag = NULL;
  opts[0].val = 1;

  opts[1].name = "hybrid_reduction_levels";
  opts[1].has_arg = 1;
  opts[1].flag = NULL;
  opts[1].val = 2;

  opts[2].name = "hybrid_reduction_factor";
  opts[2].has_arg = 1;
  opts[2].flag = NULL;
  opts[2].val = 3;

  opts[3].name = "use_rcb_coarsening";
  opts[3].has_arg = 0;
  opts[3].flag = NULL;
  opts[3].val = 4;

  opts[4].name = "use_high_order_bits";
  opts[4].has_arg = 0;
  opts[4].flag = NULL;
  opts[4].val = 5;

  opts[5].name = "use_hypergraph_queries";
  opts[5].has_arg = 0;
  opts[5].flag = NULL;
  opts[5].val = 6;

  opts[6].name = "use_graph_queries";
  opts[6].has_arg = 0;
  opts[6].flag = NULL;
  opts[6].val = 7;

  opts[7].name = "help";
  opts[7].has_arg = 0;
  opts[7].flag = NULL;
  opts[7].val = 8;

  opts[8].name = "size";
  opts[8].has_arg = 1;
  opts[8].flag = NULL;
  opts[8].val = 9;

  opts[9].name = "use_varying_weights";
  opts[9].has_arg = 0;
  opts[9].flag = NULL;
  opts[9].val = 10;

  opts[10].name = 0;
  opts[10].has_arg = 0;
  opts[10].flag = NULL;
  opts[10].val = 0;

  while (1){
    rc = getopt_long_only(argc, argv, "",  opts, NULL);

    if (rc == '?'){
      MPI_Barrier(MPI_COMM_WORLD);
      if (myRank == 0) usage();
      MPI_Finalize();
      return 1;
    }
    else if (rc == 1){
      do_rcb = 1;
    }
    else if (rc == 2){
      hybrid_reduction_levels= atoi(optarg);
    }
    else if (rc == 3){
      hybrid_reduction_factor= atof(optarg);
    }
    else if (rc == 4){
      use_rcb_with_phg = 1;
    }
    else if (rc == 5){
      if (sizeof(ZOLTAN_ID_TYPE) != 8){
        if (myRank == 0){
          printf("The ZOLTAN_ID_TYPE has only 4 bytes, so we ignore --use_high_order_bits.\n");
        }
      }
      else{
        gid_base = high_order_bit;
      }
    }
    else if (rc == 6){
      use_hg = 1;
      use_graph = 0;
    }
    else if (rc == 7){
      use_hg = 0;
      use_graph = 1;
    }
    else if (rc == 10){
      unit_weights = 0;
    }
    else if (rc == 8){
      if (myRank == 0) usage();
      MPI_Finalize();
      return 0;
    }
    else if (rc == 9){
      numGlobalVertices = (size_t)atoll(optarg);
    }
    else if (rc <= 0){
      break;
    }
  }

  sprintf(factorBuf, "%f", hybrid_reduction_factor);
  sprintf(levelBuf, "%d", hybrid_reduction_levels);

  if (!myRank){
    printf("|========================================================|\n");
    printf("Creating graph with approximately %lld vertices.\n",numGlobalVertices);
    if (do_rcb)
      printf("Using method RCB.\n");
    else if (use_rcb_with_phg){
      printf("Using method PHG with geometry.\n");
      printf("HYBRID_REDUCTION_FACTOR %s\n",factorBuf);
      printf("HYBRID_REDUCTION_LEVELS %s\n",levelBuf);
    }
    else
      printf("Using method PHG.\n");

    if (!do_rcb && use_hg)
      printf("Using hypergraph queries.\n");
    if (!do_rcb && use_graph)
      printf("Using graph queries.\n");
    if (gid_base != 0)
      printf("Using the high order 32 bits of the global ID fields.\n");
    if (unit_weights)
      printf("Vertex weights will all be the same.\n");
    else
      printf("Vertex weights will vary widely.\n");
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

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* number of weights per vertex */
  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");/* number of weights per hyperedge */
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");

  if (do_rcb){
    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");   /* partitioning method */
  }
  else{
    Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
    Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
    if (use_rcb_with_phg){
      Zoltan_Set_Param(zz, "PHG_COARSENING_METHOD", "RCB");
      Zoltan_Set_Param(zz, "HYBRID_REDUCTION_FACTOR", factorBuf);
      Zoltan_Set_Param(zz, "HYBRID_REDUCTION_LEVELS", levelBuf);
    }
  }

  /* PHG parameters  - see the Zoltan User's Guide for many more
   *   (The "REPARTITION" approach asks Zoltan to create a partitioning that is
   *    better but is not too far from the current partitioning, rather than partitioning 
   *    from scratch.  It may be faster but of lower quality that LB_APPROACH=PARTITION.)

  Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
  */

  /* Application defined query functions */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, NULL);

  if (use_hg){
    Zoltan_Set_HG_Size_CS_Fn(zz, get_local_hypergraph_size, NULL);
    Zoltan_Set_HG_CS_Fn(zz, get_local_hypergraph, NULL);
    Zoltan_Set_HG_Size_Edge_Wts_Fn(zz, get_edge_weights_size, NULL);
    Zoltan_Set_HG_Edge_Wts_Fn(zz, get_edge_weights, NULL);
  }

  if (use_graph){
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list,  NULL);
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list,  NULL);
  }

  if (do_rcb || use_rcb_with_phg){
    Zoltan_Set_Num_Geom_Fn(zz, get_coordinate_dimension,  NULL);
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry,  NULL);
  }

  Zoltan_Set_Part_Multi_Fn(zz, get_partition_list, NULL);

  rc = get_partition_quality(zz, cutn+0, cutl+0, imbalance+0);    /* partition quality */

  /******************************************************************
  ** Zoltan can now partition the vertices of hypergraph.
  ** In this simple example, we assume the number of partitions is
  ** equal to the number of processes.  Process rank 0 will own
  ** partition 0, process rank 1 will own partition 1, and so on.
  ******************************************************************/

  Zoltan_Memory_Debug(1);
  Zoltan_Memory_Reset(ZOLTAN_MEM_STAT_MAXIMUM);

  gettimeofday(&t1, NULL);

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

  gettimeofday(&t2, NULL);
  localMBytes = (double)Zoltan_Memory_Usage(ZOLTAN_MEM_STAT_MAXIMUM)/(1024.0*1024);

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    return 1;
  }

  /******************************************************************
  ** "Migrate"
  ******************************************************************/

  vertex_part = (int *)malloc(sizeof(int) * numMyVertices);

  if (!vertex_part){
    printf("sorry memory error...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    return 1;
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
    printf("\nBALANCE after running Zoltan with %d processes.\n", numProcs);
  }

  rc = get_partition_quality(zz, cutn+1, cutl+1, imbalance+1);   /* new partition quality */

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

  if (vertex_part) free(vertex_part);
  if (vwgts) free(vwgts);

  MPI_Allreduce(&localMBytes, &avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  avg /= (double)numProcs;
  MPI_Allreduce(&localMBytes, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&localMBytes, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  MPI_Finalize();

  if (myRank == 0){
    printf("Imbalance:   BEFORE %f    AFTER %f\n",imbalance[0],imbalance[1]);
    printf("CUTN :   BEFORE %f    AFTER %f\n",cutn[0],cutn[1]);
    printf("CUTL :   BEFORE %f    AFTER %f\n",cutl[0],cutl[1]);
    printf("Total MBytes in use by test while Zoltan is running: %12.3f\n",
             mbytes/(1024.0*1024));
    printf("Min/Avg/Max of maximum MBytes in use by Zoltan:    %12.3f / %12.3f / %12.3f\n",
             min, avg, max);

    startusecs = (t1.tv_sec * 1e6) + t1.tv_usec;
    endusecs = (t2.tv_sec * 1e6) + t2.tv_usec;
    diff = endusecs - startusecs;
    printf("Time spent in partitioning (s): %f\n",(double)diff/1e6);
  }

  status = 0;

  if ((imbalance[1] >= imbalance[0]) && (cutn[1] >= cutn[0])){
    if (myRank == 0)
      printf("FAILED: partition quality did not improve\n");
    status = 1;
  }

  return status;
}

/* Create a simple graph.  The vertices of the graph are the objects for Zoltan to partition.
 * The graph itself can be used to generate a hypergraph in this way: A vertex and all of its
 * neighbors represent a hyperedge.
 *
 * Each process' vertices are on a circle in the x-y plane.  The vertices connect "left" and "right" to the vertices of the processes below and above.
 */
static int create_a_graph()
{
  int i, j, nvtxs, num4, mid;
  int midProc;
  int random_weights = (vwgt_dim > 0);
  int heavyPart = (myRank % 3 == 0);
  double radius, radianStride, radians;

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

  /******************************************************************
  ** Some coordinates.
  ******************************************************************/

  radius = (double)num4/4.0;
  radianStride = M_PI/2.0/(num4 + 1);
  radians = M_PI/2.0 - radianStride;
  x_coord = (double *)malloc(sizeof(double) * num4);
  y_coord = (double *)malloc(sizeof(double) * num4);

  for (i=0; i < num4; i++){
    y_coord[i] = radius * sin(radians);
    x_coord[i] = radius * cos(radians);
    radians -= radianStride;
  }

  lid_quarter[0] = 0;
  lid_quarter[1] = num4;
  lid_quarter[2] = num4*2;
  lid_quarter[3] = num4*3;
  lid_quarter[4] = num4*4;

  /******************************************************************
  ** Check again that this test makes sense.
  ******************************************************************/

  if (gid_base > 0 && (numGlobalVertices/2) > gid_base){
    /* half of the vertex gids should be below gid_base, the other half at or above */
    if (myRank == 0){
      printf("ERROR: When using higher order bits, the global number of vertices must be < %lld\n",
                 (gid_base)*2);
    }
    MPI_Finalize();
    return 1;
  }

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

  vertexGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numProcs );

  mbytes += sizeof(ZOLTAN_ID_TYPE) * numProcs;

  if (gid_base > 0){
    midProc = numProcs / 2;
    vertexGIDs[0] = 0;
    vertexGIDs[midProc] = gid_base;

    for (i=1,j=midProc+1; (i < midProc) || (j < numProcs); i++,j++){
      if (i < midProc) vertexGIDs[i] = vertexGIDs[i-1] + nvtxs;
      if (j < numProcs) vertexGIDs[j] = vertexGIDs[j-1] + nvtxs;
    }
  }

  else{

    vertexGIDs[0] = 0;
    for (i=1; i < numProcs; i++){
      vertexGIDs[i] = vertexGIDs[i-1] + nvtxs;
    }
  }

  if (myRank == 0){
    printf("Hypergraph will have %lld hyperedges, %d on each process\n", numGlobalVertices, nvtxs);
  }

  if (vwgt_dim == 0) 
    vwgt_dim = 1;

  if (!unit_weights){
    vwgts = (float *)calloc( vwgt_dim * nvtxs, sizeof(float));
    if (!vwgts) return 1;
  
    mbytes += vwgt_dim * nvtxs * sizeof(float);
    
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
  }

  /* Make edge cut costs higher in the "center" */

  mid = (int)(numProcs / 2.0);
  edgeCutCost = ((myRank < mid) ? myRank + 1 : numProcs - myRank);

  return 0;
}
