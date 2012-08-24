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
/************************************************************
* This is called a stress test because it builds an
* arbitrarily large graph.  It tests the HIER_ASSIST
* option to hierarchical partitioning.
*
* TODO:
* Create a function that performs communication where comm
* volume is proportional to the graph edge weight.  This is
* to test the value of partitioning to network hierarchy.
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <signal.h>
#include <getopt.h>
#include "zz_const.h"

#define NUM_GLOBAL_VERTICES     2500    /* default */

static int verbose=0;
static int myRank, numProcs, numMyPins;
static long cylCount, cylSize, numGlobalVertices;
static long myFirstGID; 
static int numMyVertices;
static float heavyCommWeight = 100;
static float lowCommWeight = 1;

ZOLTAN_ID_TYPE *vtxGID = NULL;
ZOLTAN_ID_TYPE *nborGID = NULL;
int *nborIndex = NULL;
int *nborProc = NULL;
float *edgeWgt = NULL;

struct Zoltan_DD_Struct *dd=NULL;
static void debug(struct Zoltan_Struct *zz, char *s, int stop);
static void usage();

static void check_error_status(int status, char *s)
{
int gstatus;

  MPI_Allreduce(&status, &gstatus, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (gstatus > 0){
    if (myRank == 0){
      fprintf(stderr,"Error: %s\n",s);
    }
    MPI_Finalize();
    exit(1);
  }
}

static void free_graph()
{
  if (dd) Zoltan_DD_Destroy(&dd);

  if (vtxGID) free(vtxGID);
  if (nborIndex) free(nborIndex);
  if (nborGID) free(nborGID);
  if (nborProc) free(nborProc);
  if (edgeWgt) free(edgeWgt);

  vtxGID = nborGID = NULL;
  nborIndex = nborProc = NULL;
  edgeWgt = NULL;
  dd=NULL;
}

static void gid_location(long gid, long *cylID, long *ringID)
{
  ldiv_t result;
  result = ldiv(gid, cylSize);
  *cylID = result.quot;  /* cylinder ID */
  *ringID = result.rem;  /* ring ID */
}
static long gid_up(long gid)
{
  long cylID, ringID;
  gid_location(gid, &cylID, &ringID); 

  if (ringID == cylSize-1){
    return cylID * cylSize;
  }
  else{
    return gid + 1;
  }
}
static long gid_down(long gid)
{
  long cylID, ringID;
  gid_location(gid, &cylID, &ringID); 

  if (ringID == 0){
    return gid + cylSize - 1;
  }
  else{
    return gid - 1;
  }
}
static long gid_left(long gid)
{
  long cylID, ringID;
  gid_location(gid, &cylID, &ringID); 

  if (cylID == 0){
    return -1;
  }

  return gid - cylSize;
}
static long gid_right(long gid)
{
  long cylID, ringID;
  gid_location(gid, &cylID, &ringID); 

  if (cylID == cylCount - 1){
    return -1;
  }

  return gid + cylSize;
}
static int num_neighbors(long gid)
{
  long cylID, ringID;
  int nnbors = 2;   /* up and down */

  gid_location(gid, &cylID, &ringID); 

  if (cylID > 0) nnbors++;   /* left */
  if (cylID < cylCount-1) nnbors++;   /* right */

  return nnbors;
}

static int get_nbor_info(long gid, long *nbors, float *wgt)
{
  int i=0;
  long n;

  /* TODO   Vary the weights more to see more of a difference */

  nbors[i] = gid_up(gid);
  wgt[i] = lowCommWeight;

  nbors[++i] = gid_down(gid);
  wgt[i] = lowCommWeight;

  n = gid_left(gid);
  if (n >= 0){
    nbors[++i] = n;
    wgt[i] = heavyCommWeight;
  }

  n = gid_right(gid);
  if (n >= 0){
    nbors[++i] = n;
    wgt[i] = heavyCommWeight;
  }

  return i+1;
}

static int create_a_graph()
{
  int rc, i, sum, n, j;
  float wgts[4]; 
  long gid, nbors[4], count[2];
  ldiv_t result;
  float c;

  c = sqrt((float)numGlobalVertices);

  c = (c < 3.0) ? 3.0 : c;

  cylCount = cylSize = (long)c;

  numGlobalVertices = cylCount * cylSize;

  result = ldiv(numGlobalVertices, (long)numProcs);

  numMyVertices = (int)result.quot + (myRank < result.rem ? 1 : 0);

  count[0] = (long)numMyVertices;

  MPI_Scan(count+0, count+1, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  myFirstGID = count[1] - count[0];

  numMyPins = 0;

  for (i=0, gid = myFirstGID; i < numMyVertices; i++, gid++){
    numMyPins += num_neighbors(gid);
  }

  vtxGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numMyVertices);
  nborIndex = (int *)malloc(sizeof(int) * (numMyVertices + 1));
  nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numMyPins);
  nborProc = (int *)malloc(sizeof(int) * numMyPins);
  edgeWgt = (float *)malloc(sizeof(float) * numMyPins);

  if (numMyPins && !(vtxGID || nborIndex || nborGID || nborProc || edgeWgt)){
    fprintf(stderr,"%d out of memory\n",myRank);
    return 1;
  }

  nborIndex[0] = 0;

  for (i=0, gid=myFirstGID, n=0; i < numMyVertices; i++, gid++){
    vtxGID[i] = (ZOLTAN_ID_TYPE)gid;
    sum = get_nbor_info(gid, nbors, wgts);

    for (j=0; j < sum; j++, n++){
      nborGID[n] = (ZOLTAN_ID_TYPE)nbors[j];
      edgeWgt[n] = wgts[j];

      if (nborGID[n] < myFirstGID)
        nborProc[n] = myRank - 1;
      else if (nborGID[n] >= myFirstGID + numMyVertices)
        nborProc[n] = myRank + 1;
      else
        nborProc[n] = myRank;
    }
    nborIndex[i+1] = nborIndex[i] + sum; 
  }

  rc = Zoltan_DD_Create(&dd, MPI_COMM_WORLD, 1, 0, 0, (int)(numGlobalVertices / numProcs), 0);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    fprintf(stderr,"%d DD Create failure\n",myRank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  rc = Zoltan_DD_Update(dd, vtxGID, NULL, NULL, NULL, numMyVertices);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    fprintf(stderr,"%d DD Update failure in create\n",myRank);
    return 1;
  }

  return 0;
}

static int reallocate_buffers(int numNewVertices, int numNewPins)
{
  int status = 0;
  ZOLTAN_ID_TYPE *idbuf=NULL;
  int *ibuf=NULL;
  float *fbuf=NULL;

  if (verbose) MPI_Barrier(MPI_COMM_WORLD);

  if (numNewVertices > numMyVertices){   /* avoid realloc bug */
    idbuf = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numNewVertices);
    if (!idbuf) return 1;
    memcpy(idbuf, vtxGID, sizeof(ZOLTAN_ID_TYPE) * numMyVertices);
    free(vtxGID);
    vtxGID = idbuf; 
    if (verbose){
      printf("(%d) vtxGID allocated for %d vertices\n",myRank,numNewVertices);
    }

    ibuf = (int *)malloc(sizeof(int) * (numNewVertices+1));
    if (!ibuf) return 1;
    memcpy(ibuf, nborIndex, sizeof(int) * (1 +numMyVertices));
    free(nborIndex);
    nborIndex = ibuf; 
    if (verbose){
      printf("(%d) nborIndex allocated for %d indices into nbor array\n",myRank,numNewVertices+1);
    }
  }

  if (numNewPins > numMyPins){
    idbuf = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numNewPins);
    if (!idbuf) return 1;
    memcpy(idbuf, nborGID, sizeof(ZOLTAN_ID_TYPE) * numMyPins);
    free(nborGID);
    nborGID = idbuf; 
    if (verbose){
      printf("(%d) nborGID allocated for %d neighbor IDs\n",myRank,numNewPins);
    }

    ibuf = (int *)malloc(sizeof(int) * numNewPins);
    if (!ibuf) return 1;
    memcpy(ibuf, nborProc, sizeof(int) * numMyPins);
    free(nborProc);
    nborProc = ibuf; 
    if (verbose){
      printf("(%d) nborProc allocated for %d process IDs\n",myRank,numNewPins);
    }

    fbuf = (float *)malloc(sizeof(float) * numNewPins);
    if (!fbuf) return 1;
    memcpy(fbuf, edgeWgt, sizeof(float) * numMyPins);
    free(edgeWgt);
    edgeWgt = fbuf; 
    if (verbose){
      printf("(%d) edgeWgt allocated for %d edge weights\n",myRank,numNewPins);
    }
  }

  if (verbose) {
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return status;
}

static int migrate_graph(int num_exports, int num_imports, ZOLTAN_ID_TYPE *export_lids, ZOLTAN_ID_TYPE *import_gids)
{
  int i, j, k, nextv, nextp, npins, numNewVertices, numNewPins, sum;
  long nbors[4];
  float wgts[4];
  int rc, startlocv, startloce;
  int status = 0;

  numNewVertices = numMyVertices - num_exports + num_imports;

  for (i=0; i < num_exports; i++){
    vtxGID[export_lids[i]] = ZOLTAN_ID_INVALID;
  }

  for (i=0, nextv=0, nextp=0; i < numMyVertices; i++){
    npins = nborIndex[i+1] - nborIndex[i];
    if (vtxGID[i] != ZOLTAN_ID_INVALID){
      if (i > nextv){
        vtxGID[nextv] = vtxGID[i];
        for (j=nborIndex[i], k=0; j < nborIndex[i+1]; j++, k++){
          nborGID[nextp+k] = nborGID[j];
          edgeWgt[nextp+k] = edgeWgt[j];
          /* skip nborProc because we don't know what it is yet */
        }
        nborIndex[nextv+1] = nborIndex[nextv] + npins;
      }

      nextv++;
      nextp += npins;
    }
  }

  numNewPins = nextp;

  startlocv = nextv;
  startloce = nextp;

  for (i=0; i < num_imports; i++){
    numNewPins += num_neighbors(import_gids[i]);
  }

  status = reallocate_buffers(numNewVertices, numNewPins);

  if (status == 0){
    for (i=0; i < num_imports; i++, nextv++){
      vtxGID[nextv] = import_gids[i];
      sum = get_nbor_info(import_gids[i], nbors, wgts);
  
      for (j=0; j < sum; j++, nextp++){
        nborGID[nextp] = (ZOLTAN_ID_TYPE)nbors[j];
        edgeWgt[nextp] = wgts[j];
      }
      nborIndex[nextv+1] = nborIndex[nextv] + sum; 
    }
  }
  else{
    fprintf(stderr,"memory allocation failure in reallocate buffers\n");
    return 1;
  }

  numMyVertices = numNewVertices;
  numMyPins = numNewPins;

  rc = Zoltan_DD_Update(dd, vtxGID+startlocv, NULL, NULL, NULL, num_imports);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    status = 1;
  }

  rc = Zoltan_DD_Find(dd, nborGID+startloce, NULL, NULL, NULL, numNewPins - startloce, nborProc+startloce);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    status = 1;
  }

  return status;
}

static float get_edge_cut_weight(struct Zoltan_Struct *zz)
{
ZOLTAN_GRAPH_EVAL result;
int rc;

  rc = Zoltan_LB_Eval_Graph(zz, 0, &result);

  if (rc != ZOLTAN_OK){
    fprintf(stderr,"%d Failure in Zoltan_LB_Eval_Graph\n",zz->Proc);
    return -1.0;
  }

  return result.cut_wgt[EVAL_GLOBAL_SUM];
}

void time_communication(double *t)
{
   /* TODO - perform communication proportional to edge weights */
  *t = 0;
}

/* Zoltan query functions. */

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
    globalID[i] = vtxGID[i];
    localID[i] = i;
    obj_wgts[i] = nborIndex[i+1] - nborIndex[i];
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
    numEdges[i] = nborIndex[localID[i]+1] - nborIndex[localID[i]]; 
  }
}
static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nbor, int *owner, int wgt_dim, float *ewgts, int *ierr)
{
int nextv, nextp, npins, p1, i, lid;

  *ierr = ZOLTAN_OK;

  for (nextv=0, nextp=0; nextv < num_obj; nextv++){

    lid = localID[nextv];
    p1 = nborIndex[lid];
    npins = nborIndex[lid+1] - p1;

    if (num_edges[nextv] != npins){
      fprintf(stderr,"num edges != num pins\n");
      *ierr = ZOLTAN_FATAL;
      return;
    }

    for (i=0; i <  npins; i++, nextp++){
      nbor[nextp] = nborGID[p1+i];
      owner[nextp] = nborProc[p1+i];
      ewgts[nextp] = edgeWgt[p1+i];
    }
  }
}

int main(int argc, char *argv[])
{
  int rc, do_hier, status;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  int generate_files = 0;
  char *platform=NULL, *topology=NULL;
  char *graph_package=NULL;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  struct option opts[10];
  double comm_time[10];
  float cut_weight[10];
  long nvert=0;
  char *debug_level=NULL;

  status = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  Zoltan_Initialize(argc, argv, &ver);
  zz = Zoltan_Create(MPI_COMM_WORLD);

  /******************************************************************
  ** Check that this test makes sense.
  ******************************************************************/
  
  if (sizeof(long) < sizeof(ZOLTAN_ID_TYPE)){
    if (myRank == 0){
      printf("ERROR: This code assumes that a long is at least %d bytes\n",(int)sizeof(ZOLTAN_ID_TYPE));
    }
    status = 1;
  }

  check_error_status(status, "configuration error");

  /******************************************************************
  ** Initialize zoltan
  ******************************************************************/

  /* options */

  opts[0].name = "platform";
  opts[0].has_arg = 1;
  opts[0].flag = NULL;
  opts[0].val = 1;

  opts[1].name = "topology";
  opts[1].has_arg = 1;
  opts[1].flag = NULL;
  opts[1].val = 2;

  opts[2].name = "size";
  opts[2].has_arg = 1;
  opts[2].flag = NULL;
  opts[2].val = 4;

  opts[3].name = "verbose";
  opts[3].has_arg = 0;
  opts[3].flag = NULL;
  opts[3].val = 5;

  opts[4].name = "help";
  opts[4].has_arg = 0;
  opts[4].flag = NULL;
  opts[4].val = 6;

  opts[5].name = "graph_package";
  opts[5].has_arg = 1;
  opts[5].flag = NULL;
  opts[5].val = 7;

  opts[6].name = "generate_files";
  opts[6].has_arg = 0;
  opts[6].flag = NULL;
  opts[6].val = 8;

  opts[7].name = "debug_level";
  opts[7].has_arg = 1;
  opts[7].flag = NULL;
  opts[7].val = 9;

  opts[8].name = 0;
  opts[8].has_arg = 0;
  opts[8].flag = NULL;
  opts[8].val = 0;

  status = 0;

  while (1){
    rc = getopt_long_only(argc, argv, "",  opts, NULL);

    if (rc == '?'){
      MPI_Barrier(MPI_COMM_WORLD);
      if (myRank == 0) usage();
      MPI_Finalize();
      exit(0);
    }
    else if (rc == 1){
      platform = optarg;
      if (myRank == 0)
        printf( "For platform %s\n",optarg );
    }
    else if (rc == 2){
      topology = optarg;
      if (myRank == 0)
        printf( "For topology %s\n",optarg);
    }
    else if (rc == 7){
      graph_package = optarg;
      if (myRank == 0)
        printf( "Zoltan parameter GRAPH_PACKAGE = %s\n",graph_package);
    }
    else if (rc == 8){
      generate_files = 1;
      if (myRank == 0)
        printf( "Zoltan_Generate_Files will be called for each level.\n");
    }
    else if (rc == 4){
      nvert = atol(optarg);
      if (nvert < 1) status = 1;
      check_error_status(status, "--size={approximate number of vertices}");
      if (myRank == 0){
        printf( "Graph will have approximately %ld vertices.\n",nvert);
      }
    }
    else if (rc == 5){
      verbose = 1;
    }
    else if (rc == 6){
      if (myRank == 0) usage();
      MPI_Finalize();
      exit(0);
    }
    else if (rc == 9){
      debug_level = optarg;
    }
    else if (rc <= 0){
      break;
    }
  }

  if ((platform==NULL) && (topology==NULL)){
    if (myRank == 0)
      fprintf(stdout,"No platform or topology, so we'll skip hierarchical partitioning\n");
    do_hier = 0;
  }
  else if (graph_package == NULL){
    if (myRank == 0)
      fprintf(stdout,"No graph package, so we'll skip hierarchical partitioning\n");
    do_hier = 0;
  }
  else{
    do_hier = 1;
  }

  /* start */

  Zoltan_Memory_Debug(0);

  if (nvert > 0)
    numGlobalVertices = nvert;
  else
    numGlobalVertices = NUM_GLOBAL_VERTICES;

  status = create_a_graph(); 
  check_error_status(status, "creating the graph");

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "REMAP", "0");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* number of weights per vertex */
  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");/* number of weights per hyperedge */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, NULL);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, NULL);
  Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list,  NULL);
  Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list,  NULL);

  /* GRAPH PARTITION */

  Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");

  if (graph_package)
    Zoltan_Set_Param(zz, "GRAPH_PACKAGE", graph_package);

  if (verbose){
    debug(zz, "Initial graph", 0);
  }

  if (generate_files){
    rc = Zoltan_Generate_Files(zz, "flat", myRank, 0, 1, 0);
    if (rc != ZOLTAN_OK) status = 1;
    check_error_status(status, "Zoltan_Generate_Files");
  }

  /* Performance before partitioning */
  time_communication(comm_time+0); 
  cut_weight[0] = get_edge_cut_weight(zz);

  if (cut_weight[0] < 0.0) status = 1;
  check_error_status(status, "First call to get_edge_cut_weight");

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

  if (rc != ZOLTAN_OK) status = 1;
  check_error_status(status, "First call to LB_Partition");

  status = migrate_graph(numExport, numImport, exportLocalGids, importGlobalGids);
  check_error_status(status, "migration");

  if (verbose){
    debug(zz, "After flat partitioning and migration", 0);
  }

  time_communication(comm_time+1);      /* With graph partitioning */
  cut_weight[1] = get_edge_cut_weight(zz);

  if (cut_weight[1] < 0.0) status = 1;
  check_error_status(status, "Second call to get_edge_cut_weight");

  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  if (do_hier){

    /* HIERARCHICAL PARTITION */

    free_graph();
    status = create_a_graph();
    check_error_status(status, "create graph for hierarchical partitioning");
  
    Zoltan_Set_Param(zz, "LB_METHOD", "HIER");
    Zoltan_Set_Param(zz, "HIER_ASSIST", "1");
    if (generate_files){
      Zoltan_Set_Param(zz, "HIER_GENERATE_FILES", "1"); 
    }

    if (debug_level)   /* 1, 2 or 3 */
      Zoltan_Set_Param(zz, "HIER_DEBUG_LEVEL", debug_level);
    else
      Zoltan_Set_Param(zz, "HIER_DEBUG_LEVEL", "0");

    /* TODO: Suppose graph is not symmetric, and we request SYMMETRIZE.  Do we still get
     *  a "good" answer when each sub-graph in the hierarchy is symmetrized?
     */

    if (topology)
      Zoltan_Set_Param(zz, "TOPOLOGY", topology);
    else if (platform)
      Zoltan_Set_Param(zz, "PLATFORM", platform);
  
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
  
    if (rc != ZOLTAN_OK) status = 1;
    check_error_status(status, "Second call to LB_Partition");
  
    status = migrate_graph(numExport, numImport, exportLocalGids, importGlobalGids);
    check_error_status(status, "second migration");
  
    if (verbose){
      debug(zz, "After hierarchical partitioning and migration", 0);
    }
  
    time_communication(comm_time+2);      /* With hierarchical graph partitioning */
    cut_weight[2] = get_edge_cut_weight(zz);

    if (cut_weight[2] < 0.0) status = 1;
    check_error_status(status, "Third call to get_edge_cut_weight");
  
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                        &exportProcs, &exportToPart);
  }

  Zoltan_Destroy(&zz);

  free_graph();

  if (myRank == 0){
    fprintf(stdout,"Graph cut weight before partitioning: %f\n",cut_weight[0]);
    fprintf(stdout,"             after flat partitioning: %f\n",cut_weight[1]);
    if (do_hier)
      fprintf(stdout,"     after hierarchical partitioning: %f\n",cut_weight[2]);
    fflush(stdout);
  }

  if (cut_weight[1] >= cut_weight[0]){
    status = 1;
    if (zz->Proc == 0){
      fprintf(stderr,"FAILED: No improvement shown in flat partitioning");
    }
  }

  if (do_hier && (cut_weight[2] > cut_weight[0])){
    status = 1;
    if (zz->Proc == 0){
      fprintf(stderr,"FAILED: No improvement shown in hierarchical partitioning");
    }
  }


  MPI_Finalize();

  return status;
}

static void debug(struct Zoltan_Struct *zz, char *s, int stop)
{
int i,p,j, k, nedges;

  MPI_Barrier(MPI_COMM_WORLD);

  if (s && !myRank)
    fprintf(stdout,"\n\n%s\n",s);

  if (numGlobalVertices <= 100){
    MPI_Barrier(MPI_COMM_WORLD);
    for (p=0; p < numProcs; p++){
  
      if (p == myRank){
  
        if (p==0){
          fprintf(stdout,"%ld global vertices\n",numGlobalVertices); 
        }
  
        fprintf(stdout,"Partition %d, %d vertices:\n",p,numMyVertices);
        for (i=0, k=0; i < numMyVertices; i++){
          fprintf(stdout,ZOLTAN_ID_SPEC ": ",vtxGID[i]);
          nedges = nborIndex[i+1] - nborIndex[i];
           
          for (j=0; j < nedges; j++,k++){
            fprintf(stdout,ZOLTAN_ID_SPEC "/%f/%d ",nborGID[k],edgeWgt[k],nborProc[k]);
          }
          fprintf(stdout,"\n");
        
        }
        fprintf(stdout,"\n");
        fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  Zoltan_LB_Eval_Graph(zz, 1, NULL);

  MPI_Barrier(MPI_COMM_WORLD);

  if (stop){
    MPI_Finalize();
    exit(0);
  }
}

static void usage()
{
  int i;
  printf( "\nUsage: --verbose\n");
  printf( "\n       --generate_files\n");
  printf( "\n       --graph_package={parmetis|scotch|phg}\n");
  printf( "\n       --platform={glory|redsky|ctx|odin|octopi|s861036} "
                 "| --topology=desc\n");
  printf( "\n       --size={approximate global number of vertices}\n");

  printf( "\n\nA topology description is a list of integers, for example\n");
  printf( "  Dual socket, quad core: 2, 4\n");
  printf( "  Quad socket, six cores with core pairs sharing a cache: 4, 3, 2\n");
  printf( "  Dual core workstation: 2\n\n");

  printf( "The default global number of vertices is %d\n",NUM_GLOBAL_VERTICES);
}

