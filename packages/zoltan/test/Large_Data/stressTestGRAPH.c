/**************************************************************
* Stress test for a Zoltan partitioning.
*
* Create a graph with varying edge weights.
* Time communication that corresponds with the weights.
*
* Argument is approximate number of vertices in the graph
***************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <signal.h>
#include <getopt.h>
#include "zoltan.h"
#include "zz_util_const.h"

#define NUM_GLOBAL_VERTICES     2500000    /* default */

static int verbose=0;
static int myRank, numProcs, numMyPins;
static int cylCount, cylSize, numGlobalVertices;
static int myFirstGID, numMyVertices;
static float heavyCommWeight = 100;
static float midCommWeight = 50;
static float lowCommWeight = 1;

int *vtxGID = NULL;
int *nborIndex = NULL;
int *nborGID = NULL;
int *nborProc = NULL;
float *edgeWgt = NULL;

struct Zoltan_DD_Struct *dd=NULL;
static void debug(struct Zoltan_Struct *zz, char *s, int stop);
static void usage();

static void gid_location(int gid, int *cylID, int *ringID)
{
  div_t result;
  result = div(gid, cylSize);
  *cylID = result.quot;  /* cylinder ID */
  *ringID = result.rem;  /* ring ID */
}
static int gid_up(int gid)
{
  int cylID, ringID;
  gid_location(gid, &cylID, &ringID); 

  if (ringID == cylSize-1){
    return cylID * cylSize;
  }
  else{
    return gid + 1;
  }
}
static int gid_down(int gid)
{
  int cylID, ringID;
  gid_location(gid, &cylID, &ringID); 

  if (ringID == 0){
    return gid + cylSize - 1;
  }
  else{
    return gid - 1;
  }
}
static int gid_left(int gid)
{
  int cylID, ringID;
  gid_location(gid, &cylID, &ringID); 

  if (cylID == 0){
    return -1;
  }

  return gid - cylSize;
}
static int gid_right(int gid)
{
  int cylID, ringID;
  gid_location(gid, &cylID, &ringID); 

  if (cylID == cylCount - 1){
    return -1;
  }

  return gid + cylSize;
}
static int num_neighbors(int gid)
{
  int cylID, ringID;
  int nbors = 2;   /* up and down */

  gid_location(gid, &cylID, &ringID); 

  if (cylID > 0) nbors++;   /* left */
  if (cylID < cylCount-1) nbors++;   /* right */

  return nbors;
}

static int get_nbor_info(int gid, int *nbors, float *wgt)
{
  int i=0;
  int n;

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

static void create_a_graph()
{
  int rc, i, j, sum, n, gid, nbors[4]; 
  float wgts[4]; 
  div_t result;
  float c;

  c = sqrt((float)numGlobalVertices);

  c = (c < 3.0) ? 3.0 : c;

  cylCount = cylSize = (int)c;

  numGlobalVertices = cylCount * cylSize;

  result = div(numGlobalVertices, numProcs);

  numMyVertices = result.quot + (myRank < result.rem ? 1 : 0);

  MPI_Scan(&numMyVertices, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  myFirstGID = sum - numMyVertices;

  numMyPins = 0;

  for (i=0, j = myFirstGID; i < numMyVertices; i++, j++){
    numMyPins += num_neighbors(j);
  }

  vtxGID = (int *)malloc(sizeof(int) * numMyVertices);
  nborIndex = (int *)malloc(sizeof(int) * (numMyVertices + 1));
  nborGID = (int *)malloc(sizeof(int) * numMyPins);
  nborProc = (int *)malloc(sizeof(int) * numMyPins);
  edgeWgt = (float *)malloc(sizeof(float) * numMyPins);

  if (numMyPins && !(vtxGID || nborIndex || nborGID || nborProc || edgeWgt)){
    fprintf(stderr,"%d out of memory\n",myRank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  nborIndex[0] = 0;

  for (i=0, gid=myFirstGID, n=0; i < numMyVertices; i++, gid++){
    vtxGID[i] = gid;
    sum = get_nbor_info(gid, nbors, wgts);

    for (j=0; j < sum; j++, n++){
      nborGID[n] = nbors[j];
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

  rc = Zoltan_DD_Create(&dd, MPI_COMM_WORLD, 1, 1, 0, (int)(numGlobalVertices / numProcs), 0);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    fprintf(stderr,"%d DD Create failure\n",myRank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  rc = Zoltan_DD_Update(dd, (ZOLTAN_ID_PTR)vtxGID, NULL, NULL, NULL, numMyVertices);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    fprintf(stderr,"%d DD Update failure in create\n",myRank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  return;
}

static void reallocate_buffers(int numNewVertices, int numNewPins)
{
  int *ibuf=NULL;
  float *fbuf=NULL;

  if (numNewVertices > numMyVertices){   /* avoid realloc bug */
    ibuf = (int *)malloc(sizeof(int) * numNewVertices);
    memcpy(ibuf, vtxGID, sizeof(int) * numMyVertices);
    free(vtxGID);
    vtxGID = ibuf; 

    ibuf = (int *)malloc(sizeof(int) * (numNewVertices+1));
    memcpy(ibuf, nborIndex, sizeof(int) * (1 +numMyVertices));
    free(nborIndex);
    nborIndex = ibuf; 
  }

  if (numNewPins > numMyPins){
    ibuf = (int *)malloc(sizeof(int) * numNewPins);
    memcpy(ibuf, nborGID, sizeof(int) * numMyPins);
    free(nborGID);
    nborGID = ibuf; 

    ibuf = (int *)malloc(sizeof(int) * numNewPins);
    memcpy(ibuf, nborProc, sizeof(int) * numMyPins);
    free(nborProc);
    nborProc = ibuf; 

    fbuf = (float *)malloc(sizeof(float) * numNewPins);
    memcpy(fbuf, edgeWgt, sizeof(float) * numMyPins);
    free(edgeWgt);
    edgeWgt = fbuf; 
  }

  if (numMyPins && !(vtxGID || nborIndex || nborGID || nborProc || edgeWgt)){
    fprintf(stderr,"%d out of at realloc time\n",myRank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

static void migrate_graph(int num_exports, int num_imports, int *export_lids, int *import_gids)
{
  int i, j, k, nextv, nextp, npins, numNewVertices, numNewPins, sum;
  int nbors[4];
  float wgts[4];
  int rc;

  numNewVertices = numMyVertices - num_exports + num_imports;

  for (i=0; i < num_exports; i++){
    vtxGID[export_lids[i]] = -1;
  }

  for (i=0, nextv=0, nextp=0; i < numMyVertices; i++){
    npins = nborIndex[i+1] - nborIndex[i];
    if (vtxGID[i] >= 0){
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

  for (i=0; i < num_imports; i++){
    numNewPins += num_neighbors(import_gids[i]);
  }

  reallocate_buffers(numNewVertices, numNewPins);

  for (i=0; i < num_imports; i++, nextv++){
    vtxGID[nextv] = import_gids[i];
    sum = get_nbor_info(import_gids[i], nbors, wgts);

    for (j=0; j < sum; j++, nextp++){
      nborGID[nextp] = nbors[j];
      edgeWgt[nextp] = wgts[j];
    }
    nborIndex[nextv+1] = nborIndex[nextv] + sum; 
  }

  numMyVertices = numNewVertices;
  numMyPins = numNewPins;

  rc = Zoltan_DD_Update(dd, (ZOLTAN_ID_PTR)vtxGID, NULL, NULL, NULL, numMyVertices);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    fprintf(stderr,"%d DD Update failure in migrate\n",myRank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  rc = Zoltan_DD_Find(dd, (ZOLTAN_ID_PTR)nborGID, NULL, NULL, NULL, numNewPins, nborProc);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    fprintf(stderr,"%d DD find failure\n",myRank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

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
  int rc, do_hier;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  char *platform=NULL, *node_topology=NULL, *machine_topology=NULL;
  char *graph_package=NULL;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  struct option opts[10];
  double comm_times[10];
  int nvert=0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  /* options */

  opts[0].name = "platform";
  opts[0].has_arg = 1;
  opts[0].flag = NULL;
  opts[0].val = 1;

  opts[1].name = "node_topology";
  opts[1].has_arg = 1;
  opts[1].flag = NULL;
  opts[1].val = 2;

  opts[2].name = "machine_topology";
  opts[2].has_arg = 1;
  opts[2].flag = NULL;
  opts[2].val = 3;

  opts[3].name = "size";
  opts[3].has_arg = 1;
  opts[3].flag = NULL;
  opts[3].val = 4;

  opts[4].name = "verbose";
  opts[4].has_arg = 0;
  opts[4].flag = NULL;
  opts[4].val = 5;

  opts[5].name = "help";
  opts[5].has_arg = 0;
  opts[5].flag = NULL;
  opts[5].val = 6;

  opts[6].name = "graph_package";
  opts[6].has_arg = 1;
  opts[6].flag = NULL;
  opts[6].val = 7;

  opts[7].name = 0;
  opts[7].has_arg = 0;
  opts[7].flag = NULL;
  opts[7].val = 0;

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
      node_topology = optarg;
      if (myRank == 0)
        printf( "For node topology %s\n",optarg);
    }
    else if (rc == 3){
      machine_topology = optarg;
      if (myRank == 0)
        printf( "For machine with %s cpus\n",optarg);
    }
    else if (rc == 7){
      graph_package = optarg;
      if (myRank == 0)
        printf( "Zoltan parameter GRAPH_PACKAGE = %s\n",graph_package);
    }
    else if (rc == 4){
      nvert = atoi(optarg);
      if (nvert < 1){
        if (myRank == 0)
          fprintf(stderr,"ERROR: --size={approximate number of vertices}\n");
        MPI_Finalize();
        exit(1);
      }
      else{
        if (myRank == 0){
          printf( "Graph will have approximately %d vertices.\n",nvert);
        }
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
    else if (rc <= 0){
      break;
    }
  }

  if ((platform==NULL) && (node_topology==NULL) && (machine_topology==NULL)){
    if (myRank == 0)
      fprintf(stderr,"No platform or topology, so we'll skip hierarchical partitioning\n");
    do_hier = 0;
  }
  else{
    do_hier = 1;
  }

  /* start */

  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    if (myRank == 0) printf("sorry...\n");
    MPI_Finalize();
    exit(0);
  }

  Zoltan_Memory_Debug(0);

  if (nvert > 0)
    numGlobalVertices = nvert;
  else
    numGlobalVertices = NUM_GLOBAL_VERTICES;

  create_a_graph(); 
    
  time_communication(comm_times+0);      /* Without partitioning */

  zz = Zoltan_Create(MPI_COMM_WORLD);
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
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

  migrate_graph(numExport, numImport, (int *)exportLocalGids, (int *)importGlobalGids);

  if (verbose){
    debug(zz, "After flat partitioning and migration", 0);
  }

  time_communication(comm_times+1);      /* With graph partitioning */

  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  if (do_hier){

    /* HIERARCHICAL PARTITION */
  
    Zoltan_Set_Param(zz, "LB_METHOD", "HIER");
    Zoltan_Set_Param(zz, "HIER_DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "HIER_ASSIST", "1");

    /* TODO: Suppose graph is not symmetric, and we request SYMMETRIZE.  Do we still get
     *  a "good" answer when each sub-graph in the hierarchy is symmetrized?
     */

    if (machine_topology)
      Zoltan_Set_Param(zz, "MACHINE_TOPOLOGY", machine_topology);
    else if (node_topology)
      Zoltan_Set_Param(zz, "NODE_TOPOLOGY", node_topology);
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
  
    if (rc != ZOLTAN_OK){
      printf("sorry...\n");
      MPI_Finalize();
      Zoltan_Destroy(&zz);
      exit(0);
    }
  
    migrate_graph(numExport, numImport, (int *)exportLocalGids, (int *)importGlobalGids);
  
    if (verbose){
      debug(zz, "After hierarchical partitioning and migration", 0);
    }
  
    time_communication(comm_times+2);      /* With hierarchical graph partitioning */
  
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                        &exportProcs, &exportToPart);

  }

  Zoltan_Destroy(&zz);

  if (dd) Zoltan_DD_Destroy(&dd);

  if (vtxGID) free(vtxGID);
  if (nborIndex) free(nborIndex);
  if (nborGID) free(nborGID);
  if (nborProc) free(nborProc);
  if (edgeWgt) free(edgeWgt);

  MPI_Finalize();

  return 0;
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
          fprintf(stdout,"%d global vertices\n",numGlobalVertices); 
        }
  
        fprintf(stdout,"Partition %d, %d vertices:\n",p,numMyVertices);
        for (i=0, k=0; i < numMyVertices; i++){
          fprintf(stdout,"%d: ",vtxGID[i]);
          nedges = nborIndex[i+1] - nborIndex[i];
           
          for (j=0; j < nedges; j++,k++){
            fprintf(stdout,"%d/%f/%d ",nborGID[k],edgeWgt[k],nborProc[k]);
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

extern zoltan_platform_specification zoltan_hier_platform_specs[ZOLTAN_HIER_LAST_PLATFORM];

static void usage()
{
  int i;
  printf( "\nUsage: --verbose\n");
  printf( "\n       --graph_package={parmetis|scotch|phg}\n");
  printf( "\n       --platform=desc | --node_topology=desc | --machine_topology=desc\n");
  printf( "\n       --size={approximate global number of vertices}\n");

  printf( "\nPlatform Names that Zoltan knows about:");

  for (i=0; i < ZOLTAN_HIER_LAST_PLATFORM; i++){
    if (i%8 == 0) printf("\n");
    printf( "  %s", zoltan_hier_platform_specs[i].platform_name );
  }
  printf( "\n\nA node_topology description is a list of integers, for example\n");
  printf( "  Dual socket, quad core: 2, 4\n");
  printf( "  Quad socket, six cores with core pairs sharing a cache: 4, 3, 2\n");
  printf( "A machine_topology description is a list of integers (usually just 1), for example\n");
  printf( "  Dual core workstation: 2\n\n");

  printf( "The default global number of vertices is %d\n",NUM_GLOBAL_VERTICES);
}

