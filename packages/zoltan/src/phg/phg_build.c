/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"
#include "zz_const.h"
#include "third_library_const.h"
#include "third_library_tools.h"
#include "zz_util_const.h"

/*#define DEBUG_FILL_HYPERGRAPH 1*/
    
#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}
#define FATAL_ERROR(s) { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, s); \
  ierr = ZOLTAN_FATAL; \
  goto End; \
}
#define CHECK_FOR_MPI_ERROR(rc)  \
  if (rc != MPI_SUCCESS){ \
    MPI_Error_string(rc, mpi_err_str, &mpi_err_len);  \
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, mpi_err_str); \
    ierr = ZOLTAN_FATAL; \
    goto End; \
  }

static char mpi_err_str[MPI_MAX_ERROR_STRING];
static int mpi_err_len;

  
/* Macro returning the global number of the first repartition vertex or 
 * repartition edge to be added to a processor row or column, respectively. */
#define FirstRepart(myProc, nProc, Gn) \
        ((myProc) * (int)((Gn) / (nProc)) + MIN((myProc), (Gn) % (nProc)))

/* Macro returning the number of repartition vertices or repartition edges 
 * to be added to a processor row or column, respectively. */
#define NumRepart(myProc, nProc, Gn) \
        ((myProc >= 0) ? (int)((Gn) / (nProc)) + ((myProc) < ((Gn) % (nProc))) \
                       : 0)

/* Macro returning the processor row or column storing a given 
 * repartition vertex or edge, respectively. */
#define ProcForRepart(i, repart_dist, nProc) \
        Zoltan_PHG_Gno_To_Proc_Block(i, repart_dist, nProc)

/*****************************************************************************/
/* Function prototypes */

static int hash_lookup(ZZ*, ZOLTAN_ID_PTR, int, struct Hash_Node**, struct Hash_Node **);
#ifdef DEBUG_FILL_HYPERGRAPH
static void print_hypergraph(ZZ *zz, ZHG *zhg, int sumWeight);  /* for debugging */
#endif

static int Zoltan_PHG_Add_Repart_Data(ZZ *, ZHG *, HGraph *, int *,
                                      PHGPartParams *, Partition);
static int removed_cuts_local(ZZ *zz, ZHG *zhg, 
                int max_parts, int *pin_parts, double *loccuts);
static int removed_cuts_global(ZZ *zz, ZHG *zhg, 
                int max_parts, int *pin_parts, double *loccuts, int tag);
static int getObjectSizes(ZZ *zz, ZHG *zhg);
    
/*****************************************************************************/
int Zoltan_PHG_Build_Hypergraph(
  ZZ *zz,                            /* Input : Zoltan data structure */
  ZHG **zoltan_hg,                   /* Output: Hypergraph to be allocated and built.*/
  Partition *input_parts,            /* Output: Initial partition assignments for
                                         vtxs (in 2D distribution); length = 
                                         zoltan_hg->HG->nVtx.  */
  PHGPartParams *hgp                 /* Input : Parameters for PHG partitioning.*/
)
{
/* allocates and builds hypergraph data structure using callback routines */ 
ZHG *zhg;                     /* Temporary pointer to Zoltan_HGraph. */
HGraph *phgraph;             /* Temporary pointer to HG field */
int ierr = ZOLTAN_OK;
char *yo = "Zoltan_PHG_Build_Hypergraph";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Allocate a Zoltan hypergraph.  */
  zhg = *zoltan_hg = (ZHG*) ZOLTAN_MALLOC (sizeof(ZHG));
  if (zhg == NULL) MEMORY_ERROR;

  /* Initialize the Zoltan hypergraph data fields. */
  zhg->GIDs = NULL;
  zhg->LIDs = NULL;
  zhg->Input_Parts = NULL;
  zhg->Output_Parts = NULL;
  zhg->AppObjSizes = NULL;
  zhg->showMoveVol = 0;
  zhg->GnRepartVtx = 0;
  zhg->GnRepartEdge = 0;
  zhg->nObj = 0;
  zhg->GnObj = 0;
  zhg->nRemove = 0;
  zhg->Remove_EGIDs = NULL;
  zhg->Remove_ELIDs = NULL;
  zhg->Remove_Esize = NULL;
  zhg->Remove_GEsize = NULL;
  zhg->Remove_Ewgt = NULL;
  zhg->nRemovePins = 0;
  zhg->Remove_Pin_GIDs = NULL;
  zhg->Remove_Pin_Procs = NULL;
  zhg->VtxPlan = NULL;
  zhg->Recv_GNOs = NULL;
  zhg->nRecv_GNOs = 0;

  phgraph = &(zhg->HG);
  Zoltan_HG_HGraph_Init(phgraph);

  /* just set the pointer of phgraph's comm to hgp's comm */
  phgraph->comm = &hgp->globalcomm;

  /* Use callback functions to build the hypergraph. */

  ierr = Zoltan_PHG_Fill_Hypergraph(zz, zhg, hgp, input_parts);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph");
    goto End;
  }

#ifdef KDDKDD_NO_COORDINATES_FOR_NOW
  /* KDDKDD DON'T KNOW WHAT TO DO WITH THE COORDINATES WITH 2D DISTRIB ANYWAY */
  if (zz->Get_Num_Geom != NULL && 
      (zz->Get_Geom != NULL || zz->Get_Geom_Multi != NULL)) {
     /* Geometric callbacks are registered;       */
     /* get coordinates for hypergraph objects.   */
     ZOLTAN_TRACE_DETAIL(zz, yo, "Getting Coordinates.");
     ierr = Zoltan_Get_Coordinates(zz, phgraph->nVtx, zhg->GIDs,
      zhg->LIDs, &(phgraph->nDim), &(phgraph->coor));
  }
#endif

  if (hgp->check_graph) {
    ierr = Zoltan_HG_Check(zz, phgraph);
    if (ierr == ZOLTAN_WARN) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning returned from Zoltan_HG_Check");
    }
    else if (ierr != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_HG_Check");
      goto End;     
    }
  }


  if (hgp->output_level >= PHG_DEBUG_PLOT)
    Zoltan_PHG_Plot_2D_Distrib(zz, &(zhg->HG));

  if (hgp->output_level >= PHG_DEBUG_PRINT)
    Zoltan_HG_HGraph_Print(zz, zhg, &(zhg->HG), *input_parts, stdout);

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    /* Return NULL zhg */
    Zoltan_PHG_Free_Hypergraph_Data(zhg);
    ZOLTAN_FREE(&zhg);
    *zoltan_hg = NULL;
  }
    
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/
/* Zoltan_PHG_Fill_Hypergraph                                                */
/*****************************************************************************/

/* 
 * Structures to hold hypergraph data returned by query functions,
 * and hypergraph data gathered by processes to which edges/vertices
 * map to via a hash function.
 *
 * TODO - special case where one process returns all edges, all
 * vertices.  Dispense with the hashing assignment in that case.
 * TODO - streamline code for single process application
 */
 
typedef struct _myObj{  /* Vertices returned in Get_Obj_List queries */
  int    size;          /* # objects (vertices) on proc incl lone verts */
  int GnVtx;            /* number of objects across all processes */
  ZOLTAN_ID_PTR vtxGID; /* Global ID of each vertex */
  int    *vtx_gno;      /* Global numbers in range [0,GnVtx-1] */
  float  *vwgt;         /* Vertex weights for nVtx on processor objects */
  int    *fixed;        /* Vertex assignments for fixed vertices  */ 
  int    *vtxHash;      /* Process to which GID hashes, temporary owner */
  int    *numHedges;    /* Number of hyperedges containing vertex       */
  int    *numAllHedges; /* Including removed hyperedges */
}zoltan_objects;

typedef struct _myPin{      /* Pins returned by hypergraph query functions */
  int           nHedges;    /* number of (partial) hyperedges */
  ZOLTAN_ID_PTR edgeGID;    /* edge global IDs */
  int           *esizes;    /* local size in pins of each hyperedge */
  ZOLTAN_ID_PTR pinGID;     /* global ID of pin vertex */
  int           numPins;    /* sum of esizes array */

  int           *vtxHash;   /* Temp hashed owner of pinGID */
  int           *pinGNO;    /* global number (0, GnVtx-1) for pinGID */
  int           *pinProc;   /* Real owner of pinGID (ret'd it in Get_Obj_List)*/

  int           *edgeHash;  /* process assigned edgeGID by hash function */
  float         *ewgt;      /* weights for edges (nHedges * weight_dim)  */
  int           *edgeGNO;   /* edge global numbers, consecutive, 0-based */
  int           *edgeGSize; /* global number of pins in hyperedge        */

  ZOLTAN_ID_PTR eGIDs;  /* save edge comm plan, it's used more than once*/
  int           *hashIdx;
  int           numE;
  ZOLTAN_COMM_OBJ *ePlan;
}zoltan_pins;

typedef struct _myEW{     /* Values returned by edge weight query functions */
  int           size;       /* number of edges */
  ZOLTAN_ID_PTR edgeGID;   /* edge global IDs */
  int           *edgeHash;  /* process assigned this edge by hash function */
  float         *wgt;       /* weights supplied by query function for edge */
}zoltan_ews;

typedef struct _hshEdge{ /* Edges assigned to this process with hash func */
  int           size;        /* number of used edges assigned to this process */
  int           GnEdge;      /* Total size across all procs. */
  ZOLTAN_ID_PTR edgeGID;    /* edge global IDs  */
  int           *edgeGNO;    /* edge global numbers  */
  int           *numPins;    /* global number of pins for this edge  */
  float         *wgt;        /* weights for these edges  */
}zoltan_temp_edges;

typedef struct _hshVtx{ /* Vertices assigned to this process with hash func */
  int           size;      /* number of vertices assigned to this process */
  ZOLTAN_ID_PTR vtxGID;   /* vertex global IDs  */
  int           *vtxOwner; /* process that returned vtx in Get_Obj_List  */
}zoltan_temp_vertices;

/* 
 * A search structure, to find the index of a global ID in any of the
 * above structures.
 */

typedef struct _GID_lookup{
  struct Hash_Node *htTop;
  struct Hash_Node **ht;
  int table_size;
  int numGIDs;
  int lenGID;
}GID_lookup;

static void free_zoltan_objects(zoltan_objects *zo);
static void free_zoltan_pins(zoltan_pins *zp);
/*static void print_zoltan_pins(zoltan_pins *z, int me, int ewgt_dim);*/
static void free_zoltan_ews(zoltan_ews *zew);
static void free_zoltan_temp_edges(zoltan_temp_edges *zte);
static void free_zoltan_temp_vertices(zoltan_temp_vertices *ztv);
static int map_GIDs_to_processes(ZZ *zz, ZOLTAN_ID_PTR eid, int size,
  int lenGID, int **hashedProc, int nprocs);
static GID_lookup *create_GID_lookup_table(ZOLTAN_ID_PTR gids, 
  int size, int lenGID);
static GID_lookup *create_GID_lookup_table2(ZOLTAN_ID_PTR gids, 
  int ngids, int lenGID);
static int lookup_GID(GID_lookup *lu, ZOLTAN_ID_PTR gid);
static void free_GID_lookup_table(GID_lookup **lu);

int Zoltan_PHG_Fill_Hypergraph(
  ZZ *zz,        /* Input : Zoltan data structure */
  ZHG *zhg,      /* Output: Description of hypergraph provided by the application. */
  PHGPartParams *hgp,      /* Input : Parameters for PHG partitioning.*/
  Partition *input_parts   /* Output: Initial partition assignment of vtxs in 
                              2D data distribution; length = zhg->HG->nVtx. */
)
{
/* Routine to call HG query function and build HG data structure. 
 * Output is a fully functioning parallel hypergraph with 2D distribution of
 * pins (non-zeros).
 */

char *yo = "Zoltan_PHG_Fill_Hypergraph";

ZOLTAN_COMM_OBJ *plan=NULL;

int i, j, w, cnt, dim, rc;
int msg_tag = 30000;
int ierr = ZOLTAN_OK;
int nProc = zz->Num_Proc;
int nRequests;
ZOLTAN_ID_PTR pin_requests = NULL;
ZOLTAN_ID_PTR gid_requests = NULL;
int *pin_info = NULL;
int *gid_info = NULL;
float *gid_weights = NULL;
float *src, *dest;
char *have_wgt = NULL;
int edge_gno, edge_Proc_y;
int vtx_gno, vtx_Proc_x;
int nnz, idx, method_repart;
int *proclist = NULL;
int *sendbuf = NULL;
int *egno = NULL;
int *pinIdx = NULL;
int *keep_edge = NULL;
int *gtotal = NULL;  
int *mycnt = NULL;
int *gcnt = NULL;
int *nonzeros = NULL;
int *tmparray = NULL;
int *hindex = NULL, *hvertex = NULL;
int *dist_x = NULL, *dist_y = NULL;
int nEdge, nVtx, nwgt = 0;
int nrecv, *recv_gno = NULL; 
int *tmpparts = NULL;
int add_vweight;
int ew_dim = zz->Edge_Weight_Dim;
int ew_op;
int hypergraph_callbacks = 0;
int graph_callbacks = 0;
int use_all_neighbors = 1;
ZOLTAN_ID_PTR global_ids, ew_lids = NULL;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int changed, next_id;
HGraph *phg = &(zhg->HG);
int nProc_x, nProc_y, myProc_x, myProc_y;
int proc_offset;

float frac_x, frac_y;
float *tmpwgts=NULL; 
float *fromwgt, *towgt, *wgts, *calcVwgt = NULL;

float edgeSizeThreshold;
int randomizeInitDist, zoltan_lb_eval;
int final_output, add_obj_weight, removedEdges=0;
PHGPartParams *temphgp = NULL;

ZOLTAN_ID_PTR elids = NULL;

zoltan_objects       myObjs;
GID_lookup           *lookup_myObjs = NULL;
zoltan_pins          myPins;
zoltan_ews           myEWs;
zoltan_temp_edges    myHshEdges;
GID_lookup           *lookup_myHshEdges = NULL;
zoltan_temp_vertices myHshVtxs;
GID_lookup           *lookup_myHshVtxs = NULL;

int GnFixed=0, nFixed=0;                                     
int *tmpfixed = NULL, *fixedPart = NULL;              
ZOLTAN_ID_PTR fixedGIDs = NULL;

int nRepartEdge = 0, nRepartVtx = 0;


  ZOLTAN_TRACE_ENTER(zz, yo);

  memset(&myObjs, 0, sizeof(zoltan_objects));
  memset(&myPins, 0, sizeof(zoltan_pins));
  memset(&myEWs, 0, sizeof(zoltan_ews));
  memset(&myHshEdges, 0, sizeof(zoltan_temp_edges));
  memset(&myHshVtxs, 0, sizeof(zoltan_temp_vertices));

  /**************************************************/
  /* Determine parameters                           */
  /**************************************************/

  if (hgp){
    randomizeInitDist = hgp->RandomizeInitDist;
    ew_op = hgp->edge_weight_op;
    edgeSizeThreshold = hgp->EdgeSizeThreshold;
    final_output = hgp->final_output;
    add_obj_weight = hgp->add_obj_weight;
    if ((hgp->convert_str[0] == 'n') ||   /* "neighbors" */
        (hgp->convert_str[0] == 'N')){
      use_all_neighbors = 1;
    }
    else{                                 /* "pairs"     */
      use_all_neighbors = 0;
    }
    zoltan_lb_eval = 0;

    method_repart = (!strcasecmp(hgp->hgraph_method, "REPARTITION"));
  }
  else{
    /*   
     * hgp and input_parts are undefined when we are called from 
     * Zoltan_LB_Eval.  (This only happens in the pins callback case.)
     * In this case, we want all pins to be in the removed list.
     * And we don't care how global numbers are assigned.
     */
    temphgp = (PHGPartParams *)ZOLTAN_MALLOC(sizeof(PHGPartParams));
    if (!temphgp) MEMORY_ERROR;
    
    Zoltan_PHG_Initialize_Params(zz, NULL, temphgp);

    randomizeInitDist = 0;   /* faster */
    ew_op = temphgp->edge_weight_op;
    edgeSizeThreshold = 0;   /* place all edges in "removed" list */
    final_output = 1;        /* yes, compile a "removed" list     */
    add_obj_weight = temphgp->add_obj_weight;
    if ((temphgp->convert_str[0] == 'n') ||   /* "neighbors" */
        (temphgp->convert_str[0] == 'N')){
      use_all_neighbors = 1;
    }
    else{                                 /* "pairs"     */
      use_all_neighbors = 0;
    }
    method_repart = (!strcasecmp(temphgp->hgraph_method, "REPARTITION"));

    if (temphgp->globalcomm.row_comm != MPI_COMM_NULL)
      MPI_Comm_free(&(temphgp->globalcomm.row_comm));
    if (temphgp->globalcomm.col_comm != MPI_COMM_NULL)
      MPI_Comm_free(&(temphgp->globalcomm.col_comm));
    if (temphgp->globalcomm.Communicator != MPI_COMM_NULL)
      MPI_Comm_free(&(temphgp->globalcomm.Communicator));

    ZOLTAN_FREE(&temphgp);

    zoltan_lb_eval = 1;
  }
  if (zz->Get_HG_Size_CS && zz->Get_HG_CS){
    hypergraph_callbacks = 1;
  }
  if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
           (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) {
    graph_callbacks = 1;
  }

  /* Pick right callbacks if both are available. */
  if (graph_callbacks && hypergraph_callbacks){
    if (zz->LB.Method == GRAPH)
      hypergraph_callbacks = 0;
  }

  /**************************************************/
  /* Obtain vertex information from the application */
  /**************************************************/

  ierr = Zoltan_Get_Obj_List(zz, &(zhg->nObj), &(zhg->GIDs), &(zhg->LIDs), 
                             zz->Obj_Weight_Dim, &myObjs.vwgt,
                             &(zhg->Input_Parts));

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
    goto End;
  }

  MPI_Allreduce(&(zhg->nObj), &i, 1, MPI_INT, MPI_MAX, zz->Communicator);

  if (i < 1){
    if (zz->Proc == 0){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "No objects to partition on any process")
    }
    goto End;
  }

  myObjs.size = zhg->nObj;
  myObjs.vtxGID = zhg->GIDs;

  if (myObjs.size){
    myObjs.vtx_gno = (int *)ZOLTAN_MALLOC(sizeof(int) * myObjs.size);
    myObjs.numHedges = (int *)ZOLTAN_MALLOC(sizeof(int) * myObjs.size);
    myObjs.numAllHedges = (int *)ZOLTAN_MALLOC(sizeof(int) * myObjs.size);

    if (!myObjs.vtx_gno || !myObjs.numHedges || !myObjs.numAllHedges){
      MEMORY_ERROR;
    }
  }

  /*
   * Create a search structure to lookup my vertex information
   */

  lookup_myObjs = create_GID_lookup_table(myObjs.vtxGID,
                       myObjs.size, num_gid_entries);

  if (!lookup_myObjs) MEMORY_ERROR;
  

  if (zz->Get_Num_Fixed_Obj) {  /* If registered query fixed objects/vertices */
    nFixed = zz->Get_Num_Fixed_Obj (zz->Get_Num_Fixed_Obj_Data, &ierr);
    if (ierr != ZOLTAN_OK) {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Error getting Number of Fixed Objs");
       goto End;
    }
    MPI_Allreduce (&nFixed, &GnFixed, 1, MPI_INT, MPI_SUM, zz->Communicator);
    if (hgp) hgp->UseFixedVtx = GnFixed;  /* Don't need to set UseFixedVtx
                                             if called from Zoltan_LB_Eval. */
    
    if (GnFixed && zhg->nObj) {
      myObjs.fixed = (int*) ZOLTAN_MALLOC (sizeof(int) * zhg->nObj);
      if (!myObjs.fixed) MEMORY_ERROR;
      for (i = 0; i < zhg->nObj; i++)
        myObjs.fixed[i] = -1;              /* default - no fixed assignment */
    }
      
    if (GnFixed && nFixed && zhg->nObj)  {
      fixedPart = (int*) ZOLTAN_MALLOC (sizeof(int) * nFixed);
      fixedGIDs    = ZOLTAN_MALLOC_GID_ARRAY (zz, nFixed);
       
      if (!fixedPart || !fixedGIDs)
        MEMORY_ERROR;
          
      if (zz->Get_Fixed_Obj_List) {
        zz->Get_Fixed_Obj_List (zz->Get_Fixed_Obj_List_Data, nFixed,
         num_gid_entries, fixedGIDs, fixedPart, &ierr);
        if (ierr != ZOLTAN_OK) {
          ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Error getting Fixed Obj List");
          goto End;
          }         
             
        for (i = 0; i < nFixed; i++) {
          j = lookup_GID (lookup_myObjs, fixedGIDs + (i * num_gid_entries));
          if (j < 0)
            FATAL_ERROR ("Unexpected Fixed Vertex GID received");
          myObjs.fixed[j] = fixedPart[i];  /* overwrite fixed assignment */
        }
      }
      ZOLTAN_FREE(&fixedGIDs);
      ZOLTAN_FREE(&fixedPart);
    }
  }     

  if (hypergraph_callbacks){
    /*
     * Processes will be returning vertex GIDs in query functions.  But
     * they don't know which process owns those vertices (i.e. returned
     * them in Get_Obj_List).
     *
     * Use a hash function to assign vertex global IDs to processes.
     * The assigned process will learn the owner of the vertex, and
     * return that information when requested by other processes.
     */
    ierr = map_GIDs_to_processes(zz, zhg->GIDs, zhg->nObj, num_gid_entries,
                             &myObjs.vtxHash, nProc);

    if (ierr != ZOLTAN_OK){
      goto End;
    }
  
    /* 
     * Use an unstructured communication plan to send vertex global
     * IDs to their assigned process.
     */
  
    ierr = Zoltan_Comm_Create(&plan, myObjs.size, myObjs.vtxHash, 
             zz->Communicator, msg_tag, &myHshVtxs.size);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  
    if (myHshVtxs.size > 0){
      myHshVtxs.vtxGID = ZOLTAN_MALLOC_GID_ARRAY(zz, myHshVtxs.size);
      myHshVtxs.vtxOwner = (int *)ZOLTAN_MALLOC(sizeof(int) * myHshVtxs.size);

      if (!myHshVtxs.vtxGID || !myHshVtxs.vtxOwner) MEMORY_ERROR;
    }
  
    msg_tag--;
  
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myObjs.vtxGID, 
           sizeof(ZOLTAN_ID_TYPE) * num_gid_entries, (char *)myHshVtxs.vtxGID);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  
    ierr = Zoltan_Comm_Info(plan, 
            NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, myHshVtxs.vtxOwner,
            NULL);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  
    Zoltan_Comm_Destroy(&plan);
  
    /*
     * Create a search structure allowing me to look up the owner of a given 
     * vertex GID upon request.
     */
  
    lookup_myHshVtxs = create_GID_lookup_table(myHshVtxs.vtxGID,
                         myHshVtxs.size, num_gid_entries);

    if (!lookup_myHshVtxs) MEMORY_ERROR;
  }

  /*******************************************************************/
  /* Assign consecutive numbers (gnos) based on the order of the ids */
  /*******************************************************************/
  global_ids = zhg->GIDs;

  if (randomizeInitDist) { 
    /* Randomize the input vertices */
    int tmp;
    gtotal = (int *) ZOLTAN_CALLOC(3*zz->Num_Proc+1, sizeof(int));
    if (!gtotal) MEMORY_ERROR;

    mycnt  = gtotal + zz->Num_Proc + 1;
    gcnt   = mycnt + zz->Num_Proc;

    /* Compute random processor bin. */
    /* Temporarily store processor bin number in myObjs.vtx_gno. */
    /* Count how many local vtxs selected processor bin */
    Zoltan_Srand(Zoltan_Rand(NULL)+zz->Proc, NULL);
    for (i = 0; i < myObjs.size; i++) {
      myObjs.vtx_gno[i] = Zoltan_Rand_InRange(NULL, zz->Num_Proc);
      mycnt[myObjs.vtx_gno[i]]++;
    }
    /* Compute prefix of mycnt */
    rc = MPI_Scan(mycnt, gcnt, zz->Num_Proc, MPI_INT, MPI_SUM,
                  zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    rc = MPI_Allreduce(mycnt, gtotal, zz->Num_Proc, MPI_INT, MPI_SUM, 
                  zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Compute first gno for vertices going to each target bin */
    for (tmp = 0, i = 0; i < zz->Num_Proc; i++) {
      gcnt[i] -= mycnt[i];
      tmp += gtotal[i];
      gtotal[i] = tmp - gtotal[i];
    }
    zhg->GnObj = myObjs.GnVtx = gtotal[zz->Num_Proc] = tmp;

    /* Assign gnos sequential from gcnt[bin]. */
    for (i=0; i< myObjs.size; i++) {
      tmp = myObjs.vtx_gno[i];
      myObjs.vtx_gno[i] = gtotal[tmp] + gcnt[tmp];
      gcnt[tmp]++;
    }
  }
  else {
    /* Linearly order the input vertices */
    gtotal = (int *) ZOLTAN_MALLOC((nProc+1) * sizeof(int));
    if (!gtotal) MEMORY_ERROR;

    /* Construct gtotal[i] = the number of vertices on all procs < i. */
    /* Scan to compute partial sums of the number of objs */

    rc = MPI_Scan(&myObjs.size, gtotal, 1, MPI_INT, MPI_SUM, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Gather data from all procs */

    rc = MPI_Allgather (&(gtotal[0]), 1, MPI_INT,
                   &(gtotal[1]), 1, MPI_INT, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)
    gtotal[0] = 0;
    zhg->GnObj = myObjs.GnVtx = gtotal[nProc];

    for (i=0; i< myObjs.size; i++) {
      myObjs.vtx_gno[i] =  gtotal[zz->Proc]+i;
    }
  }

  /***********************************************************************/
  /* Get hyperedge information from application through query functions. */
  /***********************************************************************/

  if (hypergraph_callbacks){
    /*
     * Each processor:
     *   owns a set of pins (nonzeros)
     *   may provide some edge weights
     *
     * The edge weights supplied by the application may not be
     * for the edges represented by it's pins.  More than one
     * process may provide edge weights for the same edge.  We combine
     * them (or flag an error) according to the setting of the
     * PHG_EDGE_WEIGHT_OPERATION parameter.
     *
     * We assume than no two processes will supply the same pin.  
     * But more than one process may supply pins for the same edge.
     */

    /*
     * Call the compressed storage pin query functions. 
     */

    ierr = Zoltan_Call_Hypergraph_Pin_Query(zz, &myPins.nHedges,
                 &myPins.numPins, &myPins.edgeGID, &pinIdx,
                 &myPins.pinGID);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    myPins.esizes = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.nHedges);

    if (myPins.nHedges && !myPins.esizes){
      ZOLTAN_FREE(&pinIdx);
      MEMORY_ERROR;
    }
    
    for (i=0; i<myPins.nHedges; i++){
      myPins.esizes[i] = pinIdx[i+1] - pinIdx[i];
    }

    ZOLTAN_FREE(&pinIdx);

    /*
     * Determine the process owning each of my pin vertices.
     */

    ierr = map_GIDs_to_processes(zz, myPins.pinGID, myPins.numPins,
               num_gid_entries, &myPins.vtxHash, nProc);


    pin_requests = NULL;
    pin_info = NULL;

    ierr = Zoltan_Comm_Create(&plan, myPins.numPins, myPins.vtxHash, 
           zz->Communicator, msg_tag, &nRequests);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    if (nRequests > 0){
      pin_requests = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
      pin_info = (int *)ZOLTAN_MALLOC(sizeof(int) * nRequests);
      if (!pin_requests || !pin_info) MEMORY_ERROR;
    }
  
    msg_tag--;
  
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myPins.pinGID, 
             sizeof(ZOLTAN_ID_TYPE) * num_gid_entries, (char *)pin_requests);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  
    for (i=0; i<nRequests; i++){
      j = lookup_GID(lookup_myHshVtxs, pin_requests + ( i * num_gid_entries));

      if (j < 0) FATAL_ERROR("Unexpected vertex GID received");

      pin_info[i] = myHshVtxs.vtxOwner[j];
    }

    if (myPins.numPins > 0){
      myPins.pinProc = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.numPins);
      if (!myPins.pinProc) MEMORY_ERROR;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *)pin_info, sizeof(int), 
                  NULL, (char *)myPins.pinProc);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    Zoltan_Comm_Destroy(&plan);

    ZOLTAN_FREE(&pin_requests);
    ZOLTAN_FREE(&pin_info);

    /* We're done with these: */

    free_GID_lookup_table(&lookup_myHshVtxs);
    free_zoltan_temp_vertices(&myHshVtxs);

    /* 
     * Use a hash function to assign each edge GID to a process.  
     */

    ierr = map_GIDs_to_processes(zz, myPins.edgeGID, myPins.nHedges,
               num_gid_entries, &myPins.edgeHash, nProc);

    if (ierr != ZOLTAN_OK){
      goto End;
    }

    /* 
     * Use an unstructured communication plan to send edge global
     * IDs to their assigned process.
     */
  
    myPins.eGIDs = NULL;
    global_ids = NULL;
    msg_tag--;
    ierr = Zoltan_Comm_Create(&myPins.ePlan, myPins.nHedges, myPins.edgeHash, 
             zz->Communicator, msg_tag, &myPins.numE);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  
    if (myPins.numE > 0){
      myPins.eGIDs = ZOLTAN_MALLOC_GID_ARRAY(zz, myPins.numE);
      myPins.hashIdx = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.numE);

      if (!myPins.eGIDs || !myPins.hashIdx) MEMORY_ERROR;
    }
  
    msg_tag--;
    ierr = Zoltan_Comm_Do(myPins.ePlan, msg_tag, (char *)myPins.edgeGID, 
             sizeof(ZOLTAN_ID_TYPE) * num_gid_entries, 
             (char *)myPins.eGIDs);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    if (myPins.numE > 0){
      global_ids = ZOLTAN_MALLOC_GID_ARRAY(zz, myPins.numE);
      if (!global_ids) MEMORY_ERROR;
      memcpy(global_ids, myPins.eGIDs, 
         myPins.numE * num_gid_entries * sizeof(ZOLTAN_ID_TYPE));
    }

    /*
     * Create a search structure allowing me to look up info about
     * any of the edge GIDs that were just assigned to me.  
     * Rewrites global_ids with a list of unique GIDs.
     */
  
    lookup_myHshEdges = create_GID_lookup_table2(global_ids,
                         myPins.numE, num_gid_entries);

    myHshEdges.size = lookup_myHshEdges->numGIDs;
    myHshEdges.edgeGID = global_ids; 

    if (myHshEdges.size > 0){
      myHshEdges.edgeGNO = (int *)ZOLTAN_MALLOC(myHshEdges.size * sizeof(int));
      myHshEdges.numPins = (int *)ZOLTAN_CALLOC(myHshEdges.size , sizeof(int));
      if (ew_dim){
        myHshEdges.wgt = 
          (float *)ZOLTAN_MALLOC(myHshEdges.size * sizeof(float) * ew_dim);
      }

      if (!myHshEdges.edgeGNO || !myHshEdges.numPins ||
          (ew_dim && !myHshEdges.wgt)){
        MEMORY_ERROR;
      }
    }

    /* Process to which edge is assigned calculates # of pins in the edge */
    if (myPins.nHedges && num_lid_entries)  {
      elids = ZOLTAN_MALLOC_LID_ARRAY(zz, myPins.nHedges);

      if (!elids)  MEMORY_ERROR; 
    
      for (i=0; i<myPins.nHedges; i++){
        elids[i*num_lid_entries] = i;
      }
    }
    if (myPins.numE){
      pin_info = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.numE);
      if (!pin_info) MEMORY_ERROR;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do(myPins.ePlan, msg_tag, (char *)myPins.esizes, 
             sizeof(int), (char *)pin_info);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    for (i=0; i<myPins.numE; i++){
      j = lookup_GID(lookup_myHshEdges, myPins.eGIDs + (i*num_gid_entries));

      myPins.hashIdx[i] = j;    /* cache it for later */

      if (j < 0) FATAL_ERROR("Invalid global edge ID received");

      myHshEdges.numPins[j] += pin_info[i];
    }

    for (i=0; i<myPins.numE; i++){
      j = myPins.hashIdx[i]; 
      pin_info[i] = myHshEdges.numPins[j];
    }
    if (myPins.nHedges > 0){
      myPins.edgeGSize = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.nHedges);
      if (!myPins.edgeGSize) MEMORY_ERROR;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do_Reverse(myPins.ePlan, msg_tag, (char *)pin_info, 
                  sizeof(int), NULL, (char *)myPins.edgeGSize);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    ZOLTAN_FREE(&pin_info);

    /*
     * Process to which edge is assigned calculates edge weights. 
     */
    if (ew_dim && zz->Get_HG_Size_Edge_Wts && zz->Get_HG_Edge_Wts){
  
      /*
       * Get edge weights
       */
  
      zz->Get_HG_Size_Edge_Wts(
                   zz->Get_HG_Size_Edge_Wts_Data, &myEWs.size, &ierr);
  
  
      if ((ierr!=ZOLTAN_OK) && (ierr!=ZOLTAN_WARN)){
        FATAL_ERROR("obtaining edge weight size");
      }
  
      if (myEWs.size > 0){
  
        myEWs.edgeGID = ZOLTAN_MALLOC_GID_ARRAY(zz, myEWs.size);
        ew_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, myEWs.size);
        myEWs.wgt =
          (float *)ZOLTAN_MALLOC(sizeof(float) * myEWs.size * ew_dim);

        if (!myEWs.edgeGID || !ew_lids || !myEWs.wgt){
          MEMORY_ERROR;
        }
  
        zz->Get_HG_Edge_Wts(zz->Get_HG_Edge_Wts_Data,
                    zz->Num_GID, zz->Num_LID, myEWs.size, ew_dim,
                    myEWs.edgeGID, ew_lids, myEWs.wgt, &ierr);
  
        if ((ierr!=ZOLTAN_OK) && (ierr!=ZOLTAN_WARN)){
          FATAL_ERROR("obtaining edge weights");
        }
  
        /* 
         * Assign a process to each hyperedge using a hash function.
         */
  
        ierr = map_GIDs_to_processes(zz, myEWs.edgeGID, myEWs.size,
                 num_gid_entries, &myEWs.edgeHash, nProc);
  
        if ((ierr!=ZOLTAN_OK) && (ierr!=ZOLTAN_WARN)){
          goto End;
        }
      }
  
      /*
       * Send all edge weights to the process that was assigned that
       * edge by the hash function.  That process will combine the
       * weights according to the PHG_EDGE_WEIGHT_OPERATION parameter.
       */
  
      gid_requests = NULL;
      gid_weights = NULL;
  
      msg_tag--;
      ierr = Zoltan_Comm_Create(&plan, myEWs.size, myEWs.edgeHash,
               zz->Communicator, msg_tag, &nRequests);
  
      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }
    
      if (nRequests > 0){
        gid_requests = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
        gid_weights = (float *)ZOLTAN_MALLOC(sizeof(float)*ew_dim*nRequests);
  
        if ( (ew_dim && !gid_weights) || !gid_requests) MEMORY_ERROR;
      }
    
      msg_tag--;
      ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myEWs.edgeGID, 
               sizeof(ZOLTAN_ID_TYPE) * num_gid_entries,
               (char *)gid_requests);
  
      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }
  
      msg_tag--;
      ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myEWs.wgt,
               sizeof(float) * ew_dim,
               (char *)gid_weights);
  
      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }
  
      if (nRequests > 0){
        /*
         * Combine edge weights supplied for the same edge using
         * the PHG_EDGE_WEIGHT_OPERATION parameter.
         */
        have_wgt = (char *)ZOLTAN_CALLOC(myHshEdges.size, 1);
        if (myHshEdges.size && !have_wgt) MEMORY_ERROR;
        src = gid_weights;
  
        for (i=0; i < nRequests; i++, src += ew_dim){
          j = lookup_GID(lookup_myHshEdges, gid_requests + (i*num_gid_entries));
  
          if (j < 0){
            /* An edge for which there are no pins, ignore it. */
            continue;  
          }
  
          dest = myHshEdges.wgt + (j * ew_dim);
  
          if (!have_wgt[j]){
            for (w=0; w<ew_dim; w++){
              dest[w] = src[w];
            }
            have_wgt[j] = 1;
          }
          else{
            if (ew_op == PHG_FLAG_ERROR_EDGE_WEIGHTS){
              for (w=0; w<ew_dim; w++){
                if (src[w] != dest[w]){
                  FATAL_ERROR(
     "Different processes supplied different edge weights for the same edge");
                }
              }
            } else if (ew_op == PHG_MAX_EDGE_WEIGHTS){
              for (w=0; w<ew_dim; w++){
                if (src[w] > dest[w]){
                  dest[w] = src[w];
                }
              }
            } else if (ew_op == PHG_ADD_EDGE_WEIGHTS){
              for (w=0; w<ew_dim; w++){
                dest[w] += src[w];
              }
            }
          }
        }
        ZOLTAN_FREE(&gid_requests);
        ZOLTAN_FREE(&gid_weights);
        ZOLTAN_FREE(&have_wgt);
      }
      Zoltan_Comm_Destroy(&plan);
  
      /*
       * Obtain hyperedge weights from the process that computed them.
       */
  
      gid_weights = NULL;
  
      if (myPins.numE > 0){
        gid_weights = (float *)ZOLTAN_MALLOC(sizeof(float)*myPins.numE*ew_dim);
        if (ew_dim && !gid_weights) MEMORY_ERROR;
      }
    
      for (i=0; i<myPins.numE; i++){
        
        j = myPins.hashIdx[i];
  
        for (w=0; w < ew_dim; w++){
          gid_weights[i * ew_dim + w] = myHshEdges.wgt[j * ew_dim + w];
        }
      }
  
      if (myPins.nHedges > 0){
        myPins.ewgt = (float *)ZOLTAN_MALLOC(sizeof(float) * ew_dim
                                                           * myPins.nHedges);
        if (ew_dim && !myPins.ewgt) MEMORY_ERROR;
      }
  
      msg_tag--;
      ierr = Zoltan_Comm_Do_Reverse(myPins.ePlan, msg_tag, (char *)gid_weights, 
              sizeof(float) * ew_dim, NULL, (char *)myPins.ewgt);
  
      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }
  
      ZOLTAN_FREE(&gid_weights);
    }

    /* 
     * If we remove dense edges, do it now.  
     */
    nEdge = myPins.nHedges;

    ierr = Zoltan_HG_ignore_some_edges(zz, zhg, myObjs.GnVtx, 
               edgeSizeThreshold, final_output, 
               &myPins.nHedges,
               myPins.edgeGID, elids, myPins.esizes, myPins.edgeGSize,
               myPins.ewgt, myPins.pinGID, myPins.pinProc);

    if (ierr != ZOLTAN_OK){
      FATAL_ERROR("");
    }

    if (myPins.nHedges < nEdge){  /* some of my edges were removed */
      myPins.numPins = 0;
      for (i=0; i < myPins.nHedges; i++){
        myPins.numPins += myPins.esizes[i];
      }

      ZOLTAN_FREE(&myPins.edgeHash);

      ierr = map_GIDs_to_processes(zz, myPins.edgeGID, myPins.nHedges,
               num_gid_entries, &myPins.edgeHash, nProc);
    }

    /*
     * If edges were removed in any process, we need to recreate 
     * the edge comm plan, and rewrite the list of edges assigned to me.
     */
    nEdge -= myPins.nHedges;
    removedEdges = 0;
  
    rc = MPI_Allreduce(&nEdge, &cnt, 1, MPI_INT, MPI_SUM, zz->Communicator);
  
    if (cnt > 0){
      removedEdges = 1;
      Zoltan_Comm_Destroy(&myPins.ePlan);
  
      msg_tag--;
      ierr = Zoltan_Comm_Create(&myPins.ePlan, myPins.nHedges, myPins.edgeHash, 
               zz->Communicator, msg_tag, &myPins.numE);

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }
    
      msg_tag--;
      ierr = Zoltan_Comm_Do(myPins.ePlan, msg_tag, (char *)myPins.edgeGID, 
               sizeof(ZOLTAN_ID_TYPE) * num_gid_entries, 
               (char *)myPins.eGIDs);

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }

      keep_edge = (int *)ZOLTAN_CALLOC(sizeof(int) , myHshEdges.size);
      if (myHshEdges.size && !keep_edge) MEMORY_ERROR;

      for (i=0; i<myPins.numE; i++){
        j = lookup_GID(lookup_myHshEdges, myPins.eGIDs + (i*num_gid_entries));

        if (j < 0) FATAL_ERROR("Invalid global edge ID received");
        
        keep_edge[j] = 1;
      }

      changed = 0;
      for (i=0; i<myHshEdges.size; i++){
        if (keep_edge[i] == 0){
          changed = 1;
          break;
        }
      }

      if (changed){
        /* list of edges that hash to my process has changed because
         * some edges were removed
         */
        next_id = 0;
        for (i=0; i<myHshEdges.size; i++){
          if (keep_edge[i]){
            if (next_id < i){

              ZOLTAN_SET_GID(zz, 
                 myHshEdges.edgeGID + (next_id * num_gid_entries),
                 myHshEdges.edgeGID + (i * num_gid_entries));

              /* Don't need to re-write numPins, wgts because
               * we are done with them.
               */

            }
            next_id++;
          }
        }
        myHshEdges.size = next_id;

        free_GID_lookup_table(&lookup_myHshEdges);

        lookup_myHshEdges = 
          create_GID_lookup_table(myHshEdges.edgeGID, next_id,
                                  num_gid_entries);

        for (i=0; i<myPins.numE; i++){
          j = lookup_GID(lookup_myHshEdges, 
                         myPins.eGIDs + (i*num_gid_entries));

          if (j < 0) FATAL_ERROR("lookup_GID failure");

          myPins.hashIdx[i] = j;    /* cache it for later */
        }
      }
      ZOLTAN_FREE(&keep_edge);
    }

  } else if (graph_callbacks){

    /* Graph query functions, one hyperedge per vertex */

    ierr = Zoltan_HG_Graph_Callbacks(zz, zhg, use_all_neighbors,
                                     myObjs.GnVtx, edgeSizeThreshold,
                                     final_output, &myPins.nHedges,
                                     &myPins.edgeGID, &elids,
                                     &myPins.esizes, &myPins.ewgt,
                                     &myPins.numPins, &myPins.pinGID, 
                                     &myPins.pinProc);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
             "Error returned from Zoltan_HG_Graph_Callbacks.");
      goto End;
    }
    
    myPins.edgeGSize = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.nHedges);
    myPins.edgeGNO = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.nHedges);

    if (myPins.nHedges && (!myPins.edgeGSize || !myPins.edgeGNO)) MEMORY_ERROR;

    for (i=0; i<myPins.nHedges; i++){
      myPins.edgeGSize[i] = myPins.esizes[i];
    }

  } else {
    /* Partition without edge information?  Or return an error? */
    ZOLTAN_PRINT_WARN(zz->Proc, yo,
       "No edge information provided, partitioning vertices.");
  }

  /* Get pin vertex global number from vertex owner, and number of
   * hyperedges containing each vertex. 
   */
  myPins.pinGNO = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.numPins);

  if (myPins.numPins && !myPins.pinGNO) MEMORY_ERROR;

  if (myObjs.size > 0){
    memset(myObjs.numHedges, 0, sizeof(int) * myObjs.size);
    memset(myObjs.numAllHedges, 0, sizeof(int) * myObjs.size);
  }

  pin_requests = NULL;
  pin_info = NULL;
  msg_tag--;
  ierr = Zoltan_Comm_Create(&plan, myPins.numPins, myPins.pinProc, 
         zz->Communicator, msg_tag, &nRequests);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  if (nRequests > 0){
    pin_requests = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
    pin_info = (int *)ZOLTAN_MALLOC(sizeof(int) * nRequests);

    if (!pin_requests || !pin_info) MEMORY_ERROR;
  }

  msg_tag--;
  ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myPins.pinGID, 
           sizeof(ZOLTAN_ID_TYPE) * num_gid_entries, (char *)pin_requests);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  for (i=0; i<nRequests; i++){
    j = lookup_GID(lookup_myObjs, pin_requests + ( i * num_gid_entries));

    if (j < 0) FATAL_ERROR("Unexpected vertex GID received");

    pin_info[i] = myObjs.vtx_gno[j];
    myObjs.numHedges[j]++;
    myObjs.numAllHedges[j]++;
  }

  msg_tag--;
  ierr = Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *)pin_info, sizeof(int), 
                NULL, (char *)myPins.pinGNO);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }
  ZOLTAN_FREE(&pin_requests);
  ZOLTAN_FREE(&pin_info);

  Zoltan_Comm_Destroy(&plan);

  /***********************************************************************/
  /* We also need the total edges a vertex is in, including removed      */
  /* edges, when (add_obj_weight == PHG_ADD_PINS_WEIGHT).                */
  /***********************************************************************/
  pin_requests = NULL;
  if (removedEdges && (add_obj_weight == PHG_ADD_PINS_WEIGHT)) {
    msg_tag--;
    ierr = Zoltan_Comm_Create(&plan, zhg->nRemovePins,
           zhg->Remove_Pin_Procs,
           zz->Communicator, msg_tag, &nRequests);
  
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  
    if (nRequests > 0){
      pin_requests = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
      if (!pin_requests) MEMORY_ERROR;
    }
  
    msg_tag--;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)zhg->Remove_Pin_GIDs,
             sizeof(ZOLTAN_ID_TYPE) * num_gid_entries, (char *)pin_requests);
  
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  
    for (i=0; i<nRequests; i++){
      j = lookup_GID(lookup_myObjs, pin_requests + ( i * num_gid_entries));
  
      if (j < 0) FATAL_ERROR("Unexpected vertex GID received");
  
      myObjs.numAllHedges[j]++;
    }
  
    Zoltan_Comm_Destroy(&plan);
  }

  ZOLTAN_FREE(&pin_requests);
  ZOLTAN_FREE(&pin_info);

  /***********************************************************************/
  /* Impose a global hyperedge numbering */
  /***********************************************************************/

  nEdge = 0;
  egno = NULL;

  if (hypergraph_callbacks){
    egno = myHshEdges.edgeGNO;
    nEdge = myHshEdges.size;
  }
  else if (graph_callbacks){
    egno = myPins.edgeGNO;
    nEdge = myPins.nHedges;
  }

  if (randomizeInitDist) {  
    /* Randomize the input edges */
    int tmp;

    memset(mycnt, 0, zz->Num_Proc * sizeof(int));

    /* Compute random processor bin. */
    /* Temporarily store processor bin number in egno. */
    /* Count how many local vtxs selected processor bin */
    for (i = 0; i < nEdge; i++) {
      egno[i] = Zoltan_Rand_InRange(NULL, zz->Num_Proc);
      mycnt[egno[i]]++;
    }
    /* Compute prefix of mycnt */
    rc = MPI_Scan(mycnt, gcnt, zz->Num_Proc, MPI_INT, MPI_SUM,
                  zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)
    rc = MPI_Allreduce(mycnt, gtotal, zz->Num_Proc, MPI_INT, MPI_SUM, 
                  zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Compute first gno for vertices going to each target bin */
    for (tmp = 0, i = 0; i < zz->Num_Proc; i++) {
      gcnt[i] -= mycnt[i];
      tmp += gtotal[i];
      gtotal[i] = tmp - gtotal[i];
    }
    myHshEdges.GnEdge = gtotal[zz->Num_Proc] = tmp;


    /* Assign gnos sequential from gcnt[bin]. */
    for (i=0; i< nEdge; i++) {
      tmp = egno[i];
      egno[i] = gtotal[tmp] + gcnt[tmp];
      gcnt[tmp]++;
    }
  }
  else {
    rc = MPI_Scan (&nEdge, gtotal, 1, MPI_INT, MPI_SUM, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Gather data from all procs */

    rc = MPI_Allgather (&(gtotal[0]), 1, MPI_INT,
                   &(gtotal[1]), 1, MPI_INT, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)
    gtotal[0] = 0;
    myHshEdges.GnEdge = gtotal[nProc];

    /* Assign global numbers to edges. */
    for (i = 0; i < nEdge; i++)
      egno[i] = gtotal[zz->Proc] + i;
  }
  ZOLTAN_FREE(&gtotal);

  if (hypergraph_callbacks){

    /* Obtain edge global number for each edge for which I have pins. */

    if (myPins.numE > 0){
      gid_info = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.numE);
      if (!gid_info) MEMORY_ERROR;
    }
  
    for (i=0; i<myPins.numE; i++){
      j = myPins.hashIdx[i];
      gid_info[i] = myHshEdges.edgeGNO[j];
    }

    if (myPins.nHedges > 0){
      myPins.edgeGNO = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.nHedges);
      if (!myPins.edgeGNO) MEMORY_ERROR;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do_Reverse(myPins.ePlan, msg_tag, (char *)gid_info, 
                  sizeof(int), NULL, (char *)myPins.edgeGNO);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    ZOLTAN_FREE(&gid_info);

    /* Done with edge comm plan */

    Zoltan_Comm_Destroy(&myPins.ePlan);
    ZOLTAN_FREE(&myPins.hashIdx);
    ZOLTAN_FREE(&myPins.eGIDs);

    /* Done with consolidated edge info (except for myHshEdges.GnEdge) */

    free_GID_lookup_table(&lookup_myHshEdges);
    free_zoltan_temp_edges(&myHshEdges);
  }

  /***********************************************************************/
  /* Vertex weights                                                      */
  /***********************************************************************/

  if (add_obj_weight != PHG_ADD_NO_WEIGHT){
    add_vweight = 1;    /* may change in future to allow more than one */
    calcVwgt = (float *)ZOLTAN_CALLOC(sizeof(float), myObjs.size);
    if (myObjs.size && !calcVwgt) MEMORY_ERROR;

    if (add_obj_weight == PHG_ADD_UNIT_WEIGHT){
      for (i=0; i<myObjs.size; i++){
        calcVwgt[i] = 1.0;
      }
    }
    else if (add_obj_weight == PHG_ADD_PINS_WEIGHT){
      for (i=0; i<myObjs.size; i++){
        calcVwgt[i] = myObjs.numAllHedges[i];
      }
    }
  }
  else{
    add_vweight = 0;
    calcVwgt = NULL;
  }

  phg->VtxWeightDim =     /* 1 or greater */
    zz->Obj_Weight_Dim +   /* 0 or more application supplied weights*/
    add_vweight;           /* 0 or 1 additional calculated weights */

  if (phg->VtxWeightDim > 1) {
    /*
     * For now, only one weight per vertex will be used.  We will
     * still save multiple weights, because in the future we may
     * be able to use them.
     */
    if (zz->Proc == 0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Too many vertex weights.");
      ZOLTAN_PRINT_WARN(zz->Proc, yo, 
         "Multiple weights per vertex were supplied.");
      ZOLTAN_PRINT_WARN(zz->Proc, yo, 
        "Only the first application supplied weight per vertex will be used.");
    }
  }

  /***********************************************************************/
  /* 
   * Create list of vertex weights.  Some may have been provided by the
   * application in Zoltan_Get_Obj_List, another may have been calculated
   * above based on the value of the ADD_OBJ_WEIGHT parameter.
   */
  /***********************************************************************/

  if (calcVwgt){
    if (zz->Obj_Weight_Dim > 0){
      wgts = 
        (float *)ZOLTAN_MALLOC(sizeof(float) * myObjs.size * phg->VtxWeightDim);
      if ((myObjs.size * phg->VtxWeightDim) && !wgts){ 
        MEMORY_ERROR;
      }
  
      towgt = wgts;
      fromwgt = myObjs.vwgt;
  
      for (i=0; i<myObjs.size; i++){
        for (j=0; j < zz->Obj_Weight_Dim; j++){
          *towgt++ = *fromwgt++;
        }
        *towgt++ = calcVwgt[i];
      }
      ZOLTAN_FREE(&calcVwgt);
    }
    else{
      wgts = calcVwgt;
      calcVwgt = NULL;
    }

    ZOLTAN_FREE(&myObjs.vwgt);

    myObjs.vwgt = wgts;
  }

  if (zoltan_lb_eval){
    /* 
     * We were called from Zoltan_LB_Eval and we're done.
     * All edges, pins, pin owners, etc are in the "removed" lists.
     * We write vertex weights for vertices owned by this process 
     * to the hypergraph structure (where pin weights normally are).
     */

    phg->vwgt = myObjs.vwgt;
    myObjs.vwgt = NULL;

    goto End;
  }

  /* 
   * Compute the distribution of vertices and edges to the 2D data
   * distribution's processor columns and rows.
   * For now, these distributions are described by arrays dist_x
   * and dist_y; in the future, we may prefer a hashing function
   * mapping GIDs to processor columns and rows. KDDKDD
   */

  nProc_x = phg->comm->nProc_x;
  nProc_y = phg->comm->nProc_y;
  myProc_x = phg->comm->myProc_x;
  myProc_y = phg->comm->myProc_y;

  phg->dist_x = dist_x = (int *) ZOLTAN_CALLOC((nProc_x+1), sizeof(int));
  phg->dist_y = dist_y = (int *) ZOLTAN_CALLOC((nProc_y+1), sizeof(int));

  if (!dist_x || !dist_y) MEMORY_ERROR;

  frac_x = (float) myObjs.GnVtx / (float) nProc_x;
  for (i = 1; i < nProc_x; i++)
    dist_x[i] = (int) (i * frac_x);
  dist_x[nProc_x] = myObjs.GnVtx;
  
  frac_y = (float) myHshEdges.GnEdge / (float) nProc_y;
  for (i = 1; i < nProc_y; i++)
    dist_y[i] = (int) (i * frac_y);
  dist_y[nProc_y] = myHshEdges.GnEdge;
  
  /* myProc_y and myProc_x can be -1 when nProc is prime and we use a 2D
   * decomposition.  One processor is excluded from the 2D communicator;
   * for it, myProc_y and myProc_x == -1. */
  nEdge = (myProc_y >= 0 ? dist_y[myProc_y+1] - dist_y[myProc_y] : 0);
  nVtx  = (myProc_x >= 0 ? dist_x[myProc_x+1] - dist_x[myProc_x] : 0);

/*printf("%d) %d edges %d vertices\n",zz->Proc, nEdge, nVtx);*/

  if (method_repart){
    /* For REPARTITION, we add one vertex per partition and one edge 
     * per object in a repartition part (connecting the object with 
     * its input partition vertex).
     * Compute the number of these per processor within the 2D distribution
     * now so we can allocate relevant arrays to be large enough to include
     * this extra data.
     */

    /* Find number of repartition vertices and repartition edges to add. */
    nRepartVtx = 0;
    nRepartEdge = 0;
    for (i = 0; i < zhg->nObj; i++)
      if (zhg->Input_Parts[i] < zz->LB.Num_Global_Parts) {
        /* Add repartition vertices & edges only for parts < requested parts. */
        nRepartEdge++;  /* Will need a repartition edge for this object */
        if (zhg->Input_Parts[i] > nRepartVtx) nRepartVtx = zhg->Input_Parts[i];
      }
    nRepartVtx++;  /* # of repartition vtx == max partition + 1 */
    MPI_Allreduce(&nRepartVtx, &(zhg->GnRepartVtx), 1, MPI_INT, MPI_MAX,
                  zz->Communicator);
    MPI_Allreduce(&nRepartEdge, &(zhg->GnRepartEdge), 1, MPI_INT, MPI_SUM,
                  zz->Communicator);

    nRepartVtx = NumRepart(myProc_x, nProc_x, zhg->GnRepartVtx);
    /* For memory allocation and easy mapping of repartition edges to 
     * processor rows, compute maximum number of repart edges possible; 
     * will reduce later. */
    nRepartEdge = NumRepart(myProc_y, nProc_y, zhg->GnObj);
  }

  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution. 
   */

  proclist = (int *)ZOLTAN_MALLOC(MAX(myPins.numPins,myObjs.size)*sizeof(int));
  sendbuf = (int *) ZOLTAN_MALLOC(myPins.numPins * 2 * sizeof(int));

  cnt = 0; 
  for (i = 0; i < myPins.nHedges; i++) {
    /* processor row for the edge */
    edge_gno = myPins.edgeGNO[i];
    edge_Proc_y = EDGE_TO_PROC_Y(phg, edge_gno);

    for (j = 0; j < myPins.esizes[i]; j++) {
      /* processor column for the vertex */
      vtx_gno = myPins.pinGNO[cnt];
      vtx_Proc_x = VTX_TO_PROC_X(phg, vtx_gno);

      proclist[cnt] = edge_Proc_y * nProc_x + vtx_Proc_x;
      sendbuf[2*cnt] = edge_gno;
      sendbuf[2*cnt+1] = vtx_gno;
      cnt++;
    } 
  }

  /*
   * Send pins to their target processors.
   * They become non-zeros in the 2D data distribution.
   */

  msg_tag--;
  ierr = Zoltan_Comm_Create(&plan, cnt, proclist, zz->Communicator,
                     msg_tag, &nnz);

  if (nnz) {
    nonzeros = (int *) ZOLTAN_MALLOC(nnz * 2 * sizeof(int));
    if (!nonzeros) MEMORY_ERROR;
  }

  msg_tag--;
  Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(int),
                 (char *) nonzeros);

  Zoltan_Comm_Destroy(&plan);

  /* Unpack the non-zeros received. */

  tmparray = (int *) ZOLTAN_CALLOC(nEdge + 1, sizeof(int));
  hindex = (int *) ZOLTAN_CALLOC(nEdge + 1 + nRepartEdge, sizeof(int));
  hvertex = (int *) ZOLTAN_MALLOC((nnz + 2 * nRepartEdge) * sizeof(int)); 
                                  /* worst case size */

  if (!tmparray || !hindex || ((nnz || nRepartEdge) && !hvertex)){
    ZOLTAN_FREE(&hindex);
    ZOLTAN_FREE(&hvertex);
    MEMORY_ERROR;
  }

  /* Count the number of nonzeros per hyperedge */
  for (i = 0; i < nnz; i++) {
    idx = EDGE_GNO_TO_LNO(phg, nonzeros[2*i]); 
    tmparray[idx]++;
  }

  /* Compute prefix sum to represent hindex correctly. */
  for (i = 0; i < nEdge; i++)  {
    hindex[i+1] = hindex[i] + tmparray[i];
    tmparray[i] = 0;
  }
       
  for (i = 0; i < nnz; i++) {
    idx = EDGE_GNO_TO_LNO(phg, nonzeros[2*i]);
    hvertex[hindex[idx]+tmparray[idx]] = VTX_GNO_TO_LNO(phg, nonzeros[2*i+1]);
    tmparray[idx]++;
  }

  phg->nVtx = nVtx;
  phg->nEdge = nEdge;
  phg->nPins = nnz;
  phg->hindex = hindex;
  phg->hvertex = hvertex;

  /* Send vertex partition assignments and weights, if any. */
  /* Can use same plan for both. */

  if (phg->comm->nProc_x > 1) {

    /* Need a communication plan mapping GIDs to their GNOs processors
     * within a row communicator.  The plan is used to send vertex weights
     * and partition assignments to the 2D distribution and/or 
     * to create return lists after partitioning
     *
     * Since for 2D decomposition and prime nProc we exclude a processor,
     * we cannot use row_comm for this operation.  We'll simulate it,
     * allowing the excluded processor to send its data to row 0.
     */

    proc_offset = (myProc_y >= 0 ? myProc_y : 0) * nProc_x;
    for (i = 0; i < myObjs.size; i++)
      proclist[i] = proc_offset + VTX_TO_PROC_X(phg, myObjs.vtx_gno[i]);
      
    msg_tag++;
    ierr = Zoltan_Comm_Create(&(zhg->VtxPlan), myObjs.size, proclist, 
                              zz->Communicator, msg_tag, &nrecv);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
    zhg->nRecv_GNOs = nrecv;

    zhg->Recv_GNOs = recv_gno = (int *) ZOLTAN_MALLOC(nrecv * sizeof(int));

    if (nrecv && !recv_gno) MEMORY_ERROR;

    /* Use plan to send global numbers to the appropriate proc_x. */
    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) myObjs.vtx_gno, 
                          sizeof(int), (char *) recv_gno);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  }
  else {
    /* Save map of what needed. */
    zhg->nRecv_GNOs = nrecv = myObjs.size;
    zhg->Recv_GNOs = recv_gno = (int *) ZOLTAN_MALLOC(nrecv * sizeof(int));

    if (nrecv && !recv_gno) MEMORY_ERROR;
   
    memcpy(zhg->Recv_GNOs, myObjs.vtx_gno, nrecv * sizeof(int));
  }

  /* Send vertex partition assignments and weights to 2D distribution. */

  tmpparts = (int *) ZOLTAN_CALLOC(phg->nVtx, sizeof(int));
  *input_parts = (int *) ZOLTAN_MALLOC((phg->nVtx + nRepartVtx) * sizeof(int));

  dim = phg->VtxWeightDim;
  nwgt = (phg->nVtx + nRepartVtx) * dim;

  tmpwgts = (float *)ZOLTAN_CALLOC(sizeof(float), nwgt);
  phg->vwgt = (float *)ZOLTAN_MALLOC(sizeof(float) * nwgt);
  
  if (GnFixed) {                              
    tmpfixed   = (int*) ZOLTAN_MALLOC(phg->nVtx * sizeof(int));
    phg->fixed_part = (int*) ZOLTAN_MALLOC((phg->nVtx + nRepartVtx) * sizeof(int));

    if ((phg->nVtx || nRepartVtx) && (!tmpfixed || !phg->fixed_part))
      MEMORY_ERROR;
      
    for (i = 0 ; i < phg->nVtx; i++)
      tmpfixed[i] = -2;  
  } 

  if ((phg->nVtx || nRepartVtx) && 
       (!tmpparts || !*input_parts || !tmpwgts || !phg->vwgt)) MEMORY_ERROR;
  
  if (phg->comm->nProc_x == 1)  {
    for (i = 0; i < myObjs.size; i++) {
      idx = myObjs.vtx_gno[i];
      tmpparts[idx] = zhg->Input_Parts[i];
      for (j=0; j<dim; j++)
        tmpwgts[idx*dim + j] = myObjs.vwgt[i*dim + j];
      if (GnFixed)                                             
        tmpfixed[idx] = myObjs.fixed[i];
    }
  }
  else {
    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) zhg->Input_Parts,
                          sizeof(int), (char *) *input_parts);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) myObjs.vwgt,
                          sizeof(float) * dim, (char *) phg->vwgt);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
    
    if (GnFixed) {                                  
       msg_tag++;
       ierr = Zoltan_Comm_Do (zhg->VtxPlan, msg_tag, (char*) myObjs.fixed,
         sizeof(int), (char*) phg->fixed_part);
       if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN))
         goto End;         
    }
       
    for (i = 0; i < nrecv; i++) {
      idx = VTX_GNO_TO_LNO(phg, recv_gno[i]);
      tmpparts[idx] = (*input_parts)[i];
      for (j=0; j<dim; j++){
        tmpwgts[idx*dim + j] = phg->vwgt[i*dim + j];  
      if (GnFixed)                                       
        tmpfixed[idx] = phg->fixed_part[i];         
      }
    }
  }

  /* Reduce partition assignments and weights for all vertices within column 
   * to all processors within column.
   */

  if (phg->comm->col_comm != MPI_COMM_NULL){
    rc = MPI_Allreduce(tmpparts, *input_parts, phg->nVtx, MPI_INT, MPI_MAX, 
                  phg->comm->col_comm);
    CHECK_FOR_MPI_ERROR(rc)

    rc = MPI_Allreduce(tmpwgts, phg->vwgt, nwgt, MPI_FLOAT, MPI_MAX, 
                  phg->comm->col_comm);
    CHECK_FOR_MPI_ERROR(rc)

    if (GnFixed) {                                          
       rc = MPI_Allreduce(tmpfixed, phg->fixed_part, phg->nVtx, MPI_INT, MPI_MAX,
            phg->comm->col_comm);
       CHECK_FOR_MPI_ERROR(rc);
    }   
  }
  ZOLTAN_FREE(&tmpfixed);
  ZOLTAN_FREE(&tmpparts);
  ZOLTAN_FREE(&tmpwgts);

  /*  Send edge weights, if any */

  dim = zz->Edge_Weight_Dim;
  if (method_repart && (!dim))
    dim = 1; /* Need edge weights for REPARTITION; force malloc of ewgt array */
  phg->EdgeWeightDim = dim;
  nwgt = (phg->nEdge + nRepartEdge) * dim;

  if (dim) {
    tmpwgts   = (float *) ZOLTAN_CALLOC(nwgt, sizeof(float));
    phg->ewgt = (float *) ZOLTAN_CALLOC(nwgt, sizeof(float));

    if (nwgt && (!phg->ewgt || !tmpwgts)) MEMORY_ERROR;

    if (zz->Edge_Weight_Dim) {  /* Edge weights provided by application */
      if (phg->comm->nProc_y == 1) {
        for (i = 0; i < myPins.nHedges; i++) {
          idx = myPins.edgeGNO[i];
          for (j = 0; j < dim; j++)
            tmpwgts[idx * dim + j] = myPins.ewgt[i*dim + j];
        }
      }
      else {
        /* 
         * Since for 2D decomposition and prime nProc we exclude a processor,
         * we cannot use col_comm for this operation.  We'll simulate it,
         * allowing the excluded processor to send its data to col 0.
         */
        proc_offset = (myProc_x >= 0 ? myProc_x : 0);
        for (i = 0; i < myPins.nHedges; i++)
          proclist[i] = proc_offset 
                      + EDGE_TO_PROC_Y(phg, myPins.edgeGNO[i]) * nProc_x;
        
        msg_tag++;
  
        ierr = Zoltan_Comm_Create(&plan, myPins.nHedges, proclist, 
                                  zz->Communicator, msg_tag, &nrecv); 
  
        if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
          goto End;
        }
  
        /* Multiple processes may have weights for the same edge */
  
        recv_gno = (int *) ZOLTAN_MALLOC(nrecv * sizeof(int));
  
        gid_weights = (float *) ZOLTAN_CALLOC(nrecv*dim, sizeof(float));
  
        have_wgt = (char *)ZOLTAN_CALLOC(phg->nEdge, sizeof(char));
  
        if ((nrecv && (!recv_gno || !gid_weights)) ||
            (phg->nEdge && !have_wgt)){
  
          ZOLTAN_FREE(&recv_gno);
          MEMORY_ERROR;
        }
  
        msg_tag++;
        ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) myPins.edgeGNO,
                              sizeof(int), (char *) recv_gno);
  
        if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
          goto End;
        }
  
        msg_tag++;
        ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) myPins.ewgt, 
                              dim*sizeof(float), (char *) gid_weights);
  
        if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
          goto End;
        }
  
        Zoltan_Comm_Destroy(&plan);
  
        for (i = 0; i < nrecv; i++) {
          idx = EDGE_GNO_TO_LNO(phg, recv_gno[i]);
          if (have_wgt[idx]) continue;
          for (j = 0; j < dim; j++) 
            tmpwgts[idx * dim + j] = gid_weights[i * dim + j];
          have_wgt[idx] = 1;
        }
        ZOLTAN_FREE(&recv_gno); 
        ZOLTAN_FREE(&have_wgt); 
        ZOLTAN_FREE(&gid_weights); 
      }

      /* Need to gather weights for all edges within row 
       * to all processors within row.
       */
  
      if (phg->comm->row_comm != MPI_COMM_NULL && nwgt > 0){
        /* error here if numprocs < numrows */

        rc = MPI_Allreduce(tmpwgts, phg->ewgt, nwgt, MPI_FLOAT, MPI_MAX, 
                      phg->comm->row_comm);
        CHECK_FOR_MPI_ERROR(rc)

      }
    }
    else { /* dim > 0 but zz->Edge_Weight_Dim == 0 */
      /* Edge weights are needed for REPARTITION but are not provided by app */
      /* Set the edge weights for input vertices to 1. */
      for (i = 0; i < phg->nEdge; i++) phg->ewgt[i] = 1.;
    }
  }
  else {
    /* KDDKDD  For now, do not assume uniform edge weights.
     * KDDKDD  Can add later if, e.g., we decide to coalesce identical edges.
     */
    phg->EdgeWeightDim = 0;
  }

  if (method_repart){
    ierr = Zoltan_PHG_Add_Repart_Data(zz, zhg, phg,
                                      myObjs.vtx_gno, hgp, *input_parts);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error from Zoltan_PHG_Add_Repart_Data");
      goto End;
    }
  } else {
      if (hgp->final_output) {
          if ((ierr = getObjectSizes(zz, zhg))!=ZOLTAN_OK)
              goto End;
      }
  }

  if (!zz->LB.Remap_Flag && zz->LB.Return_Lists == ZOLTAN_LB_NO_LISTS) {
    int gnremove;
    rc = MPI_Allreduce(&(zhg->nRemove), &gnremove, 1, MPI_INT, MPI_SUM, 
                  zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)
    if (!final_output || !gnremove) {
      /* Don't need the plan long-term; destroy it now. */
      Zoltan_Comm_Destroy(&(zhg->VtxPlan));
      ZOLTAN_FREE(&(zhg->Recv_GNOs));
      zhg->nRecv_GNOs = 0;
    }
  }

    
  ierr = Zoltan_HG_Create_Mirror(zz, phg);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error from Zoltan_HG_Create_Mirror");
    goto End;
  }

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    Zoltan_PHG_Free_Hypergraph_Data(zhg);
  }

  ZOLTAN_FREE(&tmpfixed);
  
  free_zoltan_objects(&myObjs);
  free_zoltan_pins(&myPins);
  free_zoltan_ews(&myEWs);
  free_zoltan_temp_edges(&myHshEdges);
  free_zoltan_temp_vertices(&myHshVtxs);

  free_GID_lookup_table(&lookup_myObjs);
  free_GID_lookup_table(&lookup_myHshEdges);
  free_GID_lookup_table(&lookup_myHshVtxs);

  Zoltan_Comm_Destroy(&plan);

  Zoltan_Multifree(__FILE__, __LINE__, 18, 
    &gtotal,
    &pinIdx,
    &elids,
    &pin_info,
    &keep_edge,
    &pin_requests,
    &gid_info,
    &calcVwgt,
    &ew_lids,
    &gid_requests,
    &gid_weights,
    &have_wgt,
    &proclist,
    &sendbuf,
    &tmparray,
    &tmpparts,
    &tmpwgts,
    &nonzeros);

#ifdef DEBUG_FILL_HYPERGRAPH
  for (i=0; i<zz->Num_Proc; i++){
    if (i == zz->Proc){
      printf("HYPERGRAPH on process %d\n",i);

      if (add_obj_weight == PHG_ADD_PINS_WEIGHT) 
         /* print sum of pin weights */
        print_hypergraph(zz, zhg, zz->Obj_Weight_Dim);
      else
        print_hypergraph(zz, zhg, -1);

      printf("\n");
      fflush(stdout);
    }
  
    rc = MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
/*****************************************************************************/
static void free_zoltan_objects(zoltan_objects *zo)
{
  if (zo == NULL) return;
  /* Don't free zo->vtxGID, it's pointing into ZHG structure */
  ZOLTAN_FREE(&(zo->vtx_gno));
  ZOLTAN_FREE(&(zo->vwgt));
  ZOLTAN_FREE(&(zo->vtxHash));
  ZOLTAN_FREE(&(zo->numHedges));
  ZOLTAN_FREE(&(zo->numAllHedges));
  ZOLTAN_FREE(&zo->fixed);
}
/*****************************************************************************/
static void free_zoltan_pins(zoltan_pins *zp)
{
  if (zp == NULL) return;
  ZOLTAN_FREE(&(zp->edgeGID));
  ZOLTAN_FREE(&(zp->esizes));
  ZOLTAN_FREE(&(zp->pinGID));
  ZOLTAN_FREE(&(zp->vtxHash));
  ZOLTAN_FREE(&(zp->pinGNO));
  ZOLTAN_FREE(&(zp->pinProc));
  ZOLTAN_FREE(&(zp->edgeHash));
  ZOLTAN_FREE(&(zp->ewgt));
  ZOLTAN_FREE(&(zp->edgeGNO));
  ZOLTAN_FREE(&(zp->edgeGSize));
  ZOLTAN_FREE(&(zp->eGIDs));
  ZOLTAN_FREE(&(zp->hashIdx));
  Zoltan_Comm_Destroy(&(zp->ePlan));
}
/*****************************************************************************/
static void free_zoltan_ews(zoltan_ews *zew)
{
  if (zew == NULL) return;
  ZOLTAN_FREE(&(zew->edgeGID));
  ZOLTAN_FREE(&(zew->edgeHash));
  ZOLTAN_FREE(&(zew->wgt));
}
/*****************************************************************************/
static void free_zoltan_temp_edges(zoltan_temp_edges *zte)
{
  if (zte == NULL) return;
  ZOLTAN_FREE(&(zte->edgeGID));
  ZOLTAN_FREE(&(zte->edgeGNO));
  ZOLTAN_FREE(&(zte->numPins));
  ZOLTAN_FREE(&(zte->wgt));
}
/*****************************************************************************/
static void free_zoltan_temp_vertices(zoltan_temp_vertices *ztv)
{
  if (ztv == NULL) return;
  ZOLTAN_FREE(&(ztv->vtxGID));
  ZOLTAN_FREE(&(ztv->vtxOwner));
}

/*****************************************************************************/
#if 0
#include <unistd.h>
static void print_zoltan_pins(zoltan_pins *z, int me, int ewgt_dim)
{
int i, j, k;

  sleep(me);
  printf("%d) %d hyperedges\n\n",me, z->nHedges);

  if (z->nHedges == 0) return;

  k = 0;
  for (i=0; i<z->nHedges; i++){
    if (z->edgeHash){
      printf("  GID %d, hashed to %d, num pins locally %d, GNO %d, num pins globally %d\n", 
               z->edgeGID[i], z->edgeHash[i], z->esizes[i], z->edgeGNO[i], z->edgeGSize[i]);
    }
    else{
      printf("  GID %d, num pins locally %d, GNO %d, num pins globally %d\n", 
               z->edgeGID[i], z->esizes[i], z->edgeGNO[i], z->edgeGSize[i]);
    }
    if (ewgt_dim > 0){
      for (j=0; j<ewgt_dim; j++){
        if (!j) printf("      edge weight ");
        printf(" %f ", z->ewgt[i*ewgt_dim + j]);
      }
      printf("\n");
    }
    for (j=0; j<z->esizes[i]; j++,k++){
      if (z->vtxHash){
        printf("      GID %d, hashed to proc %d, GNO %d, proc %d\n", 
          z->pinGID[k],z->vtxHash[k],z->pinGNO[k], z->pinProc[k]);
      }
      else{
        printf("      GID %d, GNO %d, proc %d\n", z->pinGID[k],z->pinGNO[k], z->pinProc[k]);
      }
    }
  }
  printf("\n");
  
}
#endif
/*****************************************************************************/
static int map_GIDs_to_processes(ZZ *zz, ZOLTAN_ID_PTR eid, int size, 
                             int lenGID, int **hashedProc, int nprocs)
{
int i, j;
int *procList;
static char *yo = "map_GIDs_to_processes";

  *hashedProc = NULL;

  if (size < 1){
    return ZOLTAN_OK;
  }

  procList = (int *)ZOLTAN_MALLOC(sizeof(int) * size);

  if (!procList){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); 
    return ZOLTAN_MEMERR;
  }

  for (i=0; i<size; i++){
    j = Zoltan_Hash(eid, lenGID, nprocs);
    procList[i] = j;
    eid += lenGID;
  }

  *hashedProc = procList;

  return ZOLTAN_OK;
}


/*****************************************************************************/

int Zoltan_PHG_Removed_Cuts(
  ZZ *zz,
  ZHG *zhg,
  double *localcuts /* Array of length 2. Upon return:
                     * localcuts[0] = ConCut: Sum_over_edges( (nparts-1)*ewgt )
                     * localcuts[1] = NetCut: Sum_over_edges( (nparts>1)*ewgt )
                     */
)
{
/* Function to compute the cuts of removed hyperedges.
 * This routine is needed ONLY for accurately computing cut statistics.
 * Depending on which query functions were defined,
 * pins for removed hyperedges may not have been retrieved before.
 * They must be retrieved now.
 */
static char *yo = "Zoltan_PHG_Removed_Cuts";
int ierr = ZOLTAN_OK;
int i, j, k, cnt, ncnt, max_parts;
struct Hash_Node *hash_nodes = NULL;  /* Hash table variables for mapping   */
struct Hash_Node **hash_tab = NULL;   /* GIDs to global numbering system.   */
int npins = 0;                   /* # of pins in removed hyperedges */
ZOLTAN_ID_PTR pins = NULL;       /* pins for removed edges */
int *pin_procs = NULL;           /* procs owning pins for removed edges */
int *pin_parts = NULL;           /* parts computed for pins in removed edges */
ZOLTAN_COMM_OBJ *plan = NULL;
int nrecv;                       /* # of requests for pin info */
ZOLTAN_ID_PTR recvpins = NULL;   /* Requested pins from other procs */
int *outparts = NULL;            /* received partition info for pins */
int *parts = NULL;               /* Array of size Num_Global_Parts indicating
                                    whether the current edge has a pin in
                                    parts[i] */
double loccuts[2];               /* Local cut values: [0]=ConCut; [1]=NetCut */
int nObj = zhg->nObj;
ZOLTAN_ID_PTR gids = zhg->GIDs;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int myinfo[2], allinfo[2];
int missingPins;
int msg_tag = 23132;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Get pin information for removed hyperedges if necessary. */

  for (i = 0; i < zhg->nRemove; i++) npins += zhg->Remove_Esize[i];

  if (zhg->nRemove && npins != 0) {
    if (zhg->Remove_Pin_Procs){ 
     /* 
      * Hypergraph callbacks built graph previously.  They are
      * expensive so we don't want to repeat them if possible.
      */
      pins = zhg->Remove_Pin_GIDs;
      pin_procs = zhg->Remove_Pin_Procs;
    }
    else if (zz->Get_Num_Edges) { /* Graph callbacks used to build hypergraph */
      ZOLTAN_ID_PTR lid;
      int ewgtdim = zz->Edge_Weight_Dim;
      float *gewgts = NULL;   /* Graph edge weights */
  
      pins = ZOLTAN_MALLOC_GID_ARRAY(zz, (npins + zhg->nRemove));
      pin_procs = (int *) ZOLTAN_MALLOC((npins + zhg->nRemove) * sizeof(int));
      if (ewgtdim)
        gewgts = (float *) ZOLTAN_MALLOC(npins * ewgtdim * sizeof(float));
      if (!pins || !pin_procs || (ewgtdim && !gewgts)) MEMORY_ERROR;
  
      if (zz->Get_Edge_List_Multi) 
        zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data,
                          num_gid_entries, num_lid_entries, zhg->nRemove,
                          zhg->Remove_EGIDs, zhg->Remove_ELIDs,
                          zhg->Remove_Esize, pins, pin_procs, ewgtdim, 
                          gewgts, &ierr);
      else {
        cnt = 0;
        for (i = 0; i < zhg->nRemove; i++) {
          lid = (num_lid_entries ? &(zhg->Remove_ELIDs[i*num_lid_entries])
                                 : NULL);
          zz->Get_Edge_List(zz->Get_Edge_List_Data,
                            num_gid_entries, num_lid_entries,
                            &(zhg->Remove_EGIDs[i*num_gid_entries]), lid,
                            &(pins[cnt]), &(pin_procs[cnt]),
                            ewgtdim, &(gewgts[cnt*ewgtdim]),
                            &ierr);
          cnt += zhg->Remove_Esize[i];
        }
      }
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                           "Error returned from getting Edge_Lists");
        goto End;
      }
  
      /* Post-process edges to add vgid[i] to hedge i. */
      cnt = npins;
      npins += zhg->nRemove;  /* Will add egid to hyperedge */
      ncnt = npins;
      for (i = zhg->nRemove-1; i >= 0; i--) {
        /* Copy the existing pins for the edge */
        for (j = 0; j < zhg->Remove_Esize[i]; j++) {
          cnt--;
          ncnt--;
          ZOLTAN_SET_GID(zz, &(pins[ncnt]), &(pins[cnt]));
          pin_procs[ncnt] = pin_procs[cnt];
          for (k = 0; k < ewgtdim; k++)   /* sum the graph-edge wgts? */
            zhg->Remove_Ewgt[i*ewgtdim + k] += gewgts[cnt*ewgtdim + k];
        }
        /* Add egid[i] */
        ncnt--;
        ZOLTAN_SET_GID(zz, &(pins[ncnt]),
                           &(zhg->Remove_EGIDs[i*num_gid_entries]));
        pin_procs[ncnt] = zz->Proc;
        zhg->Remove_Esize[i]++;
      }
      ZOLTAN_FREE(&gewgts);
    }
    else{
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
        "No way to obtain processes owning pins");
      goto End;
    }
  }

  /* Communicate Output_Part info for off-processor pins */

  Zoltan_Comm_Create(&plan, npins, pin_procs, zz->Communicator,
                     msg_tag, &nrecv);

  if (nrecv) {
    recvpins = ZOLTAN_MALLOC_GID_ARRAY(zz, nrecv);
    outparts = (int *) ZOLTAN_MALLOC(nrecv * sizeof(int));
    if (!recvpins || !outparts) MEMORY_ERROR;
  }

  Zoltan_Comm_Do(plan, msg_tag, (char *) pins, 
                 num_gid_entries * sizeof(ZOLTAN_ID_TYPE), (char *) recvpins);

  if (!zhg->Remove_Pin_GIDs){
    ZOLTAN_FREE(&pins);
  }

  if (nrecv) {
    hash_nodes = (struct Hash_Node *) ZOLTAN_MALLOC(nObj * 
                                                    sizeof(struct Hash_Node));
    hash_tab = (struct Hash_Node **) ZOLTAN_MALLOC(nObj *
                                                   sizeof(struct Hash_Node *));
    if (!hash_nodes || !hash_tab) MEMORY_ERROR;

    /* Create hash table of local indices to speed searches for received
     * pin requests.
     */

    for (i = 0; i < nObj; i++) {
      hash_tab[i] = NULL;
      hash_nodes[i].gid = &(gids[i*num_gid_entries]);
      hash_nodes[i].gno = i;    /* LOCAL index here; change name for gno?  */
    }

    for (i=0; i< nObj; i++) {
      /* insert hashed elements into hash table */
      j = Zoltan_Hash(&(gids[i*num_gid_entries]), num_gid_entries,
                      (unsigned int) nObj);
      hash_nodes[i].next = hash_tab[j];
      hash_tab[j] = &hash_nodes[i];
    }

    /* Gather partition information for pin requests */
    for (i = 0; i < nrecv; i++) {
      j = hash_lookup(zz, &(recvpins[i*num_gid_entries]), nObj, hash_tab, NULL);
      outparts[i] = zhg->Output_Parts[j];
    }
    
    ZOLTAN_FREE(&hash_nodes);
    ZOLTAN_FREE(&hash_tab);
    ZOLTAN_FREE(&recvpins);
  }

  /* Send partition info back to requesting processor */
  pin_parts = (int *) ZOLTAN_MALLOC(npins * sizeof(int));
  Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *) outparts, 
                         sizeof(int), NULL, (char *) pin_parts);

  ZOLTAN_FREE(&outparts);

  Zoltan_Comm_Destroy(&plan);

  /* Compute the cut metrics using received partition info.
   *
   * Get maximum possible partition number.  (It's not always the
   * value found in the ZZ structure.)
   *
   * Determine whether edges are split across processes.  If not,
   * then each process can count the number of cuts for it's edges.
   *
   * If edges are split across processes, then hash each edge GID
   * to a process, send that process the list of partitions the
   * edge spans, and let that process report the number of cuts.
   */

  myinfo[0] = myinfo[1] = 0;
  for (cnt=0, i = 0; i < zhg->nRemove; i++) {
    for (j = 0; j < zhg->Remove_Esize[i]; j++) {
      if (pin_parts[cnt] >= myinfo[0]) myinfo[0] = pin_parts[cnt]+1;
      cnt++;
    }
    myinfo[1] += (zhg->Remove_GEsize[i] - zhg->Remove_Esize[i]);
  }

  MPI_Allreduce(myinfo, allinfo, 2, MPI_INT, MPI_MAX, zz->Communicator);

  max_parts = allinfo[0];
  missingPins = allinfo[1];

  if (missingPins > 0){
    ierr = removed_cuts_global(zz, zhg, max_parts, pin_parts, loccuts, msg_tag);
    msg_tag -= 10;
  }
  else{
    ierr = removed_cuts_local(zz, zhg, max_parts, pin_parts, loccuts);
  }

  localcuts[0] = loccuts[0];
  localcuts[1] = loccuts[1];

End:

  if (!zhg->Remove_Pin_Procs){
    ZOLTAN_FREE(&pin_procs);
  }
  ZOLTAN_FREE(&pin_parts);
  ZOLTAN_FREE(&parts);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}
/*****************************************************************************/
static int removed_cuts_global(ZZ *zz, ZHG *zhg, 
                int max_parts, int *pin_parts, double *loccuts, int tag)
{
char *yo = "removed_cuts_global";
int ierr = ZOLTAN_OK;
int *nextpart;
int i, j, idx, nparts, nedges, numPartitions;
float ewgt;
GID_lookup *lookup=NULL;
ZOLTAN_COMM_OBJ *plan=NULL;
ZOLTAN_ID_PTR global_ids=NULL, eGID=NULL;
int *edgeHash = NULL, *eSizes=NULL, *eParts=NULL;
char **parts=NULL;
float *eWgts=NULL, *weights=NULL;
int numUniqueEdges = 0;
int wdim = zz->Edge_Weight_Dim;
int lenGID = zz->Num_GID;


  /* Assign each edge GID to a process using a hash function */

  ierr = map_GIDs_to_processes(zz, zhg->Remove_EGIDs, zhg->nRemove,
               lenGID, &edgeHash, zz->Num_Proc);

  if (ierr != ZOLTAN_OK){
    goto End;
  }

  /* Send edge GIDs, sizes, and weights to their assigned processes.  */
  /* Incoming list of edges is not unique.                            */

  tag--;
  ierr = Zoltan_Comm_Create(&plan, zhg->nRemove, edgeHash,
               zz->Communicator, tag, &nedges);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  if (nedges > 0){
    eGID = ZOLTAN_MALLOC_GID_ARRAY(zz, nedges);
    eSizes = (int *)ZOLTAN_MALLOC(sizeof(int) * nedges);
    eWgts = (float *)ZOLTAN_MALLOC(sizeof(float) * nedges * wdim);

    if (!eGID || !eSizes) MEMORY_ERROR;
    if (wdim && !eWgts) MEMORY_ERROR;
  }

  tag--;
  ierr = Zoltan_Comm_Do(plan, tag, (char *)zhg->Remove_EGIDs,
           sizeof(ZOLTAN_ID_TYPE) * lenGID, (char *)eGID);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  tag--;
  ierr = Zoltan_Comm_Do(plan, tag, (char *)zhg->Remove_Esize,
           sizeof(int), (char *)eSizes);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  if (wdim > 0){
    tag--;
    ierr = Zoltan_Comm_Do(plan, tag, (char *)zhg->Remove_Ewgt,
             sizeof(float) * wdim, (char *)eWgts);
  
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  }

  /* Send the pin partitions to their assigned processes. */

  tag--;
  ierr = Zoltan_Comm_Resize(plan, zhg->Remove_Esize, tag, &numPartitions);

  eParts = (int *)ZOLTAN_MALLOC(sizeof(int) * numPartitions);

  if (numPartitions && !eParts) MEMORY_ERROR;

  tag--;
  ierr = Zoltan_Comm_Do(plan, tag, (char *)pin_parts, sizeof(int), 
                        (char *)eParts);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  /* Create unique list of edge GIDs and mapping from GID to
     index in this list.
   */

  global_ids = ZOLTAN_MALLOC_GID_ARRAY(zz, nedges);
  if (nedges && !global_ids) MEMORY_ERROR;

  memcpy(global_ids, eGID, nedges * lenGID * sizeof(ZOLTAN_ID_TYPE));

  lookup = create_GID_lookup_table2(global_ids, nedges, lenGID);

  numUniqueEdges = lookup->numGIDs;

  /* Now count up the cuts in each edge */

  parts = (char **) ZOLTAN_MALLOC(sizeof(char *) * numUniqueEdges);
  if (numUniqueEdges && !parts) MEMORY_ERROR;

  for (i=0; i<numUniqueEdges; i++){
    parts[i] = (char *)ZOLTAN_CALLOC(sizeof(char), max_parts);
    if (!parts[i]) MEMORY_ERROR;
  }

  weights = (float *) ZOLTAN_MALLOC(sizeof(float) * numUniqueEdges * wdim);
  if (numUniqueEdges && wdim && !weights) MEMORY_ERROR;

  nextpart = eParts;

  for (i = 0; i < nedges; i++) {
    idx = lookup_GID(lookup, eGID + (i * lenGID));

    for (j=0; j<eSizes[i]; j++){
      parts[idx][*nextpart++] = 1;
    }

    for (j=0; j < wdim; j++){
      weights[idx*wdim + j] = eWgts[i*wdim + j];
    }
  }

  loccuts[0] = loccuts[1] = 0.;
  for (i = 0; i < numUniqueEdges; i++) {
    nparts = 0;
    for (j = 0; j < max_parts; j++) {
      if (parts[i][j]) nparts++;
    }
    ewgt = (wdim ? weights[i*wdim] : 1.);
    if (nparts > 1) {
      loccuts[0] += (nparts-1) * ewgt;
      loccuts[1] += ewgt;
    }
  }

End:
  free_GID_lookup_table(&lookup);
  ZOLTAN_FREE(&edgeHash);
  ZOLTAN_FREE(&eGID);
  ZOLTAN_FREE(&eSizes);
  ZOLTAN_FREE(&eWgts);
  ZOLTAN_FREE(&eParts);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&weights);
  Zoltan_Comm_Destroy(&plan);

  if (parts){
    for (i=0; i<numUniqueEdges; i++){
     ZOLTAN_FREE(&(parts[i]));
    }
    ZOLTAN_FREE(&parts);
  }

  return ierr;
}
/*****************************************************************************/
static int removed_cuts_local(ZZ *zz, ZHG *zhg, 
                int max_parts, int *pin_parts, double *loccuts)
{
char *yo = "removed_cuts_local";
int i, cnt, j, ierr, nparts;
int *parts;
float ewgt;

  ierr = ZOLTAN_OK;
  parts = (int *) ZOLTAN_CALLOC(max_parts, sizeof(int));
  if (max_parts && !parts) MEMORY_ERROR;

  cnt = 0;
  loccuts[0] = loccuts[1] = 0.;
  for (i = 0; i < zhg->nRemove; i++) {
    nparts = 0;
    for (j = 0; j < zhg->Remove_Esize[i]; j++) {
      if (parts[pin_parts[cnt]] < i+1) {
        nparts++;
      }
      parts[pin_parts[cnt]] = i+1;
      cnt++;
    }
    ewgt = (zhg->Remove_Ewgt ? zhg->Remove_Ewgt[i*zz->Edge_Weight_Dim] : 1.);
    if (nparts > 1) {
      loccuts[0] += (nparts-1) * ewgt;
      loccuts[1] += ewgt;
    }
  }

End:
  ZOLTAN_FREE(&parts);
  return ierr;
}
/*****************************************************************************/

static int hash_lookup(
  ZZ *zz,
  ZOLTAN_ID_PTR key,
  int nVtx,
  struct Hash_Node **hash_tab,
  struct Hash_Node **hn
)
{
/* Looks up a key GID in the hash table; returns its gno. */
/* Based on hash_lookup in build_graph.c. */

  int i;
  struct Hash_Node *ptr;

  if (hash_tab){
    if (hn) *hn = NULL;
  
    i = Zoltan_Hash(key, zz->Num_GID, (unsigned int) nVtx);
    for (ptr=hash_tab[i]; ptr != NULL; ptr = ptr->next){
      if (ZOLTAN_EQ_GID(zz, ptr->gid, key)){
        if (hn) *hn = ptr;
        return (ptr->gno);
      }
    }
  }
  /* Key not in hash table */
  return -1;
}

/****************************************************************************/
#ifdef DEBUG_FILL_HYPERGRAPH
static void print_hypergraph(ZZ *zz, ZHG *zhg, int sumWeight)
{
  int i, j, npins, nremovedpins;
  int ewdim = zz->Edge_Weight_Dim;
  int vwdim = zhg->HG.VtxWeightDim;
  float *wgt, *vwgt, sum;
  int *pin, *owner, *lno;
  HGraph *hg = &zhg->HG;

  printf("VERTICES (%d): (gid lid inpart outpart)\n",zhg->nObj);

  for (i=0; i<zhg->nObj; i++){
    if (zhg->Input_Parts && zhg->Output_Parts){
      printf("  %d:  %d  %d  %d  %d\n", i, zhg->GIDs[i], zhg->LIDs[i],
               zhg->Input_Parts[i], zhg->Output_Parts[i]);
    }
    else{
      printf("  %d:  %d  %d \n", i, zhg->GIDs[i], zhg->LIDs[i]);
    }
  }

  printf("REMOVED EDGES (%d): (gid lid [local-size/global-size] pins/owners %d weights)\n",
                  zhg->nRemove, ewdim);

  pin = zhg->Remove_Pin_GIDs;
  owner = zhg->Remove_Pin_Procs;
  wgt = zhg->Remove_Ewgt;
  nremovedpins = 0;

  for (i=0; i<zhg->nRemove; i++){
    printf("  %d:  %d  %d  [%d/%d]  ",i, 
        zhg->Remove_EGIDs[i], zhg->Remove_ELIDs[i],
        zhg->Remove_Esize[i], zhg->Remove_GEsize[i]);
    for (j=0; j<zhg->Remove_Esize[i]; j++){
      if (j && (j%10==0)){
        printf("\n      ");
      }
      printf("%d/%d ",*pin++, *owner++);
    }
    for (j=0; j<ewdim; j++){
      printf(" %f ", *wgt++);
    }
    printf("\n");

    nremovedpins += zhg->Remove_Esize[i];
  }
  if (zhg->nRemove > 0){
    printf("  Total pins removed: %d\n",nremovedpins);
  }

  printf("%d EDGES (%d weights), %d total PINS:\n",
          hg->nEdge, ewdim, hg->nPins);

  wgt = hg->ewgt;
  lno = hg->hvertex;
  vwgt = hg->vwgt;

  for (i=0; i<hg->nEdge; i++){
    npins = hg->hindex[i+1] - hg->hindex[i];

    printf(" edge %d: ",EDGE_LNO_TO_GNO(hg, i));
    for (j=0; j<ewdim; j++){
      printf(" %f",*wgt++);
    }
    printf("\n %d pins: ", npins);
    for (j=0; j<npins; j++){
      printf("%d ", *lno++);
    }
    printf("\n");
  }

  printf("%d PIN global numbers and %d weights:\n", hg->nVtx, vwdim);

  sum = 0;

  for (i=0; i<hg->nVtx; i++){
    printf("  %d  %d: ", i, VTX_LNO_TO_GNO(hg, i));
    for (j=0; j<vwdim; j++){
      if (j==sumWeight) sum += *vwgt;
      printf("%f ", *vwgt++);
    }
    printf("\n");
  }
  if (sum > 0.0) printf("Weight %d sums to %f\n",sumWeight+1,sum);
}
#endif

/****************************************************************************/
/* 
 * Create, access and delete a hash table mapping a GID its
 * position in a list.
 */

static void free_GID_lookup_table(GID_lookup **lu)
{
  GID_lookup *l = *lu;

  if (l == NULL) return;

  ZOLTAN_FREE(&l->htTop);
  ZOLTAN_FREE(&l->ht);
  ZOLTAN_FREE(lu);
}
/****************************************************************************/
static GID_lookup *create_GID_lookup_table(ZOLTAN_ID_PTR gids, int size, int lenGID)
{
  int i, j, tsize;
  GID_lookup *lu = NULL;

  lu = (GID_lookup *)ZOLTAN_MALLOC(sizeof(GID_lookup));
  if (!lu){
    return NULL;
  }

  tsize = size * 1.25;

  lu->htTop = (struct Hash_Node *)ZOLTAN_MALLOC(sizeof(struct Hash_Node)*size);
  lu->ht = (struct Hash_Node **)ZOLTAN_CALLOC(sizeof(struct Hash_Node*), tsize);

  if (tsize && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    ZOLTAN_FREE(&lu);
    return NULL;
  }

  lu->table_size = tsize;
  lu->numGIDs = size;
  lu->lenGID = lenGID;

  for (i=0; i<size; i++){
    lu->htTop[i].gid = gids + (i * lenGID);
    lu->htTop[i].gno = i;

    j = Zoltan_Hash(lu->htTop[i].gid, lenGID, tsize);
    
    lu->htTop[i].next = lu->ht[j];
    lu->ht[j] = lu->htTop + i;
  }
 
  return lu;
}

/****************************************************************************/
/* gid list is not unique.  rewrite it as a list of unique gids */

static GID_lookup *create_GID_lookup_table2(ZOLTAN_ID_PTR gids, int ngids, int lenGID)
{
  int i, j, k, tsize, found;
  struct Hash_Node *hn;
  ZOLTAN_ID_PTR nextGID, nextUniqueGID;
  GID_lookup *lu = NULL;

  tsize = ngids;    /* actually may be larger than number of unique ids */
  
  nextGID = nextUniqueGID = gids;

  lu = (GID_lookup *)ZOLTAN_MALLOC(sizeof(GID_lookup));
  if (!lu){
    return NULL;
  }

  lu->ht = (struct Hash_Node **)ZOLTAN_CALLOC(sizeof(struct Hash_Node*) , tsize);
  hn = lu->htTop = (struct Hash_Node *)ZOLTAN_MALLOC(sizeof(struct Hash_Node) * ngids);

  if (tsize && (!lu->htTop || !lu->ht)){
    ZOLTAN_FREE(&lu);
    ZOLTAN_FREE(&lu->htTop);
    ZOLTAN_FREE(&lu->ht);
    return NULL;
  }

  lu->lenGID = lenGID;
  lu->table_size = tsize;
  lu->numGIDs = 0;

  for (i=0; i<ngids; i++, nextGID += lenGID){

    found = lookup_GID(lu, nextGID);

    if (found >= 0) continue;

    hn->gid = nextUniqueGID;
    hn->gno = lu->numGIDs;

    if (nextUniqueGID < nextGID){
      for (k=0; k<lenGID; k++){
        nextUniqueGID[k] = nextGID[k]; 
      }
    }

    j = Zoltan_Hash(nextGID, lenGID, tsize);

    hn->next = lu->ht[j];
    lu->ht[j] = hn;
   
    hn++;
    nextUniqueGID += lenGID;
    lu->numGIDs++;
  }

  return lu;
}
/****************************************************************************/
static int lookup_GID(GID_lookup *lu, ZOLTAN_ID_PTR gid)
{
  struct Hash_Node *hn;
  int i, k, match;

  if (lu->table_size < 1) return -1;
  if (lu->numGIDs < 1) return -1;

  i = Zoltan_Hash(gid, lu->lenGID, (unsigned int) lu->table_size);
  
  for (hn=lu->ht[i]; hn != NULL; hn = hn->next){
    match = 1;
    for (k=0; k<lu->lenGID; k++){
      if (hn->gid[k] != gid[k]){
        match = 0;
        break;
      }
    }
    if (match){
      return (hn->gno);
    }
  }
  return -1;
}

static int getObjectSizes(ZZ *zz, ZHG *zhg)
{
    int i, ierr=ZOLTAN_OK;
    char *yo="getObjectSizes";

    zhg->showMoveVol = 1;
    if (zhg->nObj) {
      if (!(zhg->AppObjSizes = (int *) ZOLTAN_MALLOC(zhg->nObj * sizeof(int)))) 
        MEMORY_ERROR;
      if (zz->Get_Obj_Size_Multi) {
        zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data,
                               zz->Num_GID, zz->Num_LID, zhg->nObj,
                               zhg->GIDs, zhg->LIDs, zhg->AppObjSizes, &ierr);
        if (ierr < 0) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                          "ZOLTAN_OBJ_SIZE_MULTI function.");
          goto End;
        }
      }
      else if (zz->Get_Obj_Size) {
        for (i = 0; i < zhg->nObj; i++) {
          ZOLTAN_ID_PTR lid = (zz->Num_LID ? &(zhg->LIDs[i*zz->Num_LID]):NULL);
          zhg->AppObjSizes[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data,
                                           zz->Num_GID, zz->Num_LID,
                                           &(zhg->GIDs[i*zz->Num_GID]),
                                           lid, &ierr);
          if (ierr < 0) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                            "ZOLTAN_OBJ_SIZE function.");
            goto End;
          }
        }
      }
      else {
        for (i = 0; i < zhg->nObj; i++) zhg->AppObjSizes[i] = 1;
      }
    }
 End:
    return ierr;
}

/****************************************************************************/
static int Zoltan_PHG_Add_Repart_Data(
  ZZ *zz,
  ZHG *zhg,
  HGraph *phg,              /* Input/Output:  Input hypergraph in 2D layout;
                               changed by this routine to include repartition 
                               vertices and repartition edges.  */
  int *vtxgno,              /* global vertex numbers for input vertices;
                               size is zhg->nObj. */
  PHGPartParams *hgp,       /* Input/Output:  Partitioning parameters; need
                               to turn on UseFixedVtx. */
  Partition input_parts     /* Input/Output:  Input partition assignments for 
                               vtxs in 2D distribution; entries for repartition
                               vertices are added at the end.  */
)
{
/* Function to add repartition vertices and edges to the hypergraph.
 * Repartition vertices and edges are distributed among the processor columns
 * and rows.
 * One repartition vertex is added for each part to be included in the 
 * partition (i.e., the target number of parts in the new partition).
 * One repartition edge is added for each vertex in one of the target parts.
 * Initially, edge memory is allocated large enough to hold the maximum
 * number of repartition edges; this eases the mapping from vertices to their
 * repartition edges when we create the edges.   The nRepartEdges variables
 * represent the maximum number of edges that would be added to a processor.
 * Then we remove zero-length edges resulting from vertices that are in parts
 * that will not be in the new partition (i.e., vertex v in part q with q > k.
 * We then readjust the nRepartEdges variables to represent the actual number
 * of repartition edges added to a processor and use the actual numbers to 
 * adjust phg->dist_y.
 */
char *yo = "Zoltan_PHG_Add_Repart_Data";
PHGComm *hgc = phg->comm;
int gnRepartVtx = zhg->GnRepartVtx;
int maxgnRepartEdge = zhg->GnObj;   /* max # of repartition edges that could
                                       be added to the global hypergraph. */
int myProc_x = hgc->myProc_x;
int myProc_y = hgc->myProc_y;
int nProc_x = hgc->nProc_x;
int nProc_y = hgc->nProc_y;
int nRepartVtx = 0, nRepartEdge;         /* # of repartition vertices & edges 
                                        to add to this processor. */
int nRepartPin = 0;                  /* # of repartition pins to add to this
                                        processor. */
int firstRepartVtx = 0, firstRepartEdge = 0; /* global no. of the first repartition
                                        vertex & edge to add to this processor.
                                        Repartition vertices and repartition
                                        edges are ordered and distributed
                                        linearly. */
int myStart_vtx, nextStart_vtx;      /* Each proc in column sends info about
                                        a subset of the column's vertices to
                                        procs owning the related
                                        repartition vertices.  myStart_vtx is
                                        the first vertex sent by a proc;
                                        nextStart_vtx is the first vertex
                                        sent by the next proc in the column. */

ZOLTAN_COMM_OBJ *plan;               /* Plan for communicating input part
                                        info to procs owning corresponding
                                        repartition vertices and edges. */
int *proclist = NULL;                /* Buffers to send/recv input part info */
int *sendbuf = NULL;                 /* to procs owning the corresponding   */
int *recvbuf = NULL;                 /* repartition vertices and edges.     */
int nrecv = 0;
                                 
int *tmp_hindex, *tmp_hvertex;       /* Pointers into phg->hindex and 
                                        phg->hvertex; set to start of entries
                                        for repartition edges. */
float *tmp_ewgt;                     /* Edge weights for local repartion
                                        edges */
int *pins_per_edge = NULL;           /* # of pins in each repartition edge. */
int i, j, idx;
int cnt, tmp, dim, sum, prev;
int ierr = ZOLTAN_OK;
int *objsize = NULL;                 /* repartition edge weights. */
int *tmpobjsize = NULL;
int *colProc_cnt = NULL;             /* actual nRepartEdge per column proc */
int *repart_dist_x = NULL;           /* Distribution of repartition vertices
                                        to proc rows; similar to phg->dist_x. */
int *repart_dist_y = NULL;           /* Distribution of repartition edges
                                        to proc cols; similar to phg->dist_y. */

  if (myProc_x >= 0 && myProc_y >= 0) {  /* This proc is part of 2D decomp */

    /* Compute distribution of repartition vertices and edges */
    repart_dist_x = (int *) ZOLTAN_MALLOC((2+nProc_x+nProc_y) * sizeof(int));
    repart_dist_y = repart_dist_x + nProc_x + 1;

    for (i = 0; i < nProc_x; i++)
      repart_dist_x[i] = FirstRepart(i, nProc_x, gnRepartVtx);
    repart_dist_x[nProc_x] = gnRepartVtx;

    for (i = 0; i < nProc_y; i++)
      repart_dist_y[i] = FirstRepart(i, nProc_y, maxgnRepartEdge);
    repart_dist_y[nProc_y] = maxgnRepartEdge;

    /* Compute number of repartition vertices to add to this processor column */
    phg->nRepartVtx = nRepartVtx
                    = repart_dist_x[myProc_x+1] - repart_dist_x[myProc_x];
    firstRepartVtx = repart_dist_x[myProc_x];

    /* Compute maximum number of repartition edges to add to this proc row */
    phg->nRepartEdge = nRepartEdge 
                     = repart_dist_y[myProc_y+1] - repart_dist_y[myProc_y];
    firstRepartEdge = repart_dist_y[myProc_y];

    /* If application did not fix any vertices, allocate fixed array for 
     * repartition vertices. 
     * Initialized application vertices' fixed value to -1 (not fixed).
     */
    hgp->UseFixedVtx = 1;
    if (!phg->fixed_part && (phg->nVtx + nRepartVtx)) {
      phg->fixed_part = (int *) ZOLTAN_MALLOC((phg->nVtx + nRepartVtx)
                                               * sizeof(int));
      if (!phg->fixed_part) MEMORY_ERROR;
      for (i = 0; i < phg->nVtx; i++) phg->fixed_part[i] = -1;
    }

    /* Assuming phg arrays are allocated large enough to accommodate the
       repartition vertices, edges, and pins.  */

    /***** Create repartition vertices in this processor column. *****/
    /***** Repartition vertices are distributed linearly.        *****/
    cnt = firstRepartVtx;
    dim = phg->VtxWeightDim;
    for (i = phg->nVtx; i < (phg->nVtx + nRepartVtx); i++) {
      phg->fixed_part[i] = cnt;
      input_parts[i] = cnt++;
      if (phg->vwgt) 
        for (j = 0; j < dim; j++)
          phg->vwgt[i*dim+j] = 0.;
    }

    /***** Initialize repartition edges in this processor row. *****/
    /***** Repartition edges are distributed linearly.         *****/
    dim = phg->EdgeWeightDim;
    for (i = phg->nEdge; i < (phg->nEdge + nRepartEdge); i++) {
      phg->hindex[i] = 0;  
      for (j = 0; j < dim; j++)
        phg->ewgt[i*dim+j] = 0.;
    }
  }  

  /***** Get object sizes to weight the repartition edges    *****/
  /***** with migration costs.                               *****/

  if (zz->Get_Obj_Size_Multi || zz->Get_Obj_Size) {
    if ((ierr = getObjectSizes(zz, zhg))!=ZOLTAN_OK)
        goto End;

    objsize = (int *) ZOLTAN_CALLOC(2*phg->nVtx, sizeof(int));
    if (phg->nVtx && !objsize) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    tmpobjsize = objsize + phg->nVtx;

    if (phg->comm->nProc_x > 1) {
      /* Send obj_size data to the 2D distribution */
      /* Use zhg->VtxPlan */
      ierr = Zoltan_Comm_Do(zhg->VtxPlan, 25232, (char *) zhg->AppObjSizes, 
                            sizeof(int), (char *) objsize);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
        goto End;
      }
      for (i = 0; i < zhg->nRecv_GNOs; i++) {
        idx = VTX_GNO_TO_LNO(phg, zhg->Recv_GNOs[i]);
        tmpobjsize[idx] = objsize[i];
      }
    }
    else {
      /* nProc_x == 1; have all the size data we need. */
      for (i = 0; i < zhg->nObj; i++)
        tmpobjsize[vtxgno[i]] = zhg->AppObjSizes[i];
    }

    /* Accrue each row's contribution to objsize */
    if (myProc_x >= 0 && myProc_y >= 0) {  /* This proc is part of 2D decomp */
      MPI_Allreduce(tmpobjsize, objsize, phg->nVtx, MPI_INT, MPI_SUM,
                    phg->comm->col_comm);
    }

/*    ZOLTAN_FREE(&appobjsize);*/
  }
  else {
    /* Use uniform object sizes. */
    objsize = NULL;
  }


  /***** Create repartition edges in this processor row. *****/
  
  /* Send objects' input partitions to processor owning corresponding
   * repartition edge and repartition vertex. 
   * Each proc in column sends data for a subset of the column's 
   * vertices.
   */

  if (myProc_x >= 0 && myProc_y >= 0) {  /* This proc is part of 2D decomp */

    myStart_vtx = (myProc_y * phg->nVtx) / nProc_y;
    nextStart_vtx = ((myProc_y+1) * phg->nVtx) / nProc_y;
    cnt = nextStart_vtx - myStart_vtx;
    tmp = 0;
  
#define NSEND 3  /* number of fields to be sent for each vertex */
    if (cnt) {
      proclist = (int *) ZOLTAN_MALLOC(2 * cnt * sizeof(int));
      sendbuf = (int *) ZOLTAN_MALLOC(NSEND * 2 * cnt * sizeof(int));
      if (!proclist || !sendbuf) MEMORY_ERROR;
  
      for (i = 0; i < cnt; i++){
        int vtxproc_x;            /* Processor column for actual vtx. */
        int proc_x;               /* Processor column for repartition vtx. */
        int proc_y;               /* Processor row for repartition edge. */
        int v = i + myStart_vtx;  /* Index of vertex being processed. */
        int vgno;                 /* Gno for vertex v. */
        int k = input_parts[v];   /* Repartition vtx global number */

        if (k < gnRepartVtx) {
          vgno = VTX_LNO_TO_GNO(phg,v);  /* Gno of vertex */

          /* Determine proc col for repartition vtx associated with part k */
          proc_x = ProcForRepart(k, repart_dist_x, nProc_x);

          /* Determine proc row for repartition edge associated with vtx i */
          proc_y = ProcForRepart(vgno, repart_dist_y, nProc_y);
 
          /* Fill message buffers */
  
          /* Send message to processor owning repartition vertex and edge */
          proclist[tmp] = proc_y * nProc_x + proc_x;
          sendbuf[NSEND*tmp] = vgno;            /* GNO of vertex */
          sendbuf[NSEND*tmp+1] = k;             /* Vertex's input partition */
          sendbuf[NSEND*tmp+2] = (objsize ? objsize[v] : 1);
          tmp++;

          /* Send message to proc owning actual vertex and repartition edge */
          vtxproc_x = VTX_TO_PROC_X(phg,vgno);
          if (vtxproc_x != proc_x) {
            proclist[tmp] = proc_y * nProc_x + vtxproc_x;
            sendbuf[NSEND*tmp] = vgno;            /* GNO of vertex */
            sendbuf[NSEND*tmp+1] = k;             /* Vertex's input partition */
            sendbuf[NSEND*tmp+2] = (objsize ? objsize[v] : 1);
            tmp++;
          }
        }
      }
    }

    ZOLTAN_FREE(&objsize);
    ZOLTAN_FREE(&repart_dist_x);

    /* Create plan and send info. */

    ierr = Zoltan_Comm_Create(&plan, tmp, proclist, hgc->Communicator,
                              24555, &nrecv);

    if (nrecv && !(recvbuf = (int *) ZOLTAN_MALLOC(NSEND*nrecv*sizeof(int))))
      MEMORY_ERROR;
  
    ierr = Zoltan_Comm_Do(plan, 24556, (char *) sendbuf, NSEND * sizeof(int),
                         (char *) recvbuf);

    ierr = Zoltan_Comm_Destroy(&plan);

    /* Build hindex & hvertex for repartition vertices & repartition edges.
     * Both arrays are appended to existing arrays, which are assumed to be
     * large enough.
     */

    tmp_hindex = &(phg->hindex[phg->nEdge]);
    tmp_hvertex = phg->hvertex;
  
    /* Loop to build hindex array */
    for (i = 0; i < nrecv; i++) {
      int vtx_gno;   /* Global vtx number of vtx in received repartition edge.*/
      int rEdge_lno; /* local index of repartition edge */
      int rVtx_gno;  /* global repartition vertex number */
      int rVtx_lno;  /* local index of repartition vertex */
  
      vtx_gno = recvbuf[NSEND*i];
      rEdge_lno = vtx_gno - firstRepartEdge;
      rVtx_gno = recvbuf[NSEND*i+1];
      rVtx_lno = rVtx_gno - firstRepartVtx;


      if (rVtx_gno >= firstRepartVtx && rVtx_gno < firstRepartVtx+nRepartVtx) {
        /* Add pin for repartition vertex */
        tmp_hindex[rEdge_lno]++;
        nRepartPin++;
      }
    
      if (myProc_x == VTX_TO_PROC_X(phg, vtx_gno)) {
        /* Add pin for application vertex */
        tmp_hindex[rEdge_lno]++;
        nRepartPin++;
      }
    }
    phg->nRepartPin = nRepartPin;

    if (nRepartEdge &&
        !(pins_per_edge = (int *) ZOLTAN_MALLOC(nRepartEdge * sizeof(int)))) 
      MEMORY_ERROR;
    MPI_Allreduce(tmp_hindex, pins_per_edge, nRepartEdge, MPI_INT, MPI_SUM,
                  phg->comm->row_comm);

    /* Compute prefix sum of tmp_hindex */
    prev = tmp_hindex[0];
    tmp_hindex[0] = phg->nPins;
    for (i = 1; i <= nRepartEdge; i++) {
      cnt = tmp_hindex[i];
      tmp_hindex[i] = tmp_hindex[i-1] + prev;
      prev = cnt;
    }
   
    /* Loop to build hvertex array and set edge weights */
    if (nRepartEdge &&
        !(tmp_ewgt = (float *) ZOLTAN_CALLOC(nRepartEdge * phg->EdgeWeightDim,
                                       sizeof(float))))
      MEMORY_ERROR;

    for (i = 0; i < nrecv; i++) {
      int vtx_gno;   /* Global vtx number of vtx in received repartition edge.*/
      int rEdge_lno; /* local index of repartition edge */
      int rVtx_gno;  /* global repartition vertex number */
      int rVtx_lno;  /* local index of repartition vertex */
  
      vtx_gno = recvbuf[NSEND*i];
      rEdge_lno = vtx_gno - firstRepartEdge;
      rVtx_gno = recvbuf[NSEND*i+1];
      rVtx_lno = rVtx_gno - firstRepartVtx;
      
      cnt = 0;
      if (rVtx_gno >= firstRepartVtx && rVtx_gno < firstRepartVtx+nRepartVtx) {
        /* Add pin for repartition vertex */
        tmp_hvertex[tmp_hindex[rEdge_lno]] = phg->nVtx + rVtx_lno;
        cnt++;
        for (j = 0; j < phg->EdgeWeightDim; j++) 
          tmp_ewgt[rEdge_lno*phg->EdgeWeightDim+j] = (float) recvbuf[NSEND*i+2];
      }
      if (myProc_x == VTX_TO_PROC_X(phg, vtx_gno)) {
        /* Add pin for vertex */
        tmp_hvertex[tmp_hindex[rEdge_lno]+cnt] = VTX_GNO_TO_LNO(phg, vtx_gno);
      }
    }
#undef NSEND

    /* Remove zero-length repartition edges */
    cnt = 0;
    for (i = 0; i < nRepartEdge; i++) {
      if (pins_per_edge[i] != 0) {
        if (cnt != i) {
          tmp_hindex[cnt] = tmp_hindex[i];
          for (j = 0; j < phg->EdgeWeightDim; j++)
            tmp_ewgt[cnt*phg->EdgeWeightDim+j] = 
                     tmp_ewgt[i*phg->EdgeWeightDim+j];
        }
        cnt++;
      }
    }
    tmp_hindex[cnt] = tmp_hindex[nRepartEdge];
    phg->nRepartEdge = nRepartEdge = cnt;
    ZOLTAN_FREE(&pins_per_edge);
  
    /* Accrue edge weights across processor rows */
    MPI_Allreduce(tmp_ewgt, &phg->ewgt[phg->nEdge*phg->EdgeWeightDim],
                  nRepartEdge*phg->EdgeWeightDim, MPI_FLOAT, MPI_SUM,
                  phg->comm->row_comm);
    ZOLTAN_FREE(&tmp_ewgt);

    /* Modify hypergraph edge weights to account for the number of 
     * communications done between migrations (the RepartMultiplier).
     * Do not apply the multiplier to the repartition edges.
     */
    for (i = 0; i < phg->nEdge; i++)
      phg->ewgt[i] *= hgp->RepartMultiplier;

    /* Now update hypergraph info to include the repartition vtx & edges. */
    /* We'll subtract these back out after we finish partitioning and before
     * we compute final output and build return lists. */

    phg->nVtx += nRepartVtx;
    phg->nEdge += nRepartEdge;
    phg->nPins += nRepartPin;

    /* SANITY CHECK */
    if (phg->nPins != phg->hindex[phg->nEdge]) {
      uprintf(phg->comm, 
             "%d KDDKDD SANITY CHECK FAILED %d != %d   %d %d %d %d %d %d\n", 
             zz->Proc, phg->nPins, phg->hindex[phg->nEdge], nRepartVtx, 
             phg->nVtx, nRepartEdge, phg->nEdge, nRepartPin, phg->nPins);
      exit( -1);
    }

    /* Update distx and disty to include the repartition vtx & edges. */
    /* Easy for repartition vertices. */
    for (sum = 0, i = 0; i < nProc_x; i++) {
      phg->dist_x[i] += sum;
      sum += NumRepart(i, nProc_x, gnRepartVtx);
    }
    phg->dist_x[nProc_x] += sum;

    /* Tougher for repartition edges, since we removed the empty ones.
     * Need to get actual nRepartEdge values from each proc in column. */
    
    colProc_cnt = (int *) ZOLTAN_MALLOC(nProc_y * sizeof(int));
    MPI_Allgather(&nRepartEdge, 1, MPI_INT, colProc_cnt, 1, MPI_INT, 
                  phg->comm->col_comm);
    for (sum = 0, i = 0; i < nProc_y; i++) {
      phg->dist_y[i] += sum;
      sum += colProc_cnt[i];
    }
    phg->dist_y[nProc_y] += sum;
    ZOLTAN_FREE(&colProc_cnt);
  }

End:

  ZOLTAN_FREE(&proclist);
  ZOLTAN_FREE(&sendbuf);
  ZOLTAN_FREE(&recvbuf);
  return ierr;
}
/****************************************************************************/
int Zoltan_PHG_Remove_Repart_Data(
  ZZ *zz,
  ZHG *zhg,
  HGraph *phg,
  PHGPartParams *hgp
)
{
/* Routine to remove repartition vertices, repartition edges, and 
 * repartition pins after partitioning with REPARTITION.
 * We don't have to actually remove the data; we only have to update
 * certain hypergraph fields so subsequent processing calls don't
 * look at the data.
 */
int nProc_x = phg->comm->nProc_x;
int nProc_y = phg->comm->nProc_y;
int sum;
int i;
int *colProc_cnt = NULL;
int myProc_x = phg->comm->myProc_x;
int myProc_y = phg->comm->myProc_y;

  if (myProc_x >= 0 && myProc_y >= 0) {  /* This proc is part of 2D decomp */
    /* Only processors in the 2D decomp added repartition data, so 
     * only those processors need to remove it.
     */
    phg->nVtx  -= phg->nRepartVtx;
    phg->nEdge -= phg->nRepartEdge;
    phg->nPins -= phg->nRepartPin;

    if (zz->Edge_Weight_Dim) {
      /* Remove the RepartMultiplier from edge weights */
      for (i = 0; i < phg->nEdge; i++)
        phg->ewgt[i] /= hgp->RepartMultiplier;
    }
    else {
      /* Application didn't provide edge weights; we should remove them */
      ZOLTAN_FREE(&phg->ewgt);
      phg->EdgeWeightDim = 0;
    }

    /* Update distx and disty to remove the repartition vtx & edges. */
    for (sum = 0, i = 0; i < nProc_x; i++) {
      phg->dist_x[i] -= sum;
      sum += NumRepart(i, nProc_x, zhg->GnRepartVtx);
    }
    phg->dist_x[nProc_x] -= sum;
  
    colProc_cnt = (int *) ZOLTAN_MALLOC(nProc_y * sizeof(int));
    MPI_Allgather(&(phg->nRepartEdge), 1, MPI_INT, colProc_cnt, 1, MPI_INT, 
                  phg->comm->col_comm);
    for (sum = 0, i = 0; i < nProc_y; i++) {
      phg->dist_y[i] -= sum;
      sum += colProc_cnt[i];
    }
    phg->dist_y[nProc_y] -= sum;
    ZOLTAN_FREE(&colProc_cnt);
  }

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
