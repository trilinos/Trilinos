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
#include "parmetis_jostle.h"
#include "zz_util_const.h"

/*#define DEBUG_FILL_HYPERGRAPH 1*/
    
#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
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

/*****************************************************************************/
/* Function prototypes */

static int hash_lookup(ZZ*, ZOLTAN_ID_PTR, int, struct Hash_Node**, struct Hash_Node **);
static int get_vertex_global_numbers(ZZ *zz, ZHG *zhg,
  int nVtx, struct Hash_Node **ht,
  int nPins, ZOLTAN_ID_PTR pins, int *pin_gno, int *pin_procs);
static int resolve_edge_weight_contributions( ZZ *zz, int ew_op,
  int *ew_num_edges, ZOLTAN_ID_PTR ew_gids, ZOLTAN_ID_PTR ew_lids,
  float *ew_weights, void **ht, int *htsize);
static int combine_weights_for_same_edge(ZZ *zz, int ew_op,
          int *ew_num_edges, 
          ZOLTAN_ID_PTR ew_gids, ZOLTAN_ID_PTR ew_lids, 
          float *ew_weights, void **ht, int *htsize);
static int combine_weights(int ew_op, float *dest, float *src, int n);
#ifdef DEBUG_FILL_HYPERGRAPH
static void print_hypergraph(ZZ *zz, ZHG *zhg, int sumWeight);  /* for debugging */
#endif

/*****************************************************************************/

int Zoltan_PHG_Build_Hypergraph(
  ZZ *zz,                            /* Zoltan data structure */
  ZHG **zoltan_hg,                  /* Hypergraph to be allocated and built.*/
  Partition *input_parts,            /* Initial partition assignments for
                                        vtxs in 2D distribution; length = 
                                        zoltan_hg->HG->nVtx.  */
  PHGPartParams *hgp                 /* Parameters for PHG partitioning.*/
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
  zhg->nRemove = 0;
  zhg->Remove_EGIDs = NULL;
  zhg->Remove_ELIDs = NULL;
  zhg->Remove_Esize = NULL;
  zhg->Remove_Ewgt = NULL;
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
  }
    
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/


int Zoltan_PHG_Fill_Hypergraph(
  ZZ *zz,
  ZHG *zhg,      /* Description of hypergraph provided by the application. */
  PHGPartParams *hgp,                /* Parameters for PHG partitioning.*/
  Partition *input_parts   /* Initial partition assignment of vtxs in 
                              2D data distribution; length = zhg->HG->nVtx. */
)
{
/* Routine to call HG query function and build HG data structure. 
 * Output is a fully functioning parallel hypergraph with 2D distribution of
 * pins (non-zeros).
 */

char *yo = "Zoltan_PHG_Fill_Hypergraph";
struct application_input {     /* Data provided by hypergraph callbacks. */
  int nVtx;                         /* # objects (vertices) on proc. */
  int nEdge;                        /* # hyperedges on proc. */
  int nPins;                        /* # pins (nonzeros) on proc. */
  int GnVtx;                        /* Total nVtx across all procs. */
  int GnEdge;                       /* Total nEdge across all procs. */
  ZOLTAN_ID_PTR egids;              /* GIDs of edges. */
  ZOLTAN_ID_PTR elids;              /* LIDs of edges. */
  int *esizes;                      /* # of GIDs in each hyperedge. */
  ZOLTAN_ID_PTR pins;               /* Object GIDs (vertices) belonging to 
                                       hyperedges.  */
  int *pin_procs;                   /* Processor owning each pin vertex of 
                                       hyperedges. */
  int *pin_gno;                     /* Global numbers in range [0,GnVtx-1]
                                       for pins. */
  int *vtx_gno;                     /* Global numbers in range [0,GnVtx-1]
                                       for nVtx vertices.  app.vtx_gno[i] is
                                       global number for GID[i]. */
  int *edge_gno;                    /* Global numbers in range [0,GnEdge-1]
                                       for edges.  app.edge[i] is
                                       global number for this proc's edge i. */
  float *vwgt;              /* Vertex weights for nVtx on processor objects. */
  float *ewgt;                      /* Edge weights. */
} app;

struct Hash_Node *hash_nodes = NULL;  /* Hash table variables for mapping   */
struct Hash_Node **hash_tab = NULL;   /* GIDs to global numbering system.   */
struct Hash_Node *hn;
ZOLTAN_COMM_OBJ *plan;

int i, j, cnt, dim, rc;
int msg_tag = 30000;
int ierr = ZOLTAN_OK;
int nProc = zz->Num_Proc;
int nRequests;
ZOLTAN_ID_PTR pin_requests = NULL;
int *request_gno = NULL;
int edge_gno, edge_Proc_y;
int vtx_gno, vtx_Proc_x;
int *proclist = NULL;
int *sendbuf = NULL;
int nnz, idx;
int *nonzeros = NULL;
int *tmparray = NULL;
int *hindex = NULL, *hvertex = NULL;
int *dist_x = NULL, *dist_y = NULL;
int nEdge, nVtx, nwgt = 0;
int nrecv, *recv_gno = NULL; 
int *tmpparts = NULL;
int add_vweight;
int ew_num_edges, ew_ht_size;
int ew_dim = zz->Edge_Weight_Dim;
int hypergraph_callbacks = 0;
int graph_callbacks = 0;
void *ew_ht;
float *ew_weights;
ZOLTAN_ID_PTR global_ids, ew_gids, ew_lids;
int num_gid_entries = zz->Num_GID;
HGraph *phg = &(zhg->HG);
int nProc_x, nProc_y, myProc_x, myProc_y;
int proc_offset;

float frac_x, frac_y;
float *tmpwgts=NULL; 
float *fromwgt, *towgt, *wgts, *calcVwgt;

int *gtotal = NULL;  /* Temporary arrays used for distributing vertices/edges */
int *mycnt = NULL;
int *gcnt = NULL;

float edgeSizeThreshold;
int randomizeInitDist, edge_weight_op, zoltan_lb_eval;
int final_output, add_obj_weight;
PHGPartParams *temphgp = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (hgp){
    randomizeInitDist = hgp->RandomizeInitDist;
    edge_weight_op = hgp->edge_weight_op;
    edgeSizeThreshold = hgp->EdgeSizeThreshold;
    final_output = hgp->final_output;
    add_obj_weight = hgp->add_obj_weight;
    zoltan_lb_eval = 0;
  }
  else{
    /*   
     * hgp and input_parts are undefined when we are called from 
     * Zoltan_LB_Eval.  (This only happens in the pins callback case.)
     * In this case, we want all pins to be in the removed list.
     * And we don't care how global numbers are assigned.
     */
    temphgp = (PHGPartParams *)ZOLTAN_MALLOC(sizeof(PHGPartParams));
    Zoltan_PHG_Initialize_Params(zz, NULL, temphgp);

    randomizeInitDist = 0;   /* faster */
    edge_weight_op = temphgp->edge_weight_op;
    edgeSizeThreshold = 0;   /* place all edges in "removed" list */
    final_output = 1;        /* yes, compile a "removed" list     */
    add_obj_weight = temphgp->add_obj_weight;

    ZOLTAN_FREE(&temphgp);

    zoltan_lb_eval = 1;
  }
  if (zz->Get_HG_Size_CS && zz->Get_HG_CS){
    hypergraph_callbacks = 1;
  }
  else if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
           (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) {
    graph_callbacks = 1;
  }

  /**************************************************/
  /* Obtain vertex information from the application */
  /**************************************************/

  app.egids = NULL;
  app.elids = NULL;
  app.esizes = NULL;
  app.pins = NULL;
  app.pin_procs = NULL;
  app.pin_gno = NULL;
  app.vtx_gno = NULL;
  app.edge_gno = NULL;
  app.vwgt = NULL;
  app.ewgt = NULL;

  ierr = Zoltan_Get_Obj_List(zz, &(zhg->nObj), &(zhg->GIDs), &(zhg->LIDs), 
                             zz->Obj_Weight_Dim, &app.vwgt,
                             &(zhg->Input_Parts));
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
    goto End;
  }

  app.nVtx = zhg->nObj;

  /* 
   * Correlate GIDs in edge_verts with local indexing in zhg to build the
   * input HG.
   * Use hash table to map global IDs to local position in zhg->GIDs.
   * Based on hashing code in Zoltan_Build_Graph.
   */

  /* Construct local hash table mapping GIDs to global number (gno) */
  if (app.nVtx) {
    hash_nodes = (struct Hash_Node *) ZOLTAN_MALLOC(app.nVtx * 
                                                    sizeof(struct Hash_Node));
    hash_tab = (struct Hash_Node **) ZOLTAN_MALLOC(app.nVtx *
                                                   sizeof(struct Hash_Node *));
    app.vtx_gno = (int *) ZOLTAN_MALLOC(app.nVtx * sizeof(int));
    if (!hash_nodes || !hash_tab || !app.vtx_gno) MEMORY_ERROR;
  }

  global_ids = zhg->GIDs;

  /* Assign consecutive numbers (gnos) based on the order of the ids */
  if (randomizeInitDist) { 
    /* Randomize the input vertices */
    int tmp;
    gtotal = (int *) ZOLTAN_CALLOC(3*zz->Num_Proc+1, sizeof(int));
    mycnt  = gtotal + zz->Num_Proc + 1;
    gcnt   = mycnt + zz->Num_Proc;

    /* Compute random processor bin. */
    /* Temporarily store processor bin number in app.vtx_gno. */
    /* Count how many local vtxs selected processor bin */
    Zoltan_Srand(Zoltan_Rand(NULL)+zz->Proc, NULL);
    for (i = 0; i < app.nVtx; i++) {
      app.vtx_gno[i] = Zoltan_Rand_InRange(NULL, zz->Num_Proc);
      mycnt[app.vtx_gno[i]]++;
    }
    /* Compute prefix of mycnt */
    rc = MPI_Scan(mycnt, gcnt, zz->Num_Proc, MPI_INT, MPI_SUM, zz->Communicator);
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
    app.GnVtx = gtotal[zz->Num_Proc] = tmp;

    /* Assign gnos sequential from gcnt[bin]. */
    for (i=0; i< app.nVtx; i++) {
      hash_tab[i] = NULL;
      hash_nodes[i].gid = &(global_ids[i*num_gid_entries]);
      tmp = app.vtx_gno[i];
      hash_nodes[i].gno = app.vtx_gno[i] = gtotal[tmp] + gcnt[tmp];
      gcnt[tmp]++;
    }
  }
  else {
    /* Linearly order the input vertices */
    gtotal = (int *) ZOLTAN_MALLOC((nProc+1) * sizeof(int));
    if (!gtotal) MEMORY_ERROR;

    /* Construct gtotal[i] = the number of vertices on all procs < i. */
    /* Scan to compute partial sums of the number of objs */

    rc = MPI_Scan (&app.nVtx, gtotal, 1, MPI_INT, MPI_SUM, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Gather data from all procs */

    rc = MPI_Allgather (&(gtotal[0]), 1, MPI_INT,
                   &(gtotal[1]), 1, MPI_INT, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)
    gtotal[0] = 0;
    app.GnVtx = gtotal[nProc];

    for (i=0; i< app.nVtx; i++) {
      hash_tab[i] = NULL;
      hash_nodes[i].gid = &(global_ids[i*num_gid_entries]);
      hash_nodes[i].gno = app.vtx_gno[i] =  gtotal[zz->Proc]+i;
    }
  }

  for (i=0; i< app.nVtx; i++){
    /* insert hashed elements into hash table */
    j = Zoltan_Hash(&(global_ids[i*num_gid_entries]), num_gid_entries,
                    (unsigned int) app.nVtx);
    hash_nodes[i].next = hash_tab[j];
    hash_tab[j] = &hash_nodes[i];
  }

  /***********************************************************************/
  /* Get hyperedge information from application through query functions. */
  /***********************************************************************/

  app.nEdge = 0;
  app.nPins = 0;
  app.GnEdge = 0;

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
     */

    ew_gids = NULL;
    ew_lids = NULL;
    ew_weights = NULL;
    ew_num_edges = 0;
    ew_ht = NULL;

    if (ew_dim && zz->Get_HG_Size_Edge_Weights && zz->Get_HG_Edge_Weights){

      ierr = zz->Get_HG_Size_Edge_Weights(
                   zz->Get_HG_Size_Edge_Weights_Data, &ew_num_edges);

      if (((ierr==ZOLTAN_OK)||(ierr==ZOLTAN_WARN)) && (ew_num_edges > 0)){
        ew_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, ew_num_edges);
        ew_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, ew_num_edges);
        ew_weights =
          (float *)ZOLTAN_MALLOC(sizeof(float) * ew_num_edges * ew_dim);

        if (!ew_gids || !ew_lids || !ew_weights){
          Zoltan_Multifree(__FILE__,__LINE__,3,&ew_gids,&ew_lids,&ew_weights);
          MEMORY_ERROR;
        }

        ierr = zz->Get_HG_Edge_Weights(zz->Get_HG_Edge_Weights_Data,
                    zz->Num_GID, zz->Num_LID, ew_num_edges, ew_dim,
                    ew_gids, ew_lids, ew_weights);


      }

      if ((ierr==ZOLTAN_OK)||(ierr==ZOLTAN_WARN)){
        /* 
         * Global operation: merge all weight contributions for
         * the same hyperedge using the PHG_EDGE_WEIGHT_OPERATION
         * parameter.  Create a hash table for obtaining the weight
         * of a particular local hyperedge.
         *
         * Call may rewrite with shorter lists if edges appear more than
         * once in the list.  So this call may change ew_num_edges.
         */
        ierr = resolve_edge_weight_contributions(zz, edge_weight_op,
               &ew_num_edges, ew_gids, ew_lids, ew_weights, 
               &ew_ht, &ew_ht_size);
      }

      rc = (((ierr == ZOLTAN_OK) || (ierr == ZOLTAN_WARN)) ? 0 : 1);

      MPI_Allreduce(&rc, &cnt, 1, MPI_INT, MPI_MAX, zz->Communicator);

      if (cnt && ((ierr == ZOLTAN_OK) || (ierr == ZOLTAN_WARN))){
        ierr = ZOLTAN_FATAL;   /* another process had an error */
      }

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        Zoltan_Multifree(__FILE__, __LINE__, 3,
                         &ew_gids, &ew_lids, &ew_weights);

        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
          "Error querying/processing edge weights");
        goto End;
      }
    }
    /*
     * Call the compressed storage pin callbacks.  Redistribute the
     * the pins so each process has complete rows (edges) and no
     * two processes have the same row.  If edge weights were
     * provided above, obtain edge weights for the edges owned
     * by the process.
     *
     * This call uses the ew_ht and then frees it.
     */
    ierr = Zoltan_HG_Hypergraph_Pin_Callbacks(zz, zhg,
                                app.GnVtx, edgeSizeThreshold,
                                final_output, ew_num_edges, ew_ht, ew_ht_size,
                                &app.nEdge, &app.egids,
                                &app.elids,            /* optional */
                                &app.esizes, &app.ewgt,
                                &app.nPins, &app.pins);

    ZOLTAN_FREE(&ew_gids);
    ZOLTAN_FREE(&ew_lids);
    ZOLTAN_FREE(&ew_weights);

    } else if (graph_callbacks){

    /* Graph query functions, one hyperedge per vertex */

    ierr = Zoltan_HG_Graph_Callbacks(zz, zhg,
                                     app.GnVtx, edgeSizeThreshold,
                                     final_output, &app.nEdge, 
                                     &app.egids, &app.elids,
                                     &app.esizes, &app.ewgt,
                                     &app.nPins, &app.pins, 
                                     &app.pin_procs);

  } else {
    /* Partition without edge information?  Or return an error? */
    ZOLTAN_PRINT_WARN(zz->Proc, yo,
       "No edge information provided, partitioning vertices.");
  }

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Callbacks");
    goto End;
  }

  /***********************************************************************/
  /* Impose a global hyperedge numbering */
  /***********************************************************************/

  app.edge_gno = (int *) ZOLTAN_MALLOC(app.nEdge * sizeof(int));
  app.pin_gno = (int *) ZOLTAN_MALLOC(app.nPins * sizeof(int));
  if ((app.nEdge && !app.edge_gno) || (app.nPins && !app.pin_gno)) MEMORY_ERROR;

  if (randomizeInitDist) {  
    /* Randomize the input edges */
    int tmp;
    memset(gtotal, 0, (3*zz->Num_Proc+1)*sizeof(int)); /* gtotal was alloc'ed
                                                          for vertices */

    /* Compute random processor bin. */
    /* Temporarily store processor bin number in app.edge_gno. */
    /* Count how many local vtxs selected processor bin */
    for (i = 0; i < app.nEdge; i++) {
      app.edge_gno[i] = Zoltan_Rand_InRange(NULL, zz->Num_Proc);
      mycnt[app.edge_gno[i]]++;
    }
    /* Compute prefix of mycnt */
    rc = MPI_Scan(mycnt, gcnt, zz->Num_Proc, MPI_INT, MPI_SUM, zz->Communicator);
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
    app.GnEdge = gtotal[zz->Num_Proc] = tmp;

    /* Assign gnos sequential from gcnt[bin]. */
    for (i=0; i< app.nEdge; i++) {
      tmp = app.edge_gno[i];
      app.edge_gno[i] = gtotal[tmp] + gcnt[tmp];
      gcnt[tmp]++;
    }
  }
  else {
    rc = MPI_Scan (&app.nEdge, gtotal, 1, MPI_INT, MPI_SUM, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Gather data from all procs */

    rc = MPI_Allgather (&(gtotal[0]), 1, MPI_INT,
                   &(gtotal[1]), 1, MPI_INT, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)
    gtotal[0] = 0;
    app.GnEdge = gtotal[nProc];

    /* Assign global numbers to edges. */
    for (i = 0; i < app.nEdge; i++)
      app.edge_gno[i] = gtotal[zz->Proc] + i;
  }
  ZOLTAN_FREE(&gtotal);
   
  /***********************************************************************/
  /* 
   * Obtain the global num in range 0 .. (total_num_vtx-1) 
   * for each vertex pin.
   */
  /***********************************************************************/

  if (hypergraph_callbacks){
    /*
     * Determine which process owns the vertex for each of my
     * pins (write it to app.pin_procs), and while doing so
     * obtain the vertex global number for each pin.  (In the case
     * of graph callbacks, the pin_procs were provided by the
     * application in the callback function.)
     */
    app.pin_procs = (int *) ZOLTAN_MALLOC(app.nPins * sizeof(int));

    ierr = get_vertex_global_numbers(zz, zhg, app.nVtx, hash_tab,
             app.nPins, app.pins, app.pin_gno, app.pin_procs);

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting vertex global numbers");
      goto End;
    }
  }

  /*
   * Create a plan to communicate with all owners of my pin vertices.
   *
   * Graph Callbacks case:
   * For each pin GID in app.pins, request the global number (gno) from the
   * input processor.  (Hypergraph callbacks got this already.)  
   *
   * Graph and Hypergraph callbacks case:
   * If ADD_OBJ_WEIGHT indicated that a vertex weight equal to the number
   * of pins should be added, then obtain the sum of pins for each vertex now.
   */

  if (add_obj_weight != PHG_ADD_NO_WEIGHT){
    add_vweight = 1;    /* may change in future to allow more than one */
    calcVwgt = (float *)ZOLTAN_CALLOC(sizeof(float), app.nVtx);

    if (add_obj_weight == PHG_ADD_UNIT_WEIGHT){
      for (i=0; i<app.nVtx; i++){
        calcVwgt[i] = 1.0;
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

  if (phg->VtxWeightDim > 1){
    /*
     * For now, only one weight per vertex will be used.  We will
     * still save multiple weights, because in the future we may
     * be able to use them.
     */
    if (zz->Proc == 0){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Too many vertex weights.");
      if (add_vweight){
        ZOLTAN_PRINT_WARN(zz->Proc, yo, 
         "Both application supplied *and* ADD_OBJ_WEIGHT "
         "calculated vertex weights were provided.");
      }
      else{
        ZOLTAN_PRINT_WARN(zz->Proc, yo, 
         "Multiple weights per vertex were supplied.");
      }
      ZOLTAN_PRINT_WARN(zz->Proc, yo, 
         "Only the first application supplied weight per vertex will be used.");
    }
  }

  if ((add_obj_weight == PHG_ADD_PINS_WEIGHT) || graph_callbacks){
  
    ierr = Zoltan_Comm_Create(&plan, app.nPins, app.pin_procs, zz->Communicator,
                              msg_tag, &nRequests);
  
    if (nRequests) {
      pin_requests = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);

      if (graph_callbacks){
        request_gno = (int *) ZOLTAN_MALLOC(nRequests * sizeof(int));
      }
      if (!pin_requests || (graph_callbacks && !request_gno))
        MEMORY_ERROR;
    }
  
    msg_tag--;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) app.pins, 
                          sizeof(ZOLTAN_ID_TYPE) * num_gid_entries,
                          (char *) pin_requests);
  
    for (i = 0; i < nRequests; i++){
      j = hash_lookup(zz, &(pin_requests[i*num_gid_entries]), 
                      app.nVtx, hash_tab, &hn);

      if (!hn){ 
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "invalid hash table");
        goto End;
      }
      
      if (add_obj_weight == PHG_ADD_PINS_WEIGHT){
        idx = hn - hash_nodes;
        calcVwgt[idx] = calcVwgt[idx] + 1.0;
      }
      if (graph_callbacks){
        request_gno[i] = j;
      }
    }
  
    if (graph_callbacks){
      msg_tag--;
      Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *) request_gno, sizeof(int), NULL,
                             (char *) app.pin_gno);
    }
  
    Zoltan_Comm_Destroy(&plan);
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
      wgts = (float *)ZOLTAN_MALLOC(sizeof(float) * app.nVtx * phg->VtxWeightDim);
      if (!wgts){ 
        MEMORY_ERROR
      }
  
      towgt = wgts;
      fromwgt = app.vwgt;
  
      for (i=0; i<app.nVtx; i++){
        for (j=0; j < zz->Obj_Weight_Dim; j++){
          *towgt++ = *fromwgt++;
        }
        *towgt++ = calcVwgt[i];
      }
      ZOLTAN_FREE(&calcVwgt);
    }
    else{
      wgts = calcVwgt;
    }

    ZOLTAN_FREE(&app.vwgt);

    app.vwgt = wgts;
  }

  if (zoltan_lb_eval){
    /* 
     * We were called from Zoltan_LB_Eval and we're done.
     * All edges, pins, pin owners, etc are in the "removed" lists.
     * We write vertex weights for vertices owned by this process 
     * to the hypergraph structure (where pin weights normally are).
     */

    phg->vwgt = app.vwgt;
    app.vwgt = NULL;

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

  frac_x = (float) app.GnVtx / (float) nProc_x;
  for (i = 1; i < nProc_x; i++)
    dist_x[i] = (int) (i * frac_x);
  dist_x[nProc_x] = app.GnVtx;
  
  frac_y = (float) app.GnEdge / (float) nProc_y;
  for (i = 1; i < nProc_y; i++)
    dist_y[i] = (int) (i * frac_y);
  dist_y[nProc_y] = app.GnEdge;
  
  /* myProc_y and myProc_x can be -1 when nProc is prime and we use a 2D
   * decomposition.  One processor is excluded from the 2D communicator;
   * for it, myProc_y and myProc_x == -1. */
  nEdge = (myProc_y >= 0 ? dist_y[myProc_y+1] - dist_y[myProc_y] : 0);
  nVtx  = (myProc_x >= 0 ? dist_x[myProc_x+1] - dist_x[myProc_x] : 0);

  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution. 
   */

  proclist = (int *) ZOLTAN_MALLOC(MAX(app.nPins,app.nVtx) * sizeof(int));
  sendbuf = (int *) ZOLTAN_MALLOC(app.nPins * 2 * sizeof(int));

  cnt = 0; 
  for (i = 0; i < app.nEdge; i++) {
    /* processor row for the edge */
    edge_gno = app.edge_gno[i];
    edge_Proc_y = EDGE_TO_PROC_Y(phg, edge_gno);

    for (j = 0; j < app.esizes[i]; j++) {
      /* processor column for the vertex */
      vtx_gno = app.pin_gno[cnt];
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
  hindex = (int *) ZOLTAN_CALLOC(nEdge + 1, sizeof(int));
  hvertex = (int *) ZOLTAN_MALLOC(nnz * sizeof(int));

  if (!tmparray || !hindex || (nnz && !hvertex)) MEMORY_ERROR;

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
    hvertex[hindex[idx] + tmparray[idx]] = VTX_GNO_TO_LNO(phg, nonzeros[2*i+1]);
    tmparray[idx]++;
  }

  phg->nVtx = nVtx;
  phg->nEdge = nEdge;
  phg->nPins = nnz;
  phg->hindex = hindex;
  phg->hvertex = hvertex;

  ierr = Zoltan_HG_Create_Mirror(zz, phg);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error from Zoltan_HG_Create_Mirror");
    goto End;
  }

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
    for (i = 0; i < app.nVtx; i++)
      proclist[i] = proc_offset + VTX_TO_PROC_X(phg, app.vtx_gno[i]);
      
    msg_tag++;
    ierr = Zoltan_Comm_Create(&(zhg->VtxPlan), app.nVtx, proclist, 
                              zz->Communicator, msg_tag, &nrecv);
    zhg->nRecv_GNOs = nrecv;

    zhg->Recv_GNOs = recv_gno = (int *) ZOLTAN_MALLOC(nrecv * sizeof(int));

    if (nrecv && !recv_gno) MEMORY_ERROR;

    /* Use plan to send global numbers to the appropriate proc_x. */
    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) app.vtx_gno, 
                          sizeof(int), (char *) recv_gno);
  }
  else {
    /* Save map of what needed. */
    zhg->nRecv_GNOs = nrecv = app.nVtx;
    zhg->Recv_GNOs = recv_gno = app.vtx_gno;
  }

  /* Send vertex partition assignments and weights to 2D distribution. */

  tmpparts = (int *) ZOLTAN_CALLOC(phg->nVtx, sizeof(int));
  *input_parts = (int *) ZOLTAN_MALLOC(phg->nVtx * sizeof(int));

  dim = phg->VtxWeightDim;
  nwgt = phg->nVtx * dim;

  tmpwgts = (float *)ZOLTAN_CALLOC(sizeof(float), nwgt);
  phg->vwgt = (float *)ZOLTAN_MALLOC(sizeof(float) * nwgt);

  if (phg->nVtx && 
       (!tmpparts || !*input_parts || !tmpwgts || !phg->vwgt)) MEMORY_ERROR;

  if (phg->comm->nProc_x == 1)  {
    for (i = 0; i < app.nVtx; i++) {
      idx = app.vtx_gno[i];
      tmpparts[idx] = zhg->Input_Parts[i];
      for (j=0; j<dim; j++){
        tmpwgts[idx*dim + j] = app.vwgt[i*dim + j];
      }
    }
  }
  else {
    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) zhg->Input_Parts,
                          sizeof(int), (char *) *input_parts);

    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) app.vwgt,
                          sizeof(float) * dim, (char *) phg->vwgt);

    for (i = 0; i < nrecv; i++) {
      idx = VTX_GNO_TO_LNO(phg, recv_gno[i]);
      tmpparts[idx] = (*input_parts)[i];
      for (j=0; j<dim; j++){
        tmpwgts[idx*dim + j] = phg->vwgt[i*dim + j];
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
  }

  ZOLTAN_FREE(&tmpparts);
  ZOLTAN_FREE(&tmpwgts);

  if (!zz->LB.Remap_Flag && zz->LB.Return_Lists == ZOLTAN_LB_NO_LISTS) {
    int gnremove;
    rc = MPI_Allreduce(&(zhg->nRemove), &gnremove, 1, MPI_INT, MPI_SUM, 
                  zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)
    if (!final_output || !gnremove) {
      /* Don't need the plan long-term; destroy it now. */
      Zoltan_Comm_Destroy(&(zhg->VtxPlan));
      if (zhg->Recv_GNOs == app.vtx_gno) app.vtx_gno = NULL;
      ZOLTAN_FREE(&(zhg->Recv_GNOs));
      zhg->nRecv_GNOs = 0;
    }
  }

  /*  Send edge weights, if any */

  dim = phg->EdgeWeightDim = zz->Edge_Weight_Dim;
  nwgt = phg->nEdge * dim;

  if (zz->Edge_Weight_Dim) {
    tmpwgts   = (float *) ZOLTAN_CALLOC(nwgt , sizeof(float));
    phg->ewgt = (float *) ZOLTAN_CALLOC(nwgt, sizeof(float));

    if (nwgt && (!phg->ewgt || !tmpwgts)) MEMORY_ERROR;

    if (phg->comm->nProc_y == 1) {
      for (i = 0; i < app.nEdge; i++) {
        idx = app.edge_gno[i];
        for (j = 0; j < dim; j++)
          tmpwgts[idx * dim + j] = app.ewgt[i*dim + j];
      }
    }
    else {
      /* 
       * Since for 2D decomposition and prime nProc we exclude a processor,
       * we cannot use col_comm for this operation.  We'll simulate it,
       * allowing the excluded processor to send its data to col 0.
       */
      proc_offset = (myProc_x >= 0 ? myProc_x : 0);
      for (i = 0; i < app.nEdge; i++)
        proclist[i] = proc_offset 
                    + EDGE_TO_PROC_Y(phg, app.edge_gno[i]) * nProc_x;
      
      msg_tag++;

      ierr = Zoltan_Comm_Create(&plan, app.nEdge, proclist, 
                                zz->Communicator, msg_tag, &nrecv); 

      recv_gno = (int *) ZOLTAN_MALLOC(nrecv * sizeof(int));
      if (nrecv && !recv_gno) MEMORY_ERROR;

      msg_tag++;
      ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) app.edge_gno, sizeof(int),
                            (char *) recv_gno);

      /* In using phg->ewgt for recv, assuming nrecv <= phg->nVtx */
      msg_tag++;
      ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) app.ewgt, 
                            dim*sizeof(float), (char *) phg->ewgt);

      Zoltan_Comm_Destroy(&plan);

      for (i = 0; i < nrecv; i++) {
        idx = EDGE_GNO_TO_LNO(phg, recv_gno[i]);
        for (j = 0; j < dim; j++) 
          tmpwgts[idx * dim + j] = phg->ewgt[i * dim + j];
      }
      ZOLTAN_FREE(&recv_gno); 
    }

    /* Need to gather weights for all edges within row 
     * to all processors within row.
     */

    if (phg->comm->row_comm != MPI_COMM_NULL)
      rc = MPI_Allreduce(tmpwgts, phg->ewgt, nwgt, MPI_FLOAT, MPI_MAX, 
                    phg->comm->row_comm);
      CHECK_FOR_MPI_ERROR(rc)
  }
  else {
    /* KDDKDD  For now, do not assume uniform edge weights.
     * KDDKDD  Can add later if, e.g., we decide to coalesce identical edges.
     */
    phg->EdgeWeightDim = 0;
  }

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    Zoltan_PHG_Free_Hypergraph_Data(zhg);
  }
  
  Zoltan_Multifree(__FILE__, __LINE__, 11, &app.egids,
                                           &app.elids,
                                           &app.pins, 
                                           &app.esizes, 
                                           &app.pin_gno, 
                                           &app.pin_procs, 
                                           &app.vwgt,
                                           &app.ewgt,
                                           &app.edge_gno,
                                           &hash_nodes,
                                           &hash_tab);

  if (zhg->Recv_GNOs != app.vtx_gno) 
    ZOLTAN_FREE(&app.vtx_gno);
  ZOLTAN_FREE(&tmparray);
  ZOLTAN_FREE(&tmpwgts);
  ZOLTAN_FREE(&nonzeros);
  ZOLTAN_FREE(&proclist);
  ZOLTAN_FREE(&sendbuf);
  ZOLTAN_FREE(&pin_requests);
  ZOLTAN_FREE(&request_gno);

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

int Zoltan_PHG_Removed_Cuts(
  ZZ *zz,
  ZHG *zhg,
  double *localcuts  /* Array of length 2. Upon return:
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
int i, j, k, cnt, ncnt, nparts;
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
int msg_tag = 23132;
double ewgt;

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
  pin_parts = pin_procs;    /* Reuse the memory allocated for pin_procs. */
  Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *) outparts, 
                         sizeof(int), NULL, (char *) pin_parts);

  ZOLTAN_FREE(&outparts);

  Zoltan_Comm_Destroy(&plan);

  /* Compute the cut metrics using received partition info. */

  parts = (int *) ZOLTAN_CALLOC(zz->LB.Num_Global_Parts, sizeof(int));
  if (!parts) MEMORY_ERROR;

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
  
  localcuts[0] = loccuts[0];
  localcuts[1] = loccuts[1];

End:

  if (!zhg->Remove_Pin_Procs){
    ZOLTAN_FREE(&pin_procs);
  }
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
/*****************************************************************************/
static int create_gno_list(ZZ *zz, ZOLTAN_ID_PTR rcvBufGids, int nVtx,
              struct Hash_Node **ht, int *sndBufGnos);
static int apply_new_gnos(ZZ *zz,
         ZOLTAN_ID_PTR gids, int *pinIdx, int *new_gnos, int *pin_gnos,
         int *pin_procs, int *remove_pin_procs, int npins, int rproc);

static int get_vertex_global_numbers(ZZ *zz, ZHG *zhg,
  int nVtx, struct Hash_Node **ht,
  int nPins, ZOLTAN_ID_PTR pins, int *pin_gno, int *pin_procs)
{
char *yo = "get_vertex_global_numbers";
int i, gno, maxUnSet, unSet, ngnos, ngids;
ZOLTAN_ID_PTR p, rcvBufGids=NULL, sndBufGids=NULL;
int *sndBufGnos=NULL, *rcvBufGnos=NULL, *pinIdx=NULL;
int rank, nprocs, send_to_proc, recv_from_proc;
int rcvBufGidSize, rcvBufGnoSize, sndBufGidSize;
int gids_tag, gnos_tag, rc;
int nRemovePins, totalPins;
MPI_Request gids_req, gnos_req;
MPI_Status status;
int ierr=ZOLTAN_OK;


  /*  
   * First, determine whether there are any removed pins.  For
   * these, we only need to know which process owns the vertex.
   */

  nRemovePins = 0;
  for (i=0; i<zhg->nRemove; i++){
    nRemovePins += zhg->Remove_Esize[i];
  }

  if (nRemovePins){
    ZOLTAN_FREE(&zhg->Remove_Pin_Procs);
    zhg->Remove_Pin_Procs = (int *)ZOLTAN_MALLOC(nRemovePins * sizeof(int));

    /* Note which vertices are mine */
    for (i=0, p=zhg->Remove_Pin_GIDs ; i<nRemovePins; i++, p += zz->Num_GID){
      gno = hash_lookup(zz, p, nVtx, ht, NULL);
      if (gno >= 0){
        zhg->Remove_Pin_Procs[i] = zz->Proc;
      }
      else{
        zhg->Remove_Pin_Procs[i] = -1;
      }
    }
  }

  /* For pins in our edges, we need to contact the owning process
   * to find out the global number.  Create message buffers (include
   * removed pins if we have any) and determine which of the
   * pin vertices we own.
   */

  totalPins = nPins + nRemovePins;

  pinIdx = (int *)ZOLTAN_MALLOC(totalPins * sizeof(int));

  for (i=0, p=pins ; i<totalPins; i++){
    pinIdx[i] = i; 

    if (i  < nPins){
      gno = hash_lookup(zz, p, nVtx, ht, NULL);
   
      if (gno >= 0){
        pin_procs[i] = zz->Proc;
        pin_gno[i] = gno;
      }
      else{
        pin_procs[i] = -1;
        pin_gno[i] = -1;
      }
      p += zz->Num_GID;
    }
  }
  sndBufGids = ZOLTAN_MALLOC_GID_ARRAY(zz, totalPins + 1);

  sndBufGids[0] = totalPins;
  ZOLTAN_COPY_GID_ARRAY(sndBufGids + zz->Num_GID, pins, zz, nPins);

  if (nRemovePins){
    ZOLTAN_COPY_GID_ARRAY(sndBufGids + (zz->Num_GID * (1 + nPins)),
                        zhg->Remove_Pin_GIDs, zz, nRemovePins);
  }

  /* Rewite sndBufGids with only those GIDs for which we need gnos/owners
   * and update pinIdx to map them to the arrays.
   */

  unSet = apply_new_gnos(zz, sndBufGids, pinIdx, NULL, pin_gno, 
                         pin_procs, zhg->Remove_Pin_Procs, nPins, 0);

  rc = MPI_Allreduce(&unSet, &maxUnSet, 1, MPI_INT, MPI_MAX, zz->Communicator);
  CHECK_FOR_MPI_ERROR(rc)

  if (maxUnSet == 0){
    ZOLTAN_FREE(&pinIdx);
    ZOLTAN_FREE(&sndBufGids);
    return ZOLTAN_OK;    /* all processes have all global numbers */
  }

  /* Exchange lists of vertex GIDS with other processes, each of
   * which supplies global numbers for the vertices it owns.
   */

  rcvBufGids = ZOLTAN_MALLOC_GID_ARRAY(zz, maxUnSet + 1);
  sndBufGnos = (int *)ZOLTAN_MALLOC(((2 * maxUnSet)+1) * sizeof(int));
  rcvBufGnos = (int *)ZOLTAN_MALLOC(((2 * unSet)+1) * sizeof(int));

  rank = zz->Proc;
  nprocs = zz->Num_Proc;
  rcvBufGidSize = (maxUnSet + 1) * zz->Num_GID;
  sndBufGidSize = (unSet + 1) * zz->Num_GID;
  rcvBufGnoSize = 2 * unSet + 1;
  gids_tag = 0x0101;
  gnos_tag = 0x1111;

  for (i=1; i < nprocs; i++){

    send_to_proc = (rank + i) % nprocs;
    recv_from_proc = (rank + nprocs - i) % nprocs;

    /* post a recv for neighbor to left for gids to translate */

    rc = MPI_Irecv(rcvBufGids, rcvBufGidSize, MPI_INT, recv_from_proc,
               gids_tag, zz->Communicator, &gids_req);
    CHECK_FOR_MPI_ERROR(rc)

    /* send my gids to neighbor to right, and post a receive for
     * the translated gnos                                      */

    rc = MPI_Send(sndBufGids, sndBufGidSize, MPI_INT, send_to_proc,
               gids_tag, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    rc = MPI_Irecv(rcvBufGnos, rcvBufGnoSize, MPI_INT, send_to_proc,
               gnos_tag, zz->Communicator, &gnos_req);
    CHECK_FOR_MPI_ERROR(rc)

    /* Await the gids from my left, translate the ones
       I own, and send back the global numbers                 */

    rc = MPI_Wait(&gids_req, &status);
    CHECK_FOR_MPI_ERROR(rc)

    ngnos = create_gno_list(zz, rcvBufGids, nVtx, ht, sndBufGnos);

    rc = MPI_Send(sndBufGnos, (ngnos*2)+1, MPI_INT, recv_from_proc,
             gnos_tag, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Await the translated GIDs from my right, add them to
       my pin_gno list, create a new smaller list of
       untranslated GIDs.                                   */

    rc = MPI_Wait(&gnos_req, &status);
    CHECK_FOR_MPI_ERROR(rc)

    ngids = apply_new_gnos(zz, sndBufGids, pinIdx, rcvBufGnos, pin_gno, 
                  pin_procs, zhg->Remove_Pin_Procs, nPins, send_to_proc);

    sndBufGidSize = (ngids + 1) * zz->Num_GID;
    rcvBufGnoSize = 2 * ngids + 1;
  }

End:

  ZOLTAN_FREE(&pinIdx);
  ZOLTAN_FREE(&sndBufGids);
  ZOLTAN_FREE(&rcvBufGids);
  ZOLTAN_FREE(&sndBufGnos);
  ZOLTAN_FREE(&rcvBufGnos);

  if (ngids > 0){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Failed to obtain all vertex global IDs.")
    return ZOLTAN_FATAL;
  }

  return ierr;
}

static int apply_new_gnos(ZZ *zz,
         ZOLTAN_ID_PTR gids, int *pinIdx, int *new_gnos, int *pin_gnos,
         int *pin_procs, int *remove_pin_procs, int npins, int rproc)
{
int nTranslated, nUntranslated;
int ngids, i, which, gno;
int *rbuf = NULL;
ZOLTAN_ID_PTR p1, p2;

  if (new_gnos){
    nTranslated = new_gnos[0];
    rbuf = new_gnos + 1;
  }
  else{
    nTranslated = 0;
  }

  /* Add vertex global numbers received from remote process */

  for (i=0; i<nTranslated; i++){
    which = rbuf[0];
    gno = rbuf[1];

    if (pinIdx[which] < npins){
      pin_gnos[pinIdx[which]] = gno;
      pin_procs[pinIdx[which]] = rproc;
    }
    else{
      remove_pin_procs[pinIdx[which] - npins] = rproc;
    }
    rbuf += 2;
  }

  /* Rewrite Gid list to contain only Gids we need owner for
   * and update pinIdx list to map them to the arrays.
   */

  ngids = gids[0];
  p1 = p2 = gids + zz->Num_GID;
  nUntranslated = 0;

  for (i=0; i<ngids; i++){

    if ( ((pinIdx[i] < npins) && (pin_procs[pinIdx[i]] < 0)) ||
         ((pinIdx[i] >= npins) && (remove_pin_procs[pinIdx[i] - npins] < 0))){

      if (p1 < p2){
        ZOLTAN_SET_GID(zz, p1, p2);
        pinIdx[nUntranslated] = pinIdx[i];
      }
      p1 += zz->Num_GID;
      nUntranslated++;
    }
    p2 += zz->Num_GID;
  }

  gids[0] = nUntranslated;

  return nUntranslated;
}
static int create_gno_list(ZZ *zz, ZOLTAN_ID_PTR rcvBufGids, int nVtx,
              struct Hash_Node **ht, int *sndBufGnos)
{
int numIds, numGnos, i, gno;
int *buf;
ZOLTAN_ID_PTR p;

  numGnos = 0;

  if ((nVtx > 0) && rcvBufGids && (rcvBufGids[0] > 0)){

    numIds = rcvBufGids[0];
    p = rcvBufGids + zz->Num_GID;
    buf = sndBufGnos + 1;

    for (i=0; i<numIds; i++){
      gno = hash_lookup(zz, p, nVtx, ht, NULL);
      if (gno >= 0){
        *buf++ = i;
        *buf++ = gno;
        numGnos++;
      }
      p += zz->Num_GID;
    }
  }
  sndBufGnos[0] = numGnos;

  return numGnos;
}
/*****************************************************************************/
static int resolve_edge_weight_contributions(
         ZZ *zz, int ew_op,
         int *ew_num_edges,
         ZOLTAN_ID_PTR ew_gids,
         ZOLTAN_ID_PTR ew_lids,
         float *ew_weights,  /* rewrite weights if necessary */
         void **ht,          /* returns search structure for edges */
         int *htsize)        /* size of hash table */
{
  /* Multiple processes may have supplied weights for the same   */
  /* edge.  We resolve these in a manner depending on the        */
  /* PHG_EDGE_WEIGHT_OPERATION parameter.  We also create and    */
  /* return a hash table useful for looking up edge weights.      */

char *yo = "resolve_edge_weight_contributions";
ZOLTAN_ID_PTR rcvBufGids=NULL, sndBufGids=NULL, gidptr;
float *rcvBufWeights=NULL, *sndBufWeights=NULL, *wptr;
int rank, nprocs, right_proc, left_proc, rc;
int rcvBufGidSize, sndBufGidSize, sndBufWeightSize;
int gids_tag, weights_tag, match_tag;
int match_left, match_right;
MPI_Request gids_req, weights_req;
MPI_Status status;
struct _ewht{
  ZOLTAN_ID_PTR egid;
  ZOLTAN_ID_PTR elid;
  float *weights;
  struct _ewht *next;
} *ewNodes=NULL, *en;
struct _ewht **ewht = NULL;
int i, j, k, ierr, nEdgesMax, nids, first_match=0;
int nedges, ew_ht_size;
int dim = zz->Edge_Weight_Dim;
int lenGID = zz->Num_GID;

  *ht = NULL;
  ierr = ZOLTAN_OK;

  rc = MPI_Allreduce(ew_num_edges, &nEdgesMax, 1, MPI_INT, MPI_MAX,
                 zz->Communicator);
  CHECK_FOR_MPI_ERROR(rc)

  if (nEdgesMax <= 0){
    return ZOLTAN_OK;
  }

  /* If this process listed the same gid more than once, combine
   * any weights supplied for the same gid.  (zdrive may do this
   * when reading a "matrixmarket plus" file created by n processes
   * into a zdrive application of fewer than n processes.)
   * 
   * Also create a search structure to find a weight given a gid.
   */
  
  rc = combine_weights_for_same_edge(zz, ew_op,
          ew_num_edges, ew_gids, ew_lids, ew_weights, ht, htsize);

  if (rc != ZOLTAN_OK){
    ierr = rc;
    goto End;
  }

  ewht = (struct _ewht **)*ht;
  ew_ht_size = *htsize;
  nedges = *ew_num_edges;

  if (ewht){
    ewNodes = ewht[ew_ht_size];  /* stored pointer here */
  }
  else{
    ewNodes = NULL;
  }

  /*
   * Do pairwise exchanges with all other processes to match weights
   * provided for the same edge.
   */

  sndBufWeightSize = nedges*dim;
  sndBufWeights = (float *)ZOLTAN_MALLOC(sndBufWeightSize * sizeof(float));
  rcvBufGids = ZOLTAN_MALLOC_GID_ARRAY(zz, nEdgesMax+1);
  sndBufGids = ZOLTAN_MALLOC_GID_ARRAY(zz, nedges+1);
  rcvBufWeights = (float *)ZOLTAN_MALLOC(nEdgesMax*dim*sizeof(float));

  if ((sndBufWeightSize && !sndBufWeights) ||
      !rcvBufGids || !sndBufGids || !rcvBufWeights){

    Zoltan_Multifree(__FILE__, __LINE__, 6,
       &sndBufWeights, &rcvBufGids, &sndBufGids, &rcvBufWeights,
       &ewNodes, &ewht);

    return ZOLTAN_MEMERR;
  }

  rcvBufGidSize = lenGID * (nEdgesMax+1);
  sndBufGidSize = lenGID * (nedges+1);
  sndBufGids[0] = nedges;
  ZOLTAN_COPY_GID_ARRAY(sndBufGids+lenGID, ew_gids, zz, nedges);

  if (sndBufWeightSize){
    memcpy(sndBufWeights, ew_weights, sndBufWeightSize * sizeof(float));
  }

  rank = zz->Proc;
  nprocs = zz->Num_Proc;
  gids_tag = 0x2222;
  weights_tag = 0x4444;
  match_tag = 0x1111;

  for (i=1; i < nprocs; i++){
    right_proc = (rank + i) % nprocs;
    left_proc = (rank + nprocs - i) % nprocs;

    /* post a recv for neighbor to left for edge gids */

    rc = MPI_Irecv(rcvBufGids, rcvBufGidSize, MPI_INT, left_proc,
               gids_tag, zz->Communicator, &gids_req);
    CHECK_FOR_MPI_ERROR(rc)

    /* send my gids to neighbor to right */

    rc = MPI_Send(sndBufGids, sndBufGidSize, MPI_INT, right_proc,
               gids_tag, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Await edge gids from my left.  If we have any edges in
     * common, post a receive for the weights and tell
     * proc on the left to send the weights to me.
     */

    rc = MPI_Wait(&gids_req, &status);
    CHECK_FOR_MPI_ERROR(rc)

    match_left = 0;
    nids = rcvBufGids[0];
    gidptr = rcvBufGids + lenGID;

    for (k=0; (k<nids) && nedges && !match_left; k++){
      j = Zoltan_Hash(gidptr, lenGID, ew_ht_size);
      en = ewht[j];
      while (en && !match_left){
        if (ZOLTAN_EQ_GID(zz, gidptr, en->egid)){
          match_left = 1;
          first_match = k;
        }
        en = en->next;
      }
      gidptr += lenGID;
    }

    if (match_left){
        rc = MPI_Irecv(rcvBufWeights, nids*dim, MPI_FLOAT, left_proc,
             weights_tag, zz->Communicator, &weights_req);
        CHECK_FOR_MPI_ERROR(rc)
    }

    rc = MPI_Send(&match_left, 1, MPI_INT, left_proc, match_tag, 
                  zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc)

    /* Await message from right, indicating whether it needs
     * my edge weights, and send them if so.
     */

    rc = MPI_Recv(&match_right, 1, MPI_INT, right_proc, match_tag, 
                  zz->Communicator,
     &status);
    CHECK_FOR_MPI_ERROR(rc)

    if (match_right){
      rc = MPI_Send(sndBufWeights, sndBufWeightSize, MPI_FLOAT, right_proc,
               weights_tag, zz->Communicator);
      CHECK_FOR_MPI_ERROR(rc)
    }

    /* Await edge weights from left (if I requested them) and
     * resolve common weights (add, take max, or flag error if different).
     */
    if (match_left){
      rc = MPI_Wait(&weights_req, &status);
      CHECK_FOR_MPI_ERROR(rc)

      gidptr = rcvBufGids + (first_match + 1)*lenGID;
      wptr = rcvBufWeights + (first_match * dim);

      for (k=first_match; (k<nids) && (ierr!=ZOLTAN_FATAL); k++){
        j = Zoltan_Hash(gidptr, lenGID, ew_ht_size);
        en = ewht[j];
        while (en){
          if (ZOLTAN_EQ_GID(zz, gidptr, en->egid)){
            ierr = combine_weights(ew_op, en->weights, wptr, dim);
            if (ierr == ZOLTAN_FATAL){  
              ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Inconsistent edge weights");
              /* finish nprocs loop so application doesn't hang */
            }
            break;
          }
          en = en->next;
        }
        gidptr += lenGID;
        wptr += dim;
      }
    }
  }
End:

  Zoltan_Multifree(__FILE__, __LINE__, 4,
       &sndBufWeights, &rcvBufGids, &sndBufGids, &rcvBufWeights);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    ZOLTAN_FREE(&ewNodes);
    ZOLTAN_FREE(&ewht);
  }

  return ierr;
}
static int combine_weights_for_same_edge(ZZ *zz, int ew_op,
          int *ew_num_edges, 
          ZOLTAN_ID_PTR ew_gids, ZOLTAN_ID_PTR ew_lids, 
          float *ew_weights, void **ht, int *htsize)
{
char *yo = "combine_weights_for_same_edge";
int ierr, nedges, edge_count, i, j, w;
struct _ewht{
  ZOLTAN_ID_PTR egid;
  ZOLTAN_ID_PTR elid;
  float *weights;
  struct _ewht *next;
} *ewNodes=NULL, *en;
struct _ewht **ewht = NULL;
int dim = zz->Edge_Weight_Dim;
int lenGID = zz->Num_GID;
int lenLID = zz->Num_LID;
ZOLTAN_ID_PTR rgid, rlid, wgid, wlid;
float *rwgt, *wwgt;

  *ht = NULL;
  *htsize = 0;

  if (*ew_num_edges < 1){
    return ZOLTAN_OK;
  }

  /* I have a list of edge global IDs, and weights for each edge.
   * Create a search structure to locate an edge weight(s) given
   * it's global ID.  In the process, if the same global ID is
   * listed more than once, combine the weights for the edge according
   * to the PHG_EDGE_WEIGHT_OPERATION, and only list the edge once.
   */

  nedges = *ew_num_edges;
  edge_count = 0;

  ewNodes =
    (struct _ewht *)ZOLTAN_MALLOC(nedges * sizeof(struct _ewht));
  ewht =
    (struct _ewht **)ZOLTAN_CALLOC(nedges + 1, sizeof(struct _ewht *));

  if (!ewNodes || !ewht){
    ZOLTAN_FREE(&ewNodes);
    ZOLTAN_FREE(&ewht);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "memory allocation");
    return ZOLTAN_MEMERR;
  }

  rgid = wgid = ew_gids;
  rlid = wlid = ew_lids;
  rwgt = wwgt = ew_weights;

  for (i=0; i<nedges; i++){
    j = Zoltan_Hash(rgid, lenGID, (unsigned int)nedges);

    en = ewht[j];

    while (en){
      if (ZOLTAN_EQ_GID(zz, en->egid, rgid)){
        ierr = combine_weights(ew_op, en->weights, rwgt, dim);
        if (ierr != ZOLTAN_OK){
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Inconsistent edge weights");
          return ZOLTAN_FATAL;
        }
        break;
      }
      en = en->next;
    }

    if (en == NULL){   /* edge was not a duplicate */

      if (edge_count < i){
        ZOLTAN_SET_GID(zz, wgid, rgid);
        ZOLTAN_SET_LID(zz, wlid, rlid);
        for (w=0; w<dim; w++){
          wwgt[w] = rwgt[w];
        }
      }

      ewNodes[edge_count].egid = wgid;
      ewNodes[edge_count].elid = wlid;
      ewNodes[edge_count].weights = wwgt;
      ewNodes[edge_count].next = ewht[j];

      ewht[j] = ewNodes + edge_count;

      edge_count++;

      wgid += lenGID;
      wlid += lenLID;
      wwgt += dim;
    }

    rgid += lenGID;
    rlid += lenLID;
    rwgt += dim;
  }

  ewht[nedges] = ewNodes;  /* save pointer so it can be free'd later */

  *ew_num_edges = edge_count;
  *ht = (void *)ewht;
  *htsize = nedges;

  return ZOLTAN_OK;
}
static int combine_weights(int ew_op, float *dest, float *src, int n)
{
  int w;

  if (ew_op == PHG_FLAG_ERROR_EDGE_WEIGHTS){
    for (w=0; w<n; w++){
      if (src[w] != dest[w]){
        return ZOLTAN_FATAL;
      }
    }
  } else if (ew_op == PHG_MAX_EDGE_WEIGHTS){
    for (w=0; w<n; w++){
      if (src[w] > dest[w]){
        dest[w] = src[w];
      }
    }
  } else if (ew_op == PHG_ADD_EDGE_WEIGHTS){
    for (w=0; w<n; w++){
      dest[w] += src[w];
    }
  }
  else{
    return ZOLTAN_FATAL;
  }
  return ZOLTAN_OK;
}

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

  printf("REMOVED EDGES (%d): (gid lid size pins/owners %d weights)\n",
                  zhg->nRemove, ewdim);

  pin = zhg->Remove_Pin_GIDs;
  owner = zhg->Remove_Pin_Procs;
  wgt = zhg->Remove_Ewgt;
  nremovedpins = 0;

  for (i=0; i<zhg->nRemove; i++){
    printf("  %d:  %d  %d  %d  ",i, zhg->Remove_EGIDs[i], zhg->Remove_ELIDs[i],
                zhg->Remove_Esize[i]);
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


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
