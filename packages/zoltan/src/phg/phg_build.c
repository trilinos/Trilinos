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
    
#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}

/*****************************************************************************/
/* Function prototypes */

static int hash_lookup(ZZ*, ZOLTAN_ID_PTR, int, struct Hash_Node**);
static int Zoltan_PHG_Fill_Hypergraph(ZZ*, ZHG*, PHGPartParams *, Partition*);

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
    Zoltan_HG_HGraph_Free(&(zhg->HG));
    Zoltan_Multifree(__FILE__, __LINE__, 4, &(zhg->GIDs),
     &(zhg->LIDs), &(zhg->Input_Parts), zoltan_hg);
  }
    
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/


static int Zoltan_PHG_Fill_Hypergraph(
  ZZ *zz,
  ZHG *zhg,      /* Description of hypergraph provided by the application. */
  PHGPartParams *hgp,                /* Parameters for PHG partitioning.*/
  Partition *input_parts   /* Initial partition assignment of vtxs in 
                              2D data distribution; length = zhg->HG->nVtx. */
)
{
/* Routine to call HG query function and build HG data structure. 
 * Assumes (for now) that input is given as follows from application:
 *  - Application gives hyperedges it owns; it knows GIDs and processor owner
 *    of all vertices of each owned hyperedge.
 *  - Application gives each hyperedge only once (KDDKDD -- MAY REMOVE THIS
 *    RESTRICTION LATER. KDDKDD)
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
                                       for vertices.  app.vtx_gno[i] is
                                       global number for GID[i]. */
  int *edge_gno;                    /* Global numbers in range [0,GnEdge-1]
                                       for edges.  app.edge[i] is
                                       global number for this proc's edge i. */
  float *vwgt;                      /* Vertex weights. */
  float *ewgt;                      /* Edge weights. */
} app;

struct Hash_Node *hash_nodes = NULL;  /* Hash table variables for mapping   */
struct Hash_Node **hash_tab = NULL;   /* GIDs to global numbering system.   */
ZOLTAN_COMM_OBJ *plan;

int i, j, cnt, dim;
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
int add_vweight, nvwgt;

ZOLTAN_ID_PTR global_ids;
int num_gid_entries = zz->Num_GID;
HGraph *phg = &(zhg->HG);

int nProc_x = phg->comm->nProc_x;
int nProc_y = phg->comm->nProc_y;
int myProc_x = phg->comm->myProc_x;
int myProc_y = phg->comm->myProc_y;
int proc_offset;

float frac_x, frac_y;
float *tmpwgts = NULL;
float *tmp, *final;

int *gtotal = NULL;  /* Temporary arrays used for distributing vertices/edges */
int *mycnt = NULL;
int *gcnt = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

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
  if (hgp->RandomizeInitDist) { 
    /* Randomize the input vertices */
    int tmp;
    gtotal = (int *) ZOLTAN_CALLOC(3*zz->Num_Proc+1, sizeof(int));
    mycnt  = gtotal + zz->Num_Proc + 1;
    gcnt   = mycnt + zz->Num_Proc;

    /* Compute random processor bin. */
    /* Temporarily store processor bin number in app.vtx_gno. */
    /* Count how many local vtxs selected processor bin */
    Zoltan_Srand(zz->Proc, NULL);
    for (i = 0; i < app.nVtx; i++) {
      app.vtx_gno[i] = Zoltan_Rand_InRange(NULL, zz->Num_Proc);
      mycnt[app.vtx_gno[i]]++;
    }
    /* Compute prefix of mycnt */
    MPI_Scan(mycnt, gcnt, zz->Num_Proc, MPI_INT, MPI_SUM, zz->Communicator);
    MPI_Allreduce(mycnt, gtotal, zz->Num_Proc, MPI_INT, MPI_SUM, 
                  zz->Communicator);

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

    MPI_Scan (&app.nVtx, gtotal, 1, MPI_INT, MPI_SUM, zz->Communicator);

    /* Gather data from all procs */

    MPI_Allgather (&(gtotal[0]), 1, MPI_INT,
                   &(gtotal[1]), 1, MPI_INT, zz->Communicator);
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
  if (zz->Get_Num_HG_Edges && zz->Get_HG_Edge_List && zz->Get_HG_Edge_Info) {
    ierr = Zoltan_HG_Hypergraph_Edge_Callbacks(zz, zhg,
                                          app.GnVtx, hgp->EdgeSizeThreshold,
                                          hgp->final_output, &app.nEdge, 
                                          &app.egids, &app.elids,
                                          &app.esizes, &app.ewgt,
                                          &app.nPins, &app.pins, 
                                          &app.pin_procs);
  } else if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
             (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) {
    ierr = Zoltan_HG_Graph_Callbacks(zz, zhg,
                                     app.GnVtx, hgp->EdgeSizeThreshold,
                                     hgp->final_output, &app.nEdge, 
                                     &app.egids, &app.elids,
                                     &app.esizes, &app.ewgt,
                                     &app.nPins, &app.pins, 
                                     &app.pin_procs);

  } else {
    /* Partition without edge information?  Or return an error? */
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

  if (hgp->RandomizeInitDist) {  
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
    MPI_Scan(mycnt, gcnt, zz->Num_Proc, MPI_INT, MPI_SUM, zz->Communicator);
    MPI_Allreduce(mycnt, gtotal, zz->Num_Proc, MPI_INT, MPI_SUM, 
                  zz->Communicator);

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
    MPI_Scan (&app.nEdge, gtotal, 1, MPI_INT, MPI_SUM, zz->Communicator);

    /* Gather data from all procs */

    MPI_Allgather (&(gtotal[0]), 1, MPI_INT,
                   &(gtotal[1]), 1, MPI_INT, zz->Communicator);
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
   * For each pin GID in app.pins, request the global number (gno) from the
   * input processor.
   * Fill requests (using hash table) for GIDs local to this processor.
   * Upon completion, app.pin_gno will contain the global nums.
   */
  /***********************************************************************/

  ierr = Zoltan_Comm_Create(&plan, app.nPins, app.pin_procs, zz->Communicator,
                            msg_tag, &nRequests);

  if (nRequests) {
    pin_requests = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
    request_gno = (int *) ZOLTAN_MALLOC(nRequests * sizeof(int));
    if (!pin_requests || !request_gno)
      MEMORY_ERROR;
  }

  msg_tag--;
  ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) app.pins, 
                        sizeof(ZOLTAN_ID_TYPE) * num_gid_entries,
                        (char *) pin_requests);

  for (i = 0; i < nRequests; i++)
    request_gno[i] = hash_lookup(zz, &(pin_requests[i*num_gid_entries]),
                                 app.nVtx, hash_tab);

  msg_tag--;
  Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *) request_gno, sizeof(int), NULL,
                         (char *) app.pin_gno);

  Zoltan_Comm_Destroy(&plan);

  /* 
   * Compute the distribution of vertices and edges to the 2D data
   * distribution's processor columns and rows.
   * For now, these distributions are described by arrays dist_x
   * and dist_y; in the future, we may prefer a hashing function
   * mapping GIDs to processor columns and rows. KDDKDD
   */

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

    /* Use plan to send weights to the appropriate proc_x. */
    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) app.vtx_gno, 
                          sizeof(int), (char *) recv_gno);

  }
  else {
    /* Save map of what needed. */
    zhg->nRecv_GNOs = nrecv = app.nVtx;
    zhg->Recv_GNOs = recv_gno = app.vtx_gno;
  }

  /* Send vertex partition assignments to 2D distribution. */

  tmpparts = (int *) ZOLTAN_CALLOC(phg->nVtx, sizeof(int));
  *input_parts = (int *) ZOLTAN_MALLOC(phg->nVtx * sizeof(int));
  if (phg->nVtx && (!tmpparts || !*input_parts)) MEMORY_ERROR;

  if (phg->comm->nProc_x == 1)  {
    for (i = 0; i < app.nVtx; i++) {
      idx = app.vtx_gno[i];
      tmpparts[idx] = zhg->Input_Parts[i];
    }
  }
  else {
    /* In using input_parts for recv, assuming nrecv <= phg->nVtx */
    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) zhg->Input_Parts,
                          sizeof(int), (char *) *input_parts);

    for (i = 0; i < nrecv; i++) {
      idx = VTX_GNO_TO_LNO(phg, recv_gno[i]);
      tmpparts[idx] = (*input_parts)[i];
    }
  }

  /* Reduce partition assignments for all vertices within column 
   * to all processors within column.
   */

  if (phg->comm->col_comm != MPI_COMM_NULL)
    MPI_Allreduce(tmpparts, *input_parts, phg->nVtx, MPI_INT, MPI_MAX, 
                  phg->comm->col_comm);

  ZOLTAN_FREE(&tmpparts);

  /* Create arrays of vertex and edge weights */

  if (hgp->add_obj_weight != PHG_ADD_NO_WEIGHT){
    add_vweight = 1;
  }
  else{
    add_vweight = 0;
  }

  phg->VtxWeightDim =     /* 1 or greater */
    zz->Obj_Weight_Dim +   /* 0 or more application supplied weights*/
    add_vweight;           /* 0 or 1 additional calculated weights */

  nvwgt = phg->nVtx * phg->VtxWeightDim;

  nwgt = MAX(nvwgt, phg->nEdge * zz->Edge_Weight_Dim);

  tmpwgts   = (float *) ZOLTAN_CALLOC(nwgt , sizeof(float));
  phg->vwgt = (float *) ZOLTAN_CALLOC(nvwgt, sizeof(float));

  if (!tmpwgts || !phg->vwgt) MEMORY_ERROR;

  if ((phg->comm->col_comm == MPI_COMM_NULL) ||
      (zz->Obj_Weight_Dim == 0)) {
    tmp = tmpwgts;
    final = phg->vwgt;
  }
  else {
    tmp = phg->vwgt;
    final = tmpwgts;
  }

  /* Send vertex weights to 2D distribution. */

  if (add_vweight) {
    /* Either application did not specify object weights, so we
     * are creating them, or application asked for calculated
     * weights with ADD_OBJ_WEIGHT parameter. 
     * Create uniform weights, or use #pins.
     */
    cnt = zz->Obj_Weight_Dim + 1;
 
    if (hgp->add_obj_weight == PHG_ADD_PINS_WEIGHT){ /* vertex degree */

      /* sum local degrees to global degrees and use as vtx weight */
      if (phg->comm->col_comm != MPI_COMM_NULL){
        for (i = 0, j=zz->Obj_Weight_Dim; i < phg->nVtx; i++, j+=cnt){
          tmp[j] = phg->vindex[i+1] - phg->vindex[i]; /* local degree */
        }

        MPI_Allreduce(tmp, final, nvwgt, MPI_FLOAT, MPI_SUM, 
                      phg->comm->col_comm);
      }
      else{
        for (i = 0, j=zz->Obj_Weight_Dim; i < phg->nVtx; i++, j+=cnt){
          final[j] = phg->vindex[i+1] - phg->vindex[i];
        }
      }
    }
    else {                                           /* unit weights */

      for (i = 0, j=zz->Obj_Weight_Dim; i < phg->nVtx; i++, j += cnt)
        final[j] = 1.;
    }
  }

  if (zz->Obj_Weight_Dim){

    dim = zz->Obj_Weight_Dim;

    if (phg->comm->nProc_x == 1)  {
      for (i = 0; i < app.nVtx; i++) {
        idx = app.vtx_gno[i];
        for (j = 0; j < dim; j++)
          final[idx * phg->VtxWeightDim + j] = app.vwgt[i * dim + j];
      }
    }
    else {
      
      /* In using "tmp" for recv, assuming nrecv <= length of "tmp" */
      msg_tag++;
      ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) app.vwgt, 
                            dim*sizeof(float), (char *) tmp);
      
      for (i = 0; i < nrecv; i++) {
        idx = VTX_GNO_TO_LNO(phg, recv_gno[i]);
        for (j = 0; j < dim; j++) 
          final[idx * phg->VtxWeightDim + j] = tmp[i * dim + j];
      }
    }

    /* Reduce weights for all vertices within column 
     * to all processors within column.
     */

    if (phg->comm->col_comm != MPI_COMM_NULL){
      MPI_Allreduce(final, tmp, nvwgt, MPI_FLOAT, MPI_MAX, 
                    phg->comm->col_comm);
    }
  }

  if (!zz->LB.Remap_Flag && zz->LB.Return_Lists == ZOLTAN_LB_NO_LISTS) {
    int gnremove;
    MPI_Allreduce(&(zhg->nRemove), &gnremove, 1, MPI_INT, MPI_SUM, 
                  zz->Communicator);
    if (!hgp->final_output || !gnremove) {
      /* Don't need the plan long-term; destroy it now. */
      Zoltan_Comm_Destroy(&(zhg->VtxPlan));
      if (zhg->Recv_GNOs == app.vtx_gno) app.vtx_gno = NULL;
      ZOLTAN_FREE(&(zhg->Recv_GNOs));
      zhg->nRecv_GNOs = 0;
    }
  }

  /*  Send edge weights, if any */

  if (zz->Edge_Weight_Dim) {
    dim = phg->EdgeWeightDim = zz->Edge_Weight_Dim;
    for (i = 0; i < phg->nEdge; i++) tmpwgts[i] = 0;
    nwgt = phg->nEdge * dim;
    phg->ewgt = (float *) ZOLTAN_CALLOC(nwgt, sizeof(float));
    if (nwgt && !phg->ewgt)
      MEMORY_ERROR;

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
                                phg->comm->col_comm, msg_tag, &nrecv);

      recv_gno = (int *) ZOLTAN_MALLOC(nrecv * sizeof(int));
      if (nrecv && !recv_gno) MEMORY_ERROR;

      msg_tag++;
      ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) app.edge_gno, sizeof(int),
                            (char *) recv_gno);

      /* In using phg->vwgt for recv, assuming nrecv <= phg->nVtx */
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
      MPI_Allreduce(tmpwgts, phg->ewgt, nwgt, MPI_FLOAT, MPI_MAX, 
                    phg->comm->row_comm);
  }
  else {
    /* KDDKDD  For now, do not assume uniform edge weights.
     * KDDKDD  Can add later if, e.g., we decide to coalesce identical edges.
     */
    phg->EdgeWeightDim = 0;
  }

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    Zoltan_HG_HGraph_Free(phg);
  }
  
  Zoltan_Multifree(__FILE__, __LINE__, 11, &app.egids,
                                           &app.elids,
                                           &app.pins, 
                                           &app.esizes, 
                                           &app.pin_procs, 
                                           &app.pin_gno, 
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
 * Pins for removed hyperedges were not retrieved before.
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

  /* Get pin information for removed hyperedges. */

  for (i = 0; i < zhg->nRemove; i++) npins += zhg->Remove_Esize[i];

  if (zhg->nRemove && npins != 0) {
    if (zz->Get_HG_Edge_List) {
      /* Hypergraph callbacks used to build hypergraph */
      
      pins = ZOLTAN_MALLOC_GID_ARRAY(zz, npins);
      pin_procs = (int *) ZOLTAN_MALLOC(npins * sizeof(int));
      if (!pins || !pin_procs) MEMORY_ERROR;

      ierr = zz->Get_HG_Edge_List(zz->Get_HG_Edge_List_Data,
                                  num_gid_entries, num_lid_entries, 
                                  zhg->nRemove,
                                  zhg->Remove_EGIDs, zhg->Remove_ELIDs,
                                  zhg->Remove_Esize, pins, pin_procs);
   
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                           "Error returned from getting HG_Edge_Lists");
        goto End;
      }
    }
    else { 
      /* Graph callbacks used to build hypergraph */
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

  ZOLTAN_FREE(&pins);

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
      j = hash_lookup(zz, &(recvpins[i*num_gid_entries]), nObj, hash_tab);
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

  ZOLTAN_FREE(&pin_procs);
  ZOLTAN_FREE(&parts);

  return ierr;
}
/*****************************************************************************/

static int hash_lookup(
  ZZ *zz,
  ZOLTAN_ID_PTR key,
  int nVtx,
  struct Hash_Node **hash_tab
)
{
/* Looks up a key GID in the hash table; returns its gno. */
/* Based on hash_lookup in build_graph.c. */

  int i;
  struct Hash_Node *ptr;

  i = Zoltan_Hash(key, zz->Num_GID, (unsigned int) nVtx);
  for (ptr=hash_tab[i]; ptr != NULL; ptr = ptr->next){
    if (ZOLTAN_EQ_GID(zz, ptr->gid, key))
      return (ptr->gno);
  }
  /* Key not in hash table */
  return -1;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
