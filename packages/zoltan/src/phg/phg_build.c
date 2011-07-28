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

#include <float.h>
#include "phg.h"
#include "phg_verbose.h"
#include "zz_const.h"
#include "third_library_const.h"
#include "third_library_tools.h"
#include "zz_util_const.h"

/*#define REPART_FASTER_METHOD*/

/*#define DEBUG_FILL_HYPERGRAPH 1*/
  
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

static int Zoltan_PHG_Add_Repart_Data(ZZ *, ZHG *, HGraph *,
                                      PHGPartParams *, Partition);
static int calculate_cuts(ZZ *zz, ZHG *zhg, 
                int max_parts, int *pin_parts, double *loccuts);

static int getObjectSizes(ZZ *zz, ZHG *zhg);

static int remove_dense_edges(ZZ *zz, ZHG *zhg, float esize_threshold, int save_removed,
   int *nEdge, ZOLTAN_GNO_TYPE *nGlobalEdges, int *nPins, ZOLTAN_GNO_TYPE **edgeGNO, int **edgeSize, float **edgeWeight,
   ZOLTAN_GNO_TYPE **pinGNO, int **pinProcs) ;

/*****************************************************************************/
int Zoltan_PHG_Build_Hypergraph(
  ZZ *zz,                     /* Input : Zoltan data structure */
  ZHG **zoltan_hg,            /* Output: Hypergraph to be allocated and built.*/
  Partition *input_parts,     /* Output: Initial partition assignments for
                                 vtxs (in 2D distribution); length = 
                                 zoltan_hg->HG->nVtx.  */
  PHGPartParams *hgp          /* Input : Parameters for PHG partitioning.*/
)
{
/* allocates and builds hypergraph data structure using callback routines */ 
ZHG *zhg;                     /* Temporary pointer to Zoltan_HGraph. */
HGraph *phgraph;             /* Temporary pointer to HG field */
int ierr = ZOLTAN_OK;
char *yo = "Zoltan_PHG_Build_Hypergraph";
char *input_object = "hypergraph";
char msg[128];

  ZOLTAN_TRACE_ENTER(zz, yo);

  /**************************************************************
   * Get the hypergraph specified by the queries and parameters 
  ***************************************************************/

  zhg = (ZHG*) ZOLTAN_MALLOC (sizeof(ZHG));
  if (zhg == NULL) MEMORY_ERROR;

  Zoltan_Input_HG_Init(zhg);

  if (zz->LB.Method == GRAPH){
    input_object = "graph";
  }

  ierr = Zoltan_Get_Hypergraph_From_Queries(zz, hgp, zz->LB.Method, zhg);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg,"Error getting %s from application",input_object);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    goto End;
  }

  *zoltan_hg = zhg;
  phgraph = &(zhg->HG);
  Zoltan_HG_HGraph_Init(phgraph);
  phgraph->comm = &hgp->globalcomm;
  
  /* Obtain the coordinates with 2D distribution. */
  if (zz->Get_Num_Geom != NULL &&
      (zz->Get_Geom != NULL || zz->Get_Geom_Multi != NULL)) {
     /* Geometric callbacks are registered;       */
     /* get coordinates for hypergraph objects.   */

    ZOLTAN_TRACE_DETAIL(zz, yo, "Getting Coordinates.");

     ierr = Zoltan_Get_Coordinates(zz, zhg->nObj, zhg->objGID,
      zhg->objLID, &(phgraph->nDim), &(zhg->coor));
  }

  /**************************************************************
   * Build the hypergraph for PHG from zhg.
  ***************************************************************/

  ierr = Zoltan_PHG_Fill_Hypergraph(zz, zhg, hgp, input_parts);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg,"Error building hypergraph from data supplied by your %s queries",input_object);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    goto End;
  }

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

int Zoltan_PHG_Fill_Hypergraph(
  ZZ *zz,        /* Input : Zoltan data structure */
  ZHG *zhg,      /* Input: Description of hypergraph provided by the application */
                 /* Output: zhg->HG, hypergraph for Zoltan_PHG                   */
  PHGPartParams *hgp,      /* Input : Parameters for PHG partitioning.*/
  Partition *input_parts   /* Output: Initial partition assignment of vtxs in 
                              2D data distribution; length = zhg->HG->nVtx. */
)
{
/* 
 * Create a fully functioning parallel hypergraph with 2D distribution of pins (non-zeros).
 */

char *yo = "Zoltan_PHG_Fill_Hypergraph";

ZOLTAN_COMM_OBJ *plan=NULL;

int ierr = ZOLTAN_OK;
int i, j, cnt, dim, rc;
int msg_tag = 30000;
float *gid_weights = NULL;
char *have_wgt = NULL;
ZOLTAN_GNO_TYPE edge_gno, vtx_gno;
int edge_Proc_y, vtx_Proc_x;
int nnz, idx, method_repart;
int *proclist = NULL;
ZOLTAN_GNO_TYPE *sendbuf = NULL;
int *tmparray = NULL;
int *hindex = NULL;
ZOLTAN_GNO_TYPE *dist_x = NULL, *dist_y = NULL;
ZOLTAN_GNO_TYPE *nonzeros = NULL; 
int *hvertex=NULL;
int nEdge, nVtx, nwgt = 0;
int nrecv = 0; 
ZOLTAN_GNO_TYPE *recv_gno = NULL; 
int *tmpparts = NULL;
int num_gid_entries = zz->Num_GID;
ZOLTAN_GNO_TYPE totalNumEdges, nGlobalEdges;
int nLocalEdges, nPins;
int *edgeSize = NULL, *pinProcs = NULL;
ZOLTAN_GNO_TYPE *edgeGNO = NULL, *pinGNO = NULL;
float *edgeWeight = NULL;
HGraph *phg = &(zhg->HG);
int nProc_x, nProc_y, myProc_x, myProc_y;
int proc_offset;
char msg[250];

float frac_x, frac_y;
float *tmpwgts=NULL; 

float edgeSizeThreshold;
int randomizeInitDist;
int final_output;

int nCoords;
double *tmpcoords = NULL;
 
phg_GID_lookup  *lookup_myObjs = NULL;
int GnFixed=0, nFixed=0;                                     
int *tmpfixed = NULL, *fixedPart = NULL;              
ZOLTAN_ID_PTR fixedGIDs = NULL;
int nRepartEdge = 0, nRepartVtx = 0;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /**************************************************/
  /* Determine parameters                           */
  /**************************************************/

  final_output = hgp->final_output;
  edgeSizeThreshold = hgp->EdgeSizeThreshold;
  method_repart = (!strcasecmp(hgp->hgraph_method, "REPARTITION"));
  randomizeInitDist = hgp->RandomizeInitDist;

  /**************************************************/
  /* Note any fixed vertices                        */
  /**************************************************/
  
  if (zz->Get_Num_Fixed_Obj) {  /* If registered query fixed objects/vertices */
    nFixed = zz->Get_Num_Fixed_Obj (zz->Get_Num_Fixed_Obj_Data, &ierr);
    if (ierr != ZOLTAN_OK) {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Error getting Number of Fixed Objs");
       goto End;
    }
    MPI_Allreduce (&nFixed, &GnFixed, 1, MPI_INT, MPI_SUM, zz->Communicator);
    hgp->UseFixedVtx = GnFixed;  
    
    if (GnFixed && zhg->nObj) {
      zhg->fixed = (int*) ZOLTAN_MALLOC (sizeof(int) * zhg->nObj);
      if (!zhg->fixed) MEMORY_ERROR;
      for (i = 0; i < zhg->nObj; i++)
        zhg->fixed[i] = -1;              /* default - no fixed assignment */
    }
      
    if (GnFixed && nFixed && zhg->nObj)  {

      lookup_myObjs = phg_create_GID_lookup_table(zhg->objGID, zhg->nObj, num_gid_entries);
      if (!lookup_myObjs) MEMORY_ERROR;

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
          j = phg_lookup_GID (lookup_myObjs, fixedGIDs + (i * num_gid_entries));
          if (j < 0)
            FATAL_ERROR ("Unexpected Fixed Vertex GID received");
          zhg->fixed[j] = fixedPart[i];  /* overwrite fixed assignment */
        }
      }
      ZOLTAN_FREE(&fixedGIDs);
      ZOLTAN_FREE(&fixedPart);
      phg_free_GID_lookup_table(&lookup_myObjs);
    }
  }

  /****************************************************************************************
   * If it is desired to remove dense edges, divide the list of edges into
   * two lists.  The ZHG structure will contain the removed edges (if final_output is true), 
   * and the kept edges will be returned.
   ****************************************************************************************/

  totalNumEdges = zhg->globalHedges;

  ierr = remove_dense_edges(zz, zhg, edgeSizeThreshold, final_output, 
            &nLocalEdges, &nGlobalEdges, &nPins,
            &edgeGNO, &edgeSize, &edgeWeight, &pinGNO, &pinProcs);

  if (nGlobalEdges < totalNumEdges){

    /* re-assign edge global numbers if any edges were removed */

    ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, edgeGNO, nLocalEdges, 
                                             randomizeInitDist, &totalNumEdges);

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error reassigning global numbers to edges");
      goto End;
    }
  }

  /****************************************************************************************
   * Compute the distribution of vertices and edges to the 2D data distribution's processor 
   * columns and rows. For now, these distributions are described by arrays dist_x  and dist_y; 
   * in the future, we may prefer a hashing function mapping GIDs to processor columns and rows. 
   ****************************************************************************************/

  nProc_x = phg->comm->nProc_x;
  nProc_y = phg->comm->nProc_y;
  myProc_x = phg->comm->myProc_x;
  myProc_y = phg->comm->myProc_y;

  phg->dist_x = dist_x = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC((nProc_x+1), sizeof(ZOLTAN_GNO_TYPE));
  phg->dist_y = dist_y = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC((nProc_y+1), sizeof(ZOLTAN_GNO_TYPE));

  phg->VtxWeightDim = zhg->objWeightDim;

  if (!dist_x || !dist_y) MEMORY_ERROR;

  frac_x = (float) (zhg->globalObj / (float) nProc_x);
  for (i = 1; i < nProc_x; i++)
    dist_x[i] = (ZOLTAN_GNO_TYPE) (i * frac_x);
  dist_x[nProc_x] = (ZOLTAN_GNO_TYPE)zhg->globalObj;
  
  frac_y = (float)nGlobalEdges / (float)nProc_y;
  for (i = 1; i < nProc_y; i++)
    dist_y[i] = (ZOLTAN_GNO_TYPE) (i * frac_y);
  dist_y[nProc_y] = (ZOLTAN_GNO_TYPE)nGlobalEdges;

  if (Zoltan_overflow_test(dist_x[nProc_x] - dist_x[nProc_x - 1]) ||
      Zoltan_overflow_test(dist_y[nProc_y] - dist_y[nProc_y - 1])) {

    if (zz->Proc == 0){
      ZOLTAN_GNO_TYPE n = (zhg->globalObj > nGlobalEdges) ? zhg->globalObj : nGlobalEdges;
      double sqrtminprocs = ((double)n / INT_MAX) + .5;
      int minprocs = (int)sqrtminprocs * (int)sqrtminprocs;
      sprintf(msg,
       "At least %d processes are required to achieve the 2D distribution of the problem matrix\n",
       minprocs);
    }

    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    goto End;
  }

  /* myProc_y and myProc_x can be -1 when nProc is prime and we use a 2D
   * decomposition.  One processor is excluded from the 2D communicator;
   * for it, myProc_y and myProc_x == -1. */
  nEdge = (myProc_y >= 0 ? (int)(dist_y[myProc_y+1] - dist_y[myProc_y]) : 0);
  nVtx  = (myProc_x >= 0 ? (int)(dist_x[myProc_x+1] - dist_x[myProc_x]) : 0);

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
    nRepartEdge = NumRepart(myProc_y, nProc_y, zhg->globalObj);
  }

  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution. 
   */
  proclist = (int *)ZOLTAN_MALLOC(MAX(nPins, zhg->nObj)*sizeof(int));
  sendbuf = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nPins * 2 * sizeof(ZOLTAN_GNO_TYPE));

  if ( ((nPins || zhg->nObj) && !proclist) || (nPins && !sendbuf)){
    MEMORY_ERROR;
  }

  cnt = 0; 
  for (i = 0; i < nLocalEdges; i++) {
    /* processor row for the edge */
    edge_gno = edgeGNO[i];
    edge_Proc_y = EDGE_TO_PROC_Y(phg, edge_gno);

    for (j = 0; j < edgeSize[i]; j++) {
      /* processor column for the vertex */
      vtx_gno = pinGNO[cnt];
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
  ierr = Zoltan_Comm_Create(&plan, cnt, proclist, zz->Communicator, msg_tag, &nnz);

  if (nnz) {
    nonzeros = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nnz * 2 * sizeof(ZOLTAN_GNO_TYPE));
    if (!nonzeros) MEMORY_ERROR;
  }

  msg_tag--;
  Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(ZOLTAN_GNO_TYPE), (char *) nonzeros);

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
    for (i = 0; i < zhg->nObj; i++)
      proclist[i] = proc_offset + VTX_TO_PROC_X(phg, zhg->objGNO[i]);
      
    msg_tag++;
    ierr = Zoltan_Comm_Create(&(zhg->VtxPlan), zhg->nObj, proclist, 
                              zz->Communicator, msg_tag, &nrecv);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
    zhg->nRecv_GNOs = nrecv;

    zhg->Recv_GNOs = recv_gno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nrecv * sizeof(ZOLTAN_GNO_TYPE));

    if (nrecv && !recv_gno) MEMORY_ERROR;

    /* Use plan to send global numbers to the appropriate proc_x. */
    msg_tag++;
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) zhg->objGNO, sizeof(ZOLTAN_GNO_TYPE), (char *) recv_gno);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  }
  else {
    /* Save map of what needed. */
    zhg->nRecv_GNOs = nrecv = zhg->nObj;
    zhg->Recv_GNOs = recv_gno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nrecv * sizeof(ZOLTAN_GNO_TYPE));

    if (nrecv && !recv_gno) MEMORY_ERROR;
   
    memcpy(zhg->Recv_GNOs, zhg->objGNO, nrecv * sizeof(ZOLTAN_GNO_TYPE));
  }

  /* Send vertex partition assignments and weights to 2D distribution. */

  tmpparts = (int *) ZOLTAN_CALLOC(phg->nVtx, sizeof(int));
  if (phg->nVtx && !tmpparts) MEMORY_ERROR;

  *input_parts = (int *) ZOLTAN_MALLOC((phg->nVtx + nRepartVtx) * sizeof(int));
  if ((phg->nVtx || nRepartVtx) && !*input_parts) MEMORY_ERROR;

  dim = phg->VtxWeightDim;
  nwgt = (phg->nVtx + nRepartVtx) * dim;

  tmpwgts = (float *)ZOLTAN_CALLOC(sizeof(float), nwgt);
  phg->vwgt = (float *)ZOLTAN_MALLOC(sizeof(float) * nwgt);
  if (nwgt && (!tmpwgts || !phg->vwgt)) MEMORY_ERROR;

  if (phg->nDim > 0) {
    nCoords = (phg->nVtx + nRepartVtx) * phg->nDim;

    tmpcoords = (double *)ZOLTAN_CALLOC(sizeof(double), nCoords);
    phg->coor = (double *)ZOLTAN_MALLOC(sizeof(double) * nCoords);
    if (nCoords && (!tmpcoords || !phg->coor)) MEMORY_ERROR;
  }
  
  if (GnFixed) {                              
    tmpfixed   = (int*) ZOLTAN_MALLOC(phg->nVtx * sizeof(int));
    if (phg->nVtx && !tmpfixed) MEMORY_ERROR;

    phg->fixed_part = (int*) ZOLTAN_MALLOC((phg->nVtx + nRepartVtx) * sizeof(int));
    if ((phg->nVtx || nRepartVtx) && !phg->fixed_part) MEMORY_ERROR;

    for (i = 0 ; i < phg->nVtx; i++)
      tmpfixed[i] = -2;  
  } 

  if (phg->comm->nProc_x == 1)  {
    for (i = 0; i < zhg->nObj; i++) {
      idx = zhg->objGNO[i];
      tmpparts[idx] = zhg->Input_Parts[i];
      for (j = 0; j < dim; j++)
        tmpwgts[idx*dim + j] = zhg->objWeight[i*dim + j];
      if (phg->nDim > 0)
	for (j = 0; j < phg->nDim; j++){
	  tmpcoords[idx*phg->nDim + j] = zhg->coor[i*phg->nDim + j];
	  phg->coor[i*phg->nDim + j] = -DBL_MAX; /* defined in float.h */
	}
      if (GnFixed)                                             
        tmpfixed[idx] = zhg->fixed[i];
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
    ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) zhg->objWeight,
                          sizeof(float) * dim, (char *) phg->vwgt);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    if (phg->nDim > 0) {
      msg_tag++;
      ierr = Zoltan_Comm_Do(zhg->VtxPlan, msg_tag, (char *) zhg->coor,
			    sizeof(double) * phg->nDim, (char *) phg->coor);

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
	goto End;
      }
    }
    
    if (GnFixed) {                                  
       msg_tag++;
       ierr = Zoltan_Comm_Do (zhg->VtxPlan, msg_tag, (char*) zhg->fixed,
         sizeof(int), (char*) phg->fixed_part);
       if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN))
         goto End;         
    }
       
    for (i = 0; i < nrecv; i++) {
      idx = VTX_GNO_TO_LNO(phg, recv_gno[i]);
      tmpparts[idx] = (*input_parts)[i];
      for (j=0; j<dim; j++)
        tmpwgts[idx*dim + j] = phg->vwgt[i*dim + j];  
      if (GnFixed)                                       
	tmpfixed[idx] = phg->fixed_part[i];         
      if (phg->nDim > 0) {
	for (j = 0; j < phg->nDim; j++){
	  tmpcoords[idx*phg->nDim + j] = phg->coor[i*phg->nDim + j];
	  phg->coor[i*phg->nDim + j] = -DBL_MAX;
	}
      }
    }
  }

  /* Reduce partition assignments, weights and coordinates for all
   * vertices within column to all processors within column.
   */

  if (phg->comm->col_comm != MPI_COMM_NULL){
    rc = MPI_Allreduce(tmpparts, *input_parts, phg->nVtx, MPI_INT, MPI_MAX, 
                  phg->comm->col_comm);
    CHECK_FOR_MPI_ERROR(rc);

    rc = MPI_Allreduce(tmpwgts, phg->vwgt, nwgt, MPI_FLOAT, MPI_MAX,
                  phg->comm->col_comm);
    CHECK_FOR_MPI_ERROR(rc);

    if (phg->nDim > 0) {
      rc = MPI_Allreduce(tmpcoords, phg->coor, nCoords, MPI_DOUBLE, MPI_SUM,
			 phg->comm->col_comm);
      CHECK_FOR_MPI_ERROR(rc);
    }
    
    if (GnFixed) {                                          
       rc = MPI_Allreduce(tmpfixed, phg->fixed_part, phg->nVtx, MPI_INT, MPI_MAX,
            phg->comm->col_comm);
       CHECK_FOR_MPI_ERROR(rc);
    }
  }
  
  ZOLTAN_FREE(&tmpcoords);
  ZOLTAN_FREE(&tmpfixed);
  ZOLTAN_FREE(&tmpparts);
  ZOLTAN_FREE(&tmpwgts);
  
  /*  Send edge weights, if any */

  dim = zhg->edgeWeightDim;
  if (method_repart && (!dim))
    dim = 1; /* Need edge weights for REPARTITION; force malloc of ewgt array */
  phg->EdgeWeightDim = dim;
  nwgt = (phg->nEdge + nRepartEdge) * dim;

  if (dim) {
    tmpwgts   = (float *) ZOLTAN_CALLOC(nwgt, sizeof(float));
    phg->ewgt = (float *) ZOLTAN_CALLOC(nwgt, sizeof(float));

    if (nwgt && (!phg->ewgt || !tmpwgts)) MEMORY_ERROR;

    if (zhg->edgeWeightDim) { 
      /* 
       * Edge weights provided by application or set in Zoltan_Get_Hypergraph_From_Queries
       */
      if (phg->comm->nProc_y == 1) {
        for (i = 0; i < nLocalEdges; i++) {
          idx = edgeGNO[i];
          for (j = 0; j < dim; j++)
            tmpwgts[idx * dim + j] = edgeWeight[i*dim + j];
        }
      }
      else {
        /* 
         * Since for 2D decomposition and prime nProc we exclude a processor,
         * we cannot use col_comm for this operation.  We'll simulate it,
         * allowing the excluded processor to send its data to col 0.
         */
        proc_offset = (myProc_x >= 0 ? myProc_x : 0);
        for (i = 0; i < nLocalEdges; i++)
          proclist[i] = proc_offset 
                      + EDGE_TO_PROC_Y(phg, edgeGNO[i]) * nProc_x;
        
        msg_tag++;
  
        ierr = Zoltan_Comm_Create(&plan, nLocalEdges, proclist, 
                                  zz->Communicator, msg_tag, &nrecv); 
  
        if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
          goto End;
        }
  
        /* Multiple processes may have weights for the same edge */
  
        recv_gno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nrecv * sizeof(ZOLTAN_GNO_TYPE));
  
        gid_weights = (float *) ZOLTAN_CALLOC(nrecv*dim, sizeof(float));
  
        have_wgt = (char *)ZOLTAN_CALLOC(phg->nEdge, sizeof(char));
  
        if ((nrecv && (!recv_gno || !gid_weights)) ||
            (phg->nEdge && !have_wgt)){
  
          ZOLTAN_FREE(&recv_gno);
          MEMORY_ERROR;
        }
  
        msg_tag++;
        ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) edgeGNO, sizeof(ZOLTAN_GNO_TYPE), (char *) recv_gno);
  
        if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
          goto End;
        }
  
        msg_tag++;
        ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) edgeWeight, 
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
        CHECK_FOR_MPI_ERROR(rc);
      }
    }
    else { /* dim > 0 but zhg->edgeWeightDim == 0 */
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

/* HERE phg->ewgt are all 3.0 */


  if (method_repart){
    ierr = Zoltan_PHG_Add_Repart_Data(zz, zhg, phg, hgp, *input_parts);
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
    rc = MPI_Allreduce(&(zhg->nHedges), &gnremove, 1, MPI_INT, MPI_SUM, 
                  zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc);
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

  ZOLTAN_FREE(&edgeSize);
  ZOLTAN_FREE(&edgeGNO);
  ZOLTAN_FREE(&edgeWeight);
  ZOLTAN_FREE(&pinGNO);
  ZOLTAN_FREE(&pinProcs);

  Zoltan_Comm_Destroy(&plan);

  Zoltan_Multifree(__FILE__, __LINE__, 9, 
    &proclist,
    &sendbuf,
    &nonzeros,
    &tmparray,
    &tmpparts,
    &tmpwgts,
    &tmpfixed,
    &tmpcoords,
    &gid_weights,
    &have_wgt);

#ifdef DEBUG_FILL_HYPERGRAPH
  for (i=0; i<zz->Num_Proc; i++){
    if (i == zz->Proc){
      printf("HYPERGRAPH on process %d\n",i);

      if (hgp->add_obj_weight == PHG_ADD_PINS_WEIGHT) 
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

int Zoltan_PHG_Cuts(
  ZZ *zz,
  ZHG *zhg,
  double *localcuts /* Array of length 2. Upon return:
                     * localcuts[0] = ConCut: Sum_over_edges( (nparts-1)*ewgt )
                     * localcuts[1] = NetCut: Sum_over_edges( (nparts>1)*ewgt )
                     */
)
{

/* Function to compute the cuts of hyperedges listed in zhg.  */
/* Edges in zhg are entire edges, they are not distributed across processes */

static char *yo = "Zoltan_PHG_Cuts";
int ierr = ZOLTAN_OK;
int i;
ZOLTAN_MAP *map;
int npins = 0;                   /* # of pins in hyperedges */
ZOLTAN_GNO_TYPE *pins = NULL;                   /* pins for edges */
int *pin_procs = NULL;           /* procs owning pins for edges */
int *pin_parts = NULL;           /* parts computed for pins in edges */
ZOLTAN_GNO_TYPE *myObjGNO = NULL;
ZOLTAN_COMM_OBJ *plan = NULL;
int nrecv;                       /* # of requests for pin info */
ZOLTAN_GNO_TYPE *recvpins = NULL;            /* Requested pins from other procs */
int *outparts = NULL;            /* received partition info for pins */
int num_parts, max_parts;
int msg_tag = 23132;
int index;
intptr_t iptr;                   /* an int the size of a pointer */

  ZOLTAN_TRACE_ENTER(zz, yo);

  npins = zhg->nPins;
  pins = zhg->pinGNO;
  pin_procs = zhg->Pin_Procs;
  myObjGNO = zhg->objGNO;

  /* Communicate Output_Part info for off-processor pins */

  Zoltan_Comm_Create(&plan, npins, pin_procs, zz->Communicator, msg_tag, &nrecv);

  if (nrecv) {
    recvpins = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * nrecv);
    outparts = (int *)ZOLTAN_MALLOC(sizeof(int) * nrecv);
    if (!recvpins || !outparts) MEMORY_ERROR;
  }

  Zoltan_Comm_Do(plan, msg_tag, (char *) pins, sizeof(ZOLTAN_GNO_TYPE), (char *) recvpins);

  if (nrecv) {

    map = Zoltan_Map_Create(zz, 0, sizeof(ZOLTAN_GNO_TYPE), 0, zhg->nObj);
    if (map == NULL) goto End;

    for (iptr=0; iptr < zhg->nObj; iptr++){
      ierr = Zoltan_Map_Add(zz, map, (char *)(myObjGNO + iptr), iptr + 1);
      if (ierr != ZOLTAN_OK) goto End;
    }
    
    for (i = 0; i < nrecv; i++) {
      ierr = Zoltan_Map_Find(zz, map, (char *)(recvpins + i), &iptr);
      if (ierr != ZOLTAN_OK) goto End;
      if ((void *)iptr == NULL){
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in pin map.");
         goto End;
      }

      index = (int)iptr - 1;
      outparts[i] = zhg->Output_Parts[index];
    }
    Zoltan_Map_Destroy(zz, &map);
  }

  ZOLTAN_FREE(&recvpins);

  /* Send partition info back to requesting processor */
  pin_parts = (int *) ZOLTAN_MALLOC(npins * sizeof(int));
  if (npins && !pin_parts) MEMORY_ERROR;

  Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *) outparts, sizeof(int), NULL, (char *) pin_parts);

  ZOLTAN_FREE(&outparts);

  Zoltan_Comm_Destroy(&plan);

  /* Compute the cut metrics using received partition info.
   *
   * Get maximum possible partition number.  (It's not always the
   * value found in the ZZ structure.)
   */

  num_parts = 0;
  for (i=0; i < npins; i++){
    if (pin_parts[i] >= num_parts) num_parts = pin_parts[i]+1;
  }

  MPI_Allreduce(&num_parts, &max_parts, 1, MPI_INT, MPI_MAX, zz->Communicator);

  /*
   * Calculate the cut metrics.  We assume that edges are not divided
   * across processes, so each process does this calculation locally.
   */

  ierr = calculate_cuts(zz, zhg, max_parts, pin_parts, localcuts);

End:

  ZOLTAN_FREE(&pin_parts);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}

/*****************************************************************************/

static int calculate_cuts(ZZ *zz, ZHG *zhg, 
                int max_parts, int *pin_parts, double *loccuts)
{
char *yo = "calculate_cuts";
int i, cnt, j, ierr, nparts;
int *parts;
float ewgt;

  ierr = ZOLTAN_OK;
  parts = (int *) ZOLTAN_CALLOC(max_parts, sizeof(int));
  if (max_parts && !parts) MEMORY_ERROR;

  cnt = 0;
  loccuts[0] = loccuts[1] = 0.;
  for (i = 0; i < zhg->nHedges; i++) {
    nparts = 0;
    for (j = 0; j < zhg->Esize[i]; j++) {
      if (parts[pin_parts[cnt]] < i+1) {
        nparts++;
      }
      parts[pin_parts[cnt]] = i+1;
      cnt++;
    }
    ewgt = (zhg->Ewgt ? zhg->Ewgt[i*zhg->edgeWeightDim] : 1.);
    if (nparts > 1) {
      loccuts[0] += (nparts-1) * ewgt;
      loccuts[1] += ewgt;
    }
  }

End:
  ZOLTAN_FREE(&parts);
  return ierr;
}

/****************************************************************************/

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
                               zhg->objGID, zhg->objLID, zhg->AppObjSizes, &ierr);
        if (ierr < 0) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                          "ZOLTAN_OBJ_SIZE_MULTI function.");
          goto End;
        }
      }
      else if (zz->Get_Obj_Size) {
        for (i = 0; i < zhg->nObj; i++) {
          ZOLTAN_ID_PTR lid = (zz->Num_LID ? &(zhg->objLID[i*zz->Num_LID]):NULL);
          zhg->AppObjSizes[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data,
                                           zz->Num_GID, zz->Num_LID,
                                           &(zhg->objGID[i*zz->Num_GID]),
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
ZOLTAN_GNO_TYPE maxgnRepartEdge = zhg->globalObj; /* max # of repartition edges that could
                                                     be added to the global hypergraph. */
int myProc_x = hgc->myProc_x;
int myProc_y = hgc->myProc_y;
int nProc_x = hgc->nProc_x;
int nProc_y = hgc->nProc_y;
int nRepartVtx = 0, nRepartEdge;    /* # of repartition vertices & edges 
                                        to add to this processor. */
int nRepartPin = 0;                  /* # of repartition pins to add to this
                                        processor. */
ZOLTAN_GNO_TYPE firstRepartEdge=0;
int firstRepartVtx = 0;              /* global no. (0-based) of the first repartition
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
ZOLTAN_GNO_TYPE *sendgno = NULL;     /* to procs owning the corresponding   */
ZOLTAN_GNO_TYPE *recvgno = NULL;     /* repartition vertices and edges.     */


#ifndef REPART_FASTER_METHOD
int *sendpart = NULL;
int *recvpart = NULL;
int *sendsize = NULL;
int *recvsize = NULL;
#endif

int nrecv = 0;
                                 
int *tmp_hindex;                     /* Pointers into phg->hindex and 
                                        phg->hvertex; set to start of entries
                                        for repartition edges. */
int *tmp_hvertex;
float *tmp_ewgt = NULL;              /* Edge weights for local repartion
                                        edges */
int *pins_per_edge = NULL;           /* # of pins in each repartition edge. */
int i, j, idx;
int cnt, tmp, dim, sum, prev;
int ierr = ZOLTAN_OK;
int *objsize = NULL;                 /* repartition edge weights. */
int *tmpobjsize = NULL;
int *colProc_cnt = NULL;             /* actual nRepartEdge per column proc */
ZOLTAN_GNO_TYPE *repart_dist_x = NULL; /* Distribution of repartition vertices
                                        to proc rows; similar to phg->dist_x. */

ZOLTAN_GNO_TYPE *repart_dist_y = NULL; /* Distribution of repartition edges
                                        to proc cols; similar to phg->dist_y. */

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (myProc_x >= 0 && myProc_y >= 0) {  /* This proc is part of 2D decomp */

    /* Compute distribution of repartition vertices and edges. repart_dist_x is a
     * ZOLTAN_GNO_TYPE only because ProcForRepart needs a ZOLTAN_GNO_TYPE array */

    repart_dist_x = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC((2+nProc_x+nProc_y) * sizeof(ZOLTAN_GNO_TYPE));

    repart_dist_y = repart_dist_x + 1+ nProc_x;

    if (!repart_dist_x || !repart_dist_y) MEMORY_ERROR;

    for (i = 0; i < nProc_x; i++)
      repart_dist_x[i] = FirstRepart(i, nProc_x, gnRepartVtx);
    repart_dist_x[nProc_x] = gnRepartVtx;

    for (i = 0; i < nProc_y; i++)
      repart_dist_y[i] = FirstRepart(i, nProc_y, maxgnRepartEdge);
    repart_dist_y[nProc_y] = maxgnRepartEdge;

    /* Compute number of repartition vertices to add to this processor column */
    phg->nRepartVtx = nRepartVtx
                    = (int)(repart_dist_x[myProc_x+1] - repart_dist_x[myProc_x]);
    firstRepartVtx = (int)repart_dist_x[myProc_x];

    /* Compute maximum number of repartition edges to add to this proc row */
    phg->nRepartEdge = nRepartEdge 
                     = (int)(repart_dist_y[myProc_y+1] - repart_dist_y[myProc_y]);
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
    cnt = firstRepartVtx;     /* Also the part that the first repart Vertex is associated with */
    dim = phg->VtxWeightDim;  /* since repartition vertex global numbers begin at 0 */
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
        tmpobjsize[zhg->objGNO[i]] = zhg->AppObjSizes[i];
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

    if (cnt) {
      proclist = (int *) ZOLTAN_MALLOC(2 * cnt * sizeof(int));
      if (!proclist) MEMORY_ERROR;

#ifdef REPART_FASTER_METHOD
      /* Most likely faster, but uses more memory */
#define NSEND 3  /* number of fields to be sent for each vertex */
      sendgno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(NSEND * 2 * cnt * sizeof(ZOLTAN_GNO_TYPE));
      if (!sendgno) MEMORY_ERROR;
#else
      /* uses less memory */
      sendgno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(2 * cnt * sizeof(ZOLTAN_GNO_TYPE));
      sendpart = (int *) ZOLTAN_MALLOC(2 * cnt * sizeof(int));
      if (!sendgno || !sendpart) MEMORY_ERROR;

      if (zz->Get_Obj_Size_Multi || zz->Get_Obj_Size) {
        sendsize = (int *) ZOLTAN_MALLOC(2 * cnt * sizeof(int));
        if (!sendsize) MEMORY_ERROR;
      }
#endif
  
      for (i = 0; i < cnt; i++){
        int vtxproc_x;            /* Processor column for actual vtx. */
        int proc_x;               /* Processor column for repartition vtx. */
        int proc_y;               /* Processor row for repartition edge. */
        int v = i + myStart_vtx;  /* Index of vertex being processed. */
        ZOLTAN_GNO_TYPE vgno;                 /* Gno for vertex v. */
        ZOLTAN_GNO_TYPE k = (ZOLTAN_GNO_TYPE)input_parts[v];   /* Repartition vtx global number */

        if (k < gnRepartVtx) {
          vgno = VTX_LNO_TO_GNO(phg,v);  /* Gno of vertex */

          /* Determine proc col for repartition vtx associated with part k */
          proc_x = ProcForRepart(k, repart_dist_x, nProc_x);

          /* Determine proc row for repartition edge associated with vtx i */
          proc_y = ProcForRepart(vgno, repart_dist_y, nProc_y);
 
          /* Fill message buffers */
  
          /* Send message to processor owning repartition vertex and edge */
          proclist[tmp] = proc_y * nProc_x + proc_x;

#ifdef REPART_FASTER_METHOD
          sendgno[NSEND*tmp] = vgno;            /* GNO of vertex */
          sendgno[NSEND*tmp+1] = k;             /* Vertex's input partition */
          sendgno[NSEND*tmp+2] = (objsize ? objsize[v] : 1);
#else
          sendgno[tmp] = vgno;
          sendpart[tmp] = (int)k;
          if (zz->Get_Obj_Size_Multi || zz->Get_Obj_Size) {
            sendsize[tmp] = (objsize ? objsize[v] : 1);
          }
#endif
          tmp++;

          /* Send message to proc owning actual vertex and repartition edge */
          vtxproc_x = VTX_TO_PROC_X(phg,vgno);
          if (vtxproc_x != proc_x) {
            proclist[tmp] = proc_y * nProc_x + vtxproc_x;

#ifdef REPART_FASTER_METHOD
            sendgno[NSEND*tmp] = vgno;            /* GNO of vertex */
            sendgno[NSEND*tmp+1] = k;             /* Vertex's input partition */
            sendgno[NSEND*tmp+2] = (objsize ? objsize[v] : 1);
#else
            sendgno[tmp] = vgno;
            sendpart[tmp] = (int)k;
            if (zz->Get_Obj_Size_Multi || zz->Get_Obj_Size) {
              sendsize[tmp] = (objsize ? objsize[v] : 1);
            }
#endif
            tmp++;
          }
        }
      }
    }

    ZOLTAN_FREE(&objsize);
    ZOLTAN_FREE(&repart_dist_x);

    /* Create plan and send info. */

    ierr = Zoltan_Comm_Create(&plan, tmp, proclist, hgc->Communicator, 24555, &nrecv);

#ifdef REPART_FASTER_METHOD
    if (nrecv && !(recvgno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(NSEND*nrecv*sizeof(ZOLTAN_GNO_TYPE))))
      MEMORY_ERROR;
  
    ierr = Zoltan_Comm_Do(plan, 24556, (char *) sendgno, NSEND * sizeof(ZOLTAN_GNO_TYPE),
                         (char *) recvgno);

    ZOLTAN_FREE(&sendgno);
#else

    if (nrecv && !(recvgno = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nrecv*sizeof(ZOLTAN_GNO_TYPE))))
      MEMORY_ERROR;
  
    ierr = Zoltan_Comm_Do(plan, 24556, (char *) sendgno, sizeof(ZOLTAN_GNO_TYPE), (char *) recvgno);

    ZOLTAN_FREE(&sendgno);

    if (nrecv && !(recvpart = (int *) ZOLTAN_MALLOC(nrecv*sizeof(int))))
      MEMORY_ERROR;
  
    ierr = Zoltan_Comm_Do(plan, 24557, (char *) sendpart, sizeof(int), (char *) recvpart);

    ZOLTAN_FREE(&sendpart);

    if (zz->Get_Obj_Size_Multi || zz->Get_Obj_Size) {
      if (nrecv && !(recvsize = (int *) ZOLTAN_MALLOC(nrecv*sizeof(int))))
        MEMORY_ERROR;

      ierr = Zoltan_Comm_Do(plan, 24558, (char *) sendsize, sizeof(int), (char *) recvsize);

      ZOLTAN_FREE(&sendsize);
    }
#endif

    ierr = Zoltan_Comm_Destroy(&plan);

    /* Build hindex & hvertex for repartition vertices & repartition edges.
     * Both arrays are appended to existing arrays, which are assumed to be
     * large enough.
     */

    tmp_hindex = &(phg->hindex[phg->nEdge]);
    tmp_hvertex = phg->hvertex;
  
    /* Loop to build hindex array */
    for (i = 0; i < nrecv; i++) {
      ZOLTAN_GNO_TYPE vtx_gno;   /* Global vtx number of vtx in received repartition edge.*/
      int rEdge_lno; /* local index of repartition edge */
      ZOLTAN_GNO_TYPE rVtx_gno;  /* global repartition vertex number */
      int rVtx_lno;  /* local index of repartition vertex */
  
#ifdef REPART_FASTER_METHOD
      vtx_gno = recvgno[NSEND*i];
      rEdge_lno = (int)(vtx_gno - firstRepartEdge);
      rVtx_gno = recvgno[NSEND*i+1];
      rVtx_lno = (int)(rVtx_gno - firstRepartVtx);
#else
      vtx_gno = recvgno[i];
      rEdge_lno = (int)(vtx_gno - firstRepartEdge);
      rVtx_gno = recvpart[i];
      rVtx_lno = (int)(rVtx_gno - firstRepartVtx);
#endif

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
      ZOLTAN_GNO_TYPE vtx_gno;   /* Global vtx number of vtx in received repartition edge.*/
      int rEdge_lno; /* local index of repartition edge */
      ZOLTAN_GNO_TYPE rVtx_gno;  /* global repartition vertex number */
      int rVtx_lno;  /* local index of repartition vertex */
  
#ifdef REPART_FASTER_METHOD
      vtx_gno = recvgno[NSEND*i];
      rEdge_lno = (int)(vtx_gno - firstRepartEdge);
      rVtx_gno = recvgno[NSEND*i+1];
      rVtx_lno = (int)(rVtx_gno - firstRepartVtx);
#else
      vtx_gno = recvgno[i];
      rEdge_lno = (int)(vtx_gno - firstRepartEdge);
      rVtx_gno = recvpart[i];
      rVtx_lno = (int)(rVtx_gno - firstRepartVtx);
#endif
      
      cnt = 0;
      if (rVtx_gno >= firstRepartVtx && rVtx_gno < firstRepartVtx+nRepartVtx) {
        /* Add pin for repartition vertex */
        tmp_hvertex[tmp_hindex[rEdge_lno]] = phg->nVtx + rVtx_lno;
        cnt++;
#ifdef REPART_FASTER_METHOD
        for (j = 0; j < phg->EdgeWeightDim; j++) 
          tmp_ewgt[rEdge_lno*phg->EdgeWeightDim+j] = (float) recvgno[NSEND*i+2];
#else
        if (zz->Get_Obj_Size_Multi || zz->Get_Obj_Size) {
          for (j = 0; j < phg->EdgeWeightDim; j++) 
            tmp_ewgt[rEdge_lno*phg->EdgeWeightDim+j] = (float) recvsize[i];
        }
        else{
          for (j = 0; j < phg->EdgeWeightDim; j++) 
            tmp_ewgt[rEdge_lno*phg->EdgeWeightDim+j] = 1.0;
        }
#endif
      }
      if (myProc_x == VTX_TO_PROC_X(phg, vtx_gno)) {
        /* Add pin for vertex */
        tmp_hvertex[tmp_hindex[rEdge_lno]+cnt] = VTX_GNO_TO_LNO(phg, vtx_gno);
      }
    }

    ZOLTAN_FREE(&recvgno);

#ifdef REPART_FASTER_METHOD
#undef NSEND
#else
    ZOLTAN_FREE(&recvpart);
    ZOLTAN_FREE(&recvsize);
#endif

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
    if (tmp_ewgt) ZOLTAN_FREE(&tmp_ewgt);

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

    if (!colProc_cnt) MEMORY_ERROR;

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
  ZOLTAN_FREE(&sendgno);
  ZOLTAN_FREE(&recvgno);

#ifndef REPART_FASTER_METHOD
  ZOLTAN_FREE(&sendpart);
  ZOLTAN_FREE(&recvpart);
  ZOLTAN_FREE(&sendsize);
  ZOLTAN_FREE(&recvsize);
#endif

  ZOLTAN_TRACE_EXIT(zz, yo);
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
ZOLTAN_GNO_TYPE sum;
int i;
int *colProc_cnt = NULL;
int myProc_x = phg->comm->myProc_x;
int myProc_y = phg->comm->myProc_y;
int ierr = ZOLTAN_OK;

  if (myProc_x >= 0 && myProc_y >= 0) {  /* This proc is part of 2D decomp */
    /* Only processors in the 2D decomp added repartition data, so 
     * only those processors need to remove it.
     */
    phg->nVtx  -= phg->nRepartVtx;
    phg->nEdge -= phg->nRepartEdge;
    phg->nPins -= phg->nRepartPin;

    if (zhg->edgeWeightDim) {
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
      phg->dist_x[i] -= (int)sum;
      sum += NumRepart(i, nProc_x, zhg->GnRepartVtx);
    }
    phg->dist_x[nProc_x] -= (int)sum;
  
    colProc_cnt = (int *) ZOLTAN_MALLOC(nProc_y * sizeof(int));
    if (!colProc_cnt){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    MPI_Allgather(&(phg->nRepartEdge), 1, MPI_INT, colProc_cnt, 1, MPI_INT, 
                  phg->comm->col_comm);
    for (sum = 0, i = 0; i < nProc_y; i++) {
      phg->dist_y[i] -= sum;
      sum += colProc_cnt[i];
    }
    phg->dist_y[nProc_y] -= sum;
    ZOLTAN_FREE(&colProc_cnt);
  }

End:
  return ierr;
}
/****************************************************************************/

static int remove_dense_edges(ZZ *zz, 
   ZHG *zhg,                            /* on output may contain removed edges */
   float esize_threshold,               /* %age of gnVtx considered a dense edge */
   int save_removed,                    /* save the removed edges in ZHG */
   int *nEdge,                          /* number of edges remaining */
   ZOLTAN_GNO_TYPE *nGlobalEdges,       /* number of global edges remaining */
   int *nPins,                          /* number of local pins remainint */
   ZOLTAN_GNO_TYPE **edgeGNO,           /* global numbers of remaining edges */
   int **edgeSize,                      /* number pins in each remaining edge */
   float **edgeWeight,                  /* weight of remaining edges */
   ZOLTAN_GNO_TYPE **pinGNO,            /* global numbers of pins in remaining edges */
   int **pinProcs)                     /* processes owning pin vertices */
{
char *yo = "remove_dense_edges";
int ierr = ZOLTAN_OK;
int i, w, esize, k, l, kpin, lpin, pin;
int ew_dim = zhg->edgeWeightDim;
ZOLTAN_GNO_TYPE global_nremove = 0; 
ZOLTAN_GNO_TYPE global_nremove_pins = 0; 
ZOLTAN_GNO_TYPE nremove = 0;
ZOLTAN_GNO_TYPE nremove_size = 0;
ZOLTAN_GNO_TYPE nkeep = 0;
ZOLTAN_GNO_TYPE nkeep_size = 0;
ZOLTAN_GNO_TYPE numGlobalEdges = 0;
double gesize_threshold;
int *goEdgeSize = NULL;
int *goPinProc = NULL;
float *goEdgeWeight = NULL;
int *keepEdgeSize = NULL;
int *keepPinProc = NULL;
float *keepEdgeWeight = NULL;
ZOLTAN_GNO_TYPE *keepPinGNO = NULL, *keepEdgeGNO = NULL, *goPinGNO = NULL, *goEdgeGNO = NULL, *gnoPtr;
int *procPtr;
float *wgtPtr;
ZOLTAN_GNO_TYPE localval[2], globalval[2];
MPI_Datatype zoltan_gno_mpi_type;

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  /* Remove dense edges and zero-sized edges from input list */

  if (esize_threshold <= 1)
    gesize_threshold = esize_threshold * (double)zhg->globalObj; /* relative */
  else
    gesize_threshold = (double)esize_threshold; /* absolute */

  for (i = 0; i < zhg->nHedges; i++) {
    if ((zhg->Esize[i] > gesize_threshold) || (zhg->Esize[i] == 0)){
      nremove++;
      nremove_size += zhg->Esize[i];
    }
    else{
      nkeep++;
      nkeep_size += zhg->Esize[i];
    }
  }

  localval[0] = nremove;
  localval[1] = nremove_size;

  MPI_Allreduce(localval, globalval, 2, zoltan_gno_mpi_type, MPI_SUM, zz->Communicator);


  global_nremove = globalval[0];
  global_nremove_pins = globalval[1];

  numGlobalEdges = zhg->globalHedges - global_nremove;

  if ((zhg->globalHedges > zz->Num_Proc) && (numGlobalEdges < zz->Num_Proc)){
    if (zz->Proc == 0) fprintf(stderr,
     "\nWARNING: PHG_EDGE_SIZE_THRESHOLD is low (%f), resulting in only " ZOLTAN_GNO_SPEC " edges\n remaining.",
        esize_threshold, numGlobalEdges);
  }

  *nGlobalEdges = numGlobalEdges;
  *nEdge = nkeep;
  *nPins = nkeep_size;

  if (nremove == 0){

    /* kept edges */

    *edgeGNO = zhg->edgeGNO;
    *edgeSize = zhg->Esize;
    *edgeWeight = zhg->Ewgt;
    *pinGNO = zhg->pinGNO;
    *pinProcs = zhg->Pin_Procs;

    /* removed edges */

    zhg->nHedges = 0;
    zhg->globalHedges = 0;
    zhg->nPins = 0;
    zhg->globalPins = 0;
    zhg->edgeGNO = NULL;
    zhg->Esize = NULL;
    zhg->Ewgt = NULL;
    zhg->pinGNO = NULL;
    zhg->Pin_Procs = NULL;

    return ZOLTAN_OK;
  }
  else if (nkeep == 0){

    *edgeGNO = NULL;
    *edgeSize = NULL;
    *edgeWeight = NULL;
    *pinGNO = NULL;
    *pinProcs = NULL;

    if (!save_removed){
      ZOLTAN_FREE(&zhg->edgeGNO);
      ZOLTAN_FREE(&zhg->Esize);
      ZOLTAN_FREE(&zhg->Ewgt);
      ZOLTAN_FREE(&zhg->pinGNO);
      ZOLTAN_FREE(&zhg->Pin_Procs);
      zhg->nHedges = 0;
      zhg->globalHedges = 0;
      zhg->nPins = 0;
      zhg->globalPins = 0;
    }

    return ZOLTAN_OK;
  }

  if (save_removed){
    goEdgeGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * nremove);
    goEdgeSize = (int *)ZOLTAN_MALLOC(sizeof(int) * nremove);
    goPinGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * nremove_size);
    goPinProc = (int *)ZOLTAN_MALLOC(sizeof(int) * nremove_size);

    if (!goEdgeGNO || !goEdgeSize || !goPinGNO || !goPinProc){
      MEMORY_ERROR;
    }

    goEdgeWeight = (float *)ZOLTAN_MALLOC(sizeof(float) * nremove * ew_dim);

    if (ew_dim && !goEdgeWeight){
      MEMORY_ERROR;
    }
  }

  keepEdgeGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * nkeep);
  keepEdgeSize = (int *)ZOLTAN_MALLOC(sizeof(int) * nkeep);
  keepPinGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * nkeep_size);
  keepPinProc = (int *)ZOLTAN_MALLOC(sizeof(int) * nkeep_size);

  if (!keepEdgeGNO || !keepEdgeSize || !keepPinGNO || !keepPinProc){
    MEMORY_ERROR;
  }

  keepEdgeWeight = (float *)ZOLTAN_MALLOC(sizeof(float) * nkeep * ew_dim);

  if (ew_dim && !keepEdgeWeight){
    MEMORY_ERROR;
  }

  gnoPtr = zhg->pinGNO;
  procPtr = zhg->Pin_Procs;
  wgtPtr = zhg->Ewgt;

  for (i=0, k=0, l=0, kpin=0, lpin=0; i < zhg->nHedges; i++){
    esize = zhg->Esize[i];
    if ((esize > gesize_threshold) || (esize == 0)){
      if (save_removed){
        goEdgeGNO[k] = zhg->edgeGNO[i];
        goEdgeSize[k] = esize;
        for (w=0; w < ew_dim; w++){
          goEdgeWeight[k*ew_dim + w] = *wgtPtr++;
        }
        for (pin=0; pin < esize; pin++){
          goPinGNO[kpin + pin] = *gnoPtr++;
          goPinProc[kpin + pin] = *procPtr++;
        }
        kpin += esize;
        k++;
      }
      else{
        wgtPtr += ew_dim;
        gnoPtr += esize;
        procPtr += esize;
      }
    }
    else{
      keepEdgeGNO[l] = zhg->edgeGNO[i];
      keepEdgeSize[l] = esize;
      for (w=0; w < ew_dim; w++){
        keepEdgeWeight[l*ew_dim + w] = *wgtPtr++;
      }
      for (pin=0; pin < esize; pin++){
        keepPinGNO[lpin + pin] = *gnoPtr++;
        keepPinProc[lpin + pin] = *procPtr++;
      }
      lpin += esize;
      l++;
    }
  }

  /* kept edges */

  *edgeGNO = keepEdgeGNO;
  *edgeSize = keepEdgeSize;
  *edgeWeight = keepEdgeWeight;
  *pinGNO = keepPinGNO;
  *pinProcs = keepPinProc;

  ZOLTAN_FREE(&zhg->edgeGNO);
  ZOLTAN_FREE(&zhg->Esize);
  ZOLTAN_FREE(&zhg->Ewgt);
  ZOLTAN_FREE(&zhg->pinGNO);
  ZOLTAN_FREE(&zhg->Pin_Procs);
  zhg->nHedges = 0;
  zhg->globalHedges = 0;
  zhg->nPins = 0;
  zhg->globalPins = 0;

  /* removed edges */

  if (save_removed){
    zhg->nHedges = nremove;
    zhg->globalHedges = global_nremove;
    zhg->nPins = nremove_size;
    zhg->globalPins = global_nremove_pins;
    zhg->edgeGNO = goEdgeGNO;
    zhg->Esize = goEdgeSize;
    zhg->Ewgt = goEdgeWeight;
    zhg->pinGNO = goPinGNO;
    zhg->Pin_Procs = goPinProc;
  }

End:

  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
