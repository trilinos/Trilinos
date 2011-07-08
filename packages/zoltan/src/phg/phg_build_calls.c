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

#include <math.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "phg.h"
#include "phg_lookup.h"
#include "third_library_const.h"
#include "third_library_tools.h"

/* #define CEDRIC_PRINT */

static int edge_weight_operation(ZZ *zz, float *dest, float *src, int ew_dim, int ew_op, int len);

static int Convert_To_CSR( ZZ *zz, int num_pins, int *col_ptr,
    int *num_lists, ZOLTAN_ID_PTR *vtx_GID, int **row_ptr, ZOLTAN_ID_PTR *edg_GID);

/*****************************************************************************/
int Zoltan_Get_Hypergraph_From_Queries(
                    ZZ *zz,                      /* input zoltan struct */
                    PHGPartParams *hgp,          /* input phg parameters */
                    int hgraph_model,            /* input model (graph/hgraph)
                                                    to use in construction */
                    ZHG *zhg)                    /* output hypergraph */
{
static char *yo = "Zoltan_Get_Hypergraph_From_Queries";
MPI_Comm comm = zz->Communicator;
int gid_size = zz->Num_GID;
int gid_chars = zz->Num_GID * sizeof(ZOLTAN_ID_TYPE);
int lid_size = zz->Num_LID;
int nProc = zz->Num_Proc;
int ew_dim = zz->Edge_Weight_Dim;

int hypergraph_callbacks = 0;
int graph_callbacks = 0;
int use_all_neighbors = -1;
int need_pin_weights = 0;
ZOLTAN_MAP *map1=NULL, *map2=NULL;
int msg_tag = 30000;
int ierr = ZOLTAN_OK;

ZOLTAN_COMM_OBJ *plan=NULL;

ZOLTAN_GNO_TYPE *sendGnoBuf = NULL, *recvGnoBuf = NULL;
ZOLTAN_ID_TYPE *sendIdBuf = NULL, *recvIdBuf = NULL;

int *procBuf= NULL;
ZOLTAN_GNO_TYPE *edgeBuf= NULL;
int *pinIdx = NULL;

float *gid_weights = NULL;
float *wgts = NULL;
float *calcVwgt = NULL;
float *sendFloatBuf= NULL;
float *recvFloatBuf= NULL;

char *flag = NULL;

ZOLTAN_ID_PTR gid_buf = NULL;
ZOLTAN_ID_PTR pin_gid_buf = NULL;
ZOLTAN_ID_PTR global_ids = NULL;
ZOLTAN_ID_PTR ew_lids = NULL;
ZOLTAN_ID_PTR gid_requests = NULL;

int add_vweight, nEdge;
int ew_op;
int randomizeInitDist, add_obj_weight;
int i, j, w, k, cnt, dim, rc;
int nRequests;

float *src, *dest;
float *fromwgt, *towgt;

float weight_val;
intptr_t index, iptr;   /* integers the same size as a pointer */
ZOLTAN_GNO_TYPE gnos[2];
ZOLTAN_GNO_TYPE tmpgno;
MPI_Datatype zoltan_gno_mpi_type;

int gno_size_for_dd;

ZOLTAN_ID_PTR fromID, toID;

zoltan_objects       myObjs;
zoltan_pins          myPins;
zoltan_ews           myEWs;
zoltan_temp_edges    myHshEdges;
zoltan_temp_vertices myHshVtxs;

phg_GID_lookup       *lookup_myObjs = NULL;
phg_GID_lookup       *lookup_myHshEdges = NULL;
phg_GID_lookup       *lookup_myHshVtxs = NULL;

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  /* Use the graph or hypergraph query functions to build the hypergraph.
   *
   * Assign a consecutive global numbering system to edges and
   * vertices (objects).  If ADD_OBJ_WEIGHT is set, add the requested weight.
   * When combining hyperedge weights from more than one process, use the
   * PHG_EDGE_WEIGHT_OPERATION parameter.
   *
   * This is the hypergraph intended by the application, before alterations
   * for the phg algorithm.  This hypergraph is suitable for Zoltan_LB_Eval.
   *
   * The hyperedges returned are whole hyperedges - they are not distributed 
   * over more than one process.
   */

  ZOLTAN_TRACE_ENTER(zz, yo);

  gno_size_for_dd = sizeof(ZOLTAN_GNO_TYPE) / sizeof(ZOLTAN_ID_TYPE);

  /* initialize temporary search structures */

  memset(&myObjs, 0, sizeof(zoltan_objects));
  memset(&myPins, 0, sizeof(zoltan_pins));
  memset(&myEWs, 0, sizeof(zoltan_ews));
  memset(&myHshEdges, 0, sizeof(zoltan_temp_edges));
  memset(&myHshVtxs, 0, sizeof(zoltan_temp_vertices));

#ifdef CEDRIC_2D_PARTITIONS
  zhg->ddHedge = NULL;
#endif

  /* get parameters */

  randomizeInitDist = hgp->RandomizeInitDist;
  ew_op = hgp->edge_weight_op;
  add_obj_weight = hgp->add_obj_weight;

  if (add_obj_weight == PHG_ADD_PINS_WEIGHT){
    need_pin_weights = 1;
  }

  if (hgraph_model == HYPERGRAPH){      /* "neighbors" */
    use_all_neighbors = 1;
  }
  else if (hgraph_model == GRAPH){      /* "pairs"     */
    use_all_neighbors = 0;
  }
  else{
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid hgraph_model option");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->Get_HG_Size_CS && zz->Get_HG_CS){
    hypergraph_callbacks = 1;
  }
  if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
           (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) {
    graph_callbacks = 1;
  }

  if (graph_callbacks && hypergraph_callbacks){
    if (hgraph_model == GRAPH)
      hypergraph_callbacks = 0;
  }

  /**************************************************/
  /* Obtain vertex information from the application */
  /**************************************************/

  ierr = Zoltan_Get_Obj_List(zz, &zhg->nObj, &zhg->objGID, &zhg->objLID,
                             zz->Obj_Weight_Dim, &zhg->objWeight,
                             &zhg->Input_Parts);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
    goto End;
  }

  if (zhg->nObj){
    /* Zoltan internal hypergraph will use sequential global numbers */
    zhg->objGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * zhg->nObj);
    if (!zhg->objGNO) MEMORY_ERROR;
  }

  /*******************************************************************/
  /* Assign vertex consecutive numbers (gnos)                        */
  /*******************************************************************/

  ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, zhg->objGNO, zhg->nObj,
                     randomizeInitDist, &zhg->globalObj);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error assigning global numbers to vertices");
    goto End;
  }

  /***********************************************************************/
  /* Get hyperedge information from application through query functions. */
  /***********************************************************************/

  if (hypergraph_callbacks){

    /* Hash each vertex GID to a process.  */

    ierr = phg_map_GIDs_to_processes(zz, zhg->objGID, zhg->nObj, gid_size, &myObjs.vtxHash, nProc);

    if (ierr != ZOLTAN_OK){
      goto End;
    }

    /* Get the vertices hashed to me, with global numbers and owner process IDs */

    ierr = Zoltan_Comm_Create(&plan, zhg->nObj, myObjs.vtxHash, comm, msg_tag, &myHshVtxs.size);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    phg_free_objects(&myObjs);

    if (myHshVtxs.size > 0){
      myHshVtxs.vtxGID = ZOLTAN_MALLOC_GID_ARRAY(zz, myHshVtxs.size);
      myHshVtxs.vtxOwner = (int *)ZOLTAN_MALLOC(sizeof(int) * myHshVtxs.size);
      myHshVtxs.vtxGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * myHshVtxs.size);

      if (!myHshVtxs.vtxGID || !myHshVtxs.vtxOwner || !myHshVtxs.vtxGNO) {
        MEMORY_ERROR;
      }
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)zhg->objGID, gid_chars, (char *)myHshVtxs.vtxGID);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)zhg->objGNO, sizeof(ZOLTAN_GNO_TYPE), (char *)myHshVtxs.vtxGNO);

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

    /* Create search structure for looking up vertex GIDs that were hashed to me */

    lookup_myHshVtxs = phg_create_GID_lookup_table(myHshVtxs.vtxGID, myHshVtxs.size, gid_size);

    if (!lookup_myHshVtxs) MEMORY_ERROR;

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
     * We assume that no two processes will supply the same pin.
     * But more than one process may supply pins for the same edge.
     */

    ierr = Zoltan_Hypergraph_Queries(zz, &myPins.nHedges,
                 &myPins.numPins, &myPins.edgeGID, &pinIdx,
                 &myPins.pinGID);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    if (myPins.nHedges){
      myPins.esizes = (int *)ZOLTAN_MALLOC(sizeof(int) * myPins.nHedges);
      if (!myPins.esizes) MEMORY_ERROR;
    }

    for (i=0; i<myPins.nHedges; i++){
      myPins.esizes[i] = pinIdx[i+1] - pinIdx[i];
    }

    ZOLTAN_FREE(&pinIdx);

    /*
     * Assign each edge global ID to a process using a hash function.
     * Let that process gather all the pin information for the edge.
     * This process will also return that edge in the ZHG structure.
     */

    ierr = phg_map_GIDs_to_processes(zz, myPins.edgeGID, myPins.nHedges,
               gid_size, &myPins.edgeHash, nProc);

    if (ierr != ZOLTAN_OK){
      goto End;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Create(&plan, myPins.nHedges, myPins.edgeHash, comm, msg_tag, &nRequests);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    /* send edge size and edge GID together */

    cnt = 1 + gid_size;

    if (nRequests > 0){
      recvIdBuf = (ZOLTAN_ID_TYPE *)ZOLTAN_MALLOC(cnt * sizeof(ZOLTAN_ID_TYPE) * nRequests);
      if (!recvIdBuf) MEMORY_ERROR;
    }
    if (myPins.nHedges > 0){
      sendIdBuf = (ZOLTAN_ID_TYPE *)ZOLTAN_MALLOC(cnt * sizeof(ZOLTAN_ID_TYPE) * myPins.nHedges);
      if (!sendIdBuf) MEMORY_ERROR;

      for (i=0, j=0, k=0; i < myPins.nHedges; i++, j+= cnt, k += gid_size){

        sendIdBuf[j] = (ZOLTAN_ID_TYPE)myPins.esizes[i];
        ZOLTAN_SET_GID(zz, sendIdBuf + j + 1, myPins.edgeGID + k);
      }
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)sendIdBuf, sizeof(ZOLTAN_ID_TYPE) * cnt, (char *)recvIdBuf);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    ZOLTAN_FREE(&sendIdBuf);

    if (nRequests > 0){
      global_ids = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
      gid_buf = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);

      if (!global_ids || !gid_buf) MEMORY_ERROR;

      for (i=0, j=0, k=0; i < nRequests; i++, j+= cnt, k += gid_size){
        ZOLTAN_SET_GID(zz, global_ids + k, recvIdBuf + j + 1);
        if (i){
          recvIdBuf[i] = recvIdBuf[j];  /* edge size */
        }
      }
      /* need to save a copy of original global_ids array */
      memcpy(gid_buf, global_ids, sizeof(ZOLTAN_ID_TYPE)*nRequests*gid_size);
    }

    /*
     * Create a search structure allowing me to look up info about
     * any of the edge GIDs that were just assigned to me.
     * Rewrites global_ids with a list of unique GIDs.
     */

    lookup_myHshEdges = phg_create_GID_lookup_table2(global_ids, nRequests, gid_size);

    zhg->nHedges = lookup_myHshEdges->numGIDs;
    zhg->nPins = 0;
    myHshEdges.edgeGID = global_ids;

    if (zhg->nHedges > 0){
      zhg->Esize  = (int *)ZOLTAN_CALLOC(zhg->nHedges , sizeof(int));
      zhg->edgeGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(zhg->nHedges * sizeof(ZOLTAN_GNO_TYPE));
      if (!zhg->Esize || !zhg->edgeGNO) MEMORY_ERROR;
    }

    for (i=0; i<nRequests; i++){
      j = phg_lookup_GID(lookup_myHshEdges, gid_buf + (i*gid_size));
      if (j < 0) FATAL_ERROR("Invalid global edge ID received");
      zhg->Esize[j] += (int)recvIdBuf[i];
    }

    for (j=0; j < zhg->nHedges; j++){
      zhg->nPins += zhg->Esize[j];
    }

#ifdef CEDRIC_2D_PARTITIONS
    if (hgp->keep_tree) {
      ZOLTAN_GNO_TYPE offset;
      ZOLTAN_GNO_TYPE *egno = NULL;
      MPI_Scan(&zhg->nHedges, &offset, 1, zoltan_gno_mpi_type, MPI_SUM, zz->Communicator);
      offset -= zhg->nHedges;
#ifdef CEDRIC_PRINT
      for (j=0; j < zhg->nHedges; j++){
	fprintf (stderr, "EDGEGID " ZOLTAN_ID_SPEC "\t%zd\n", global_ids[j], offset + j);
      }
#endif /* CEDRIC_PRINT */

      egno = (ZOLTAN_GNO_TYPE*)ZOLTAN_MALLOC(zhg->nHedges*sizeof(ZOLTAN_GNO_TYPE));
      if (zhg->nHedges && !egno) MEMORY_ERROR;

      for (j=0; j < zhg->nHedges; j++){
	egno[j] = j + offset;
      }
      Zoltan_DD_Create (&zhg->ddHedge, zz->Communicator, gno_size_for_dd, zz->Num_GID,
			0, zhg->nHedges, 0);
      Zoltan_DD_Update (zhg->ddHedge, (ZOLTAN_ID_PTR) egno, global_ids, NULL,
			  NULL, zhg->nHedges);
      ZOLTAN_FREE(&egno);
    }
#endif /* CEDRIC_2D_PARTITIONS */

    /* Get edge pins.  */

    if (zhg->nHedges){
      pinIdx = (int *)ZOLTAN_MALLOC( zhg->nHedges * sizeof(int));
      if (!pinIdx) MEMORY_ERROR;
    }
    if (zhg->nPins){
      pin_gid_buf = ZOLTAN_MALLOC_GID_ARRAY(zz, zhg->nPins);
      if (!pin_gid_buf) MEMORY_ERROR;
      myHshEdges.pinGID = ZOLTAN_MALLOC_GID_ARRAY(zz, zhg->nPins);
      if (!myHshEdges.pinGID) MEMORY_ERROR;
      zhg->pinGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(zhg->nPins * sizeof(ZOLTAN_GNO_TYPE));
      if (!zhg->pinGNO) MEMORY_ERROR;
      zhg->Pin_Procs = (int *)ZOLTAN_MALLOC(zhg->nPins * sizeof(int));
      if (!zhg->Pin_Procs) MEMORY_ERROR;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Resize(plan, myPins.esizes, msg_tag, &cnt);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myPins.pinGID, gid_chars, (char *)pin_gid_buf);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    Zoltan_Comm_Destroy(&plan);
    phg_free_pins(&myPins);

    if (zhg->nHedges > 0){
      pinIdx[0] = 0;
      for (i=1; i < zhg->nHedges; i++){
        pinIdx[i] = pinIdx[i-1] + zhg->Esize[i-1];
      }
    }

    fromID = pin_gid_buf;

    for (i=0; i < nRequests; i++){

      j = phg_lookup_GID(lookup_myHshEdges, gid_buf + (i*gid_size));
      if (j < 0) FATAL_ERROR("Invalid global edge ID received");

      toID = myHshEdges.pinGID + (pinIdx[j] * gid_size);

      for (k=0; k < recvIdBuf[i]; k++){
        ZOLTAN_SET_GID(zz, toID , fromID);
        toID += gid_size;
        fromID += gid_size;
      }

      pinIdx[j] += (int)recvIdBuf[i];
    }

    ZOLTAN_FREE(&recvIdBuf);
    ZOLTAN_FREE(&pinIdx);
    ZOLTAN_FREE(&pin_gid_buf);
    ZOLTAN_FREE(&gid_buf);

    /* For each edge, get the pin vertex global number and the process that owns the pin. */

    ierr = phg_map_GIDs_to_processes(zz, myHshEdges.pinGID, zhg->nPins, gid_size, &myHshEdges.pinHash, nProc);

    if ((ierr!=ZOLTAN_OK) && (ierr!=ZOLTAN_WARN)){
       goto End;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Create(&plan, zhg->nPins, myHshEdges.pinHash, comm, msg_tag, &nRequests);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    if (nRequests){
      gid_buf = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
      if (!gid_buf) MEMORY_ERROR;
      sendGnoBuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(nRequests * sizeof(ZOLTAN_GNO_TYPE) * 2);
      if (!sendGnoBuf) MEMORY_ERROR;
    }

    msg_tag--;

    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myHshEdges.pinGID, gid_chars, (char *)gid_buf);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    for (i=0; i<nRequests; i++){
      j = phg_lookup_GID(lookup_myHshVtxs, gid_buf + ( i * gid_size));
      if (j < 0) FATAL_ERROR("Unexpected vertex GID received");

      sendGnoBuf[2*i] = (ZOLTAN_GNO_TYPE)myHshVtxs.vtxOwner[j];
      sendGnoBuf[2*i + 1] = myHshVtxs.vtxGNO[j];
    }

    ZOLTAN_FREE(&gid_buf);
    phg_free_temp_vertices(&myHshVtxs);
    phg_free_GID_lookup_table(&lookup_myHshVtxs);

    if (zhg->nPins > 0){
      recvGnoBuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * 2 * zhg->nPins);
      if (!recvGnoBuf) MEMORY_ERROR;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *)sendGnoBuf, 
                     sizeof(ZOLTAN_GNO_TYPE) * 2, NULL, (char *)recvGnoBuf);

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    Zoltan_Comm_Destroy(&plan);
    ZOLTAN_FREE(&sendGnoBuf);

    for (i =0; i < zhg->nPins; i++){
      zhg->Pin_Procs[i] = (int)recvGnoBuf[2*i];
      zhg->pinGNO[i] = recvGnoBuf[2*i + 1];
    }

    ZOLTAN_FREE(&recvGnoBuf);

    if (need_pin_weights){

      /* Send to vertex owner the number of edges containing that vertex */

      map1 = Zoltan_Map_Create(zz, 0, sizeof(ZOLTAN_GNO_TYPE), 0, zhg->nObj);
      if (map1 == NULL) goto End;

      for (iptr=0; iptr < zhg->nObj; iptr++){
        ierr = Zoltan_Map_Add(zz, map1, (char *)(zhg->objGNO + iptr), iptr + 1);
        if (ierr != ZOLTAN_OK) goto End;
      }

      msg_tag--;
      ierr = Zoltan_Comm_Create(&plan, zhg->nPins, zhg->Pin_Procs, comm, msg_tag, &nRequests);

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }

      if (nRequests){
        recvGnoBuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(nRequests * sizeof(ZOLTAN_GNO_TYPE));
        if (!recvGnoBuf) MEMORY_ERROR;
      }

      msg_tag--;
      ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)zhg->pinGNO, sizeof(ZOLTAN_GNO_TYPE), (char *)recvGnoBuf);

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }

      Zoltan_Comm_Destroy(&plan);

      zhg->numHEdges = (int *)ZOLTAN_CALLOC(sizeof(int), zhg->nObj);
      if (zhg->nObj && !zhg->numHEdges) MEMORY_ERROR;

      for (i=0; i < nRequests; i++){
        ierr = Zoltan_Map_Find(zz, map1, (char *)(recvGnoBuf + i), &iptr);
        if (ierr != ZOLTAN_OK) goto End;
        if (iptr == ZOLTAN_NOT_FOUND) FATAL_ERROR("Unexpected vertex global number received");

        index = iptr - 1;
        zhg->numHEdges[index]++;
      }

      ZOLTAN_FREE(&recvGnoBuf);

      Zoltan_Map_Destroy(zz, &map1);
    }

    zhg->edgeWeightDim = ew_dim;

    if (ew_dim && zz->Get_HG_Size_Edge_Wts && zz->Get_HG_Edge_Wts){

      /* Get edge weights */

      zz->Get_HG_Size_Edge_Wts(zz->Get_HG_Size_Edge_Wts_Data, &myEWs.size, &ierr);

      if ((ierr!=ZOLTAN_OK) && (ierr!=ZOLTAN_WARN)){
        FATAL_ERROR("obtaining edge weight size");
      }

      cnt = sizeof(float) * ew_dim;

      if (myEWs.size > 0){

        myEWs.edgeGID = ZOLTAN_MALLOC_GID_ARRAY(zz, myEWs.size);
        ew_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, myEWs.size);
        myEWs.wgt = (float *)ZOLTAN_MALLOC(myEWs.size * cnt);

        if (!myEWs.edgeGID || !ew_lids || !myEWs.wgt){
          MEMORY_ERROR;
        }

        zz->Get_HG_Edge_Wts(zz->Get_HG_Edge_Wts_Data, gid_size, lid_size, myEWs.size, ew_dim,
                    myEWs.edgeGID, ew_lids, myEWs.wgt, &ierr);

        if ((ierr!=ZOLTAN_OK) && (ierr!=ZOLTAN_WARN)){
          FATAL_ERROR("obtaining edge weights");
        }

        ZOLTAN_FREE(&ew_lids);

        /* Get process assigned to each hyperedge.  */

        ierr = phg_map_GIDs_to_processes(zz, myEWs.edgeGID, myEWs.size, gid_size, 
                  &myEWs.edgeHash, nProc);

        if ((ierr!=ZOLTAN_OK) && (ierr!=ZOLTAN_WARN)){
          goto End;
        }
      }

      /* Send all edge weights to the process that was assigned that edge */

      gid_requests = NULL;
      gid_weights = NULL;

      msg_tag--;
      ierr = Zoltan_Comm_Create(&plan, myEWs.size, myEWs.edgeHash, comm, msg_tag, &nRequests);

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }

      if (nRequests > 0){
        gid_requests = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
        gid_weights = (float *)ZOLTAN_MALLOC(cnt*nRequests);

        if ( (ew_dim && !gid_weights) || !gid_requests) MEMORY_ERROR;
      }

      msg_tag--;
      ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myEWs.edgeGID, gid_chars, (char *)gid_requests);

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }

      msg_tag--;
      ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myEWs.wgt, cnt, (char *)gid_weights);

      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }

      Zoltan_Comm_Destroy(&plan);
      phg_free_ews(&myEWs);

      if (nRequests > 0){

        /* Combine edge weights supplied for the same edge using the PHG_EDGE_WEIGHT_OPERATION */

        zhg->Ewgt = (float *)ZOLTAN_MALLOC(cnt * zhg->nHedges);

        flag = (char *)ZOLTAN_MALLOC(sizeof(char) * zhg->nHedges);

        if (zhg->nHedges && (!flag || !zhg->Ewgt))  MEMORY_ERROR;

        memset(flag, 0, zhg->nHedges * sizeof(char));

        src = gid_weights;

        for (i=0; i < nRequests; i++, src += ew_dim){
          j = phg_lookup_GID(lookup_myHshEdges, gid_requests + (i*gid_size));

          if (j < 0){
            continue; /* An edge for which there are no pins, ignore it. */
          }

          dest = zhg->Ewgt + (j * ew_dim);

          if (!flag[j]){
            for (w=0; w<ew_dim; w++){
              dest[w] = src[w];
            }
            flag[j] = 1;
          }
          else{
            ierr = edge_weight_operation(zz, dest, src, ew_dim, ew_op, 1);
            if (ierr != ZOLTAN_OK) goto End;
          }
        }
        ZOLTAN_FREE(&gid_requests);
        ZOLTAN_FREE(&gid_weights);
        ZOLTAN_FREE(&flag);
      }
    } /* end of get hyperedge weights */

    phg_free_GID_lookup_table(&lookup_myHshEdges);
    phg_free_temp_edges(&myHshEdges);             /* this frees global_ids array */

  } else if (graph_callbacks){

    ierr = Zoltan_Graph_Queries(zz, zhg->nObj, zhg->objGID, zhg->objLID,
      &zhg->nPins,        /* total of number of neighbors of each vertex */
      &zhg->Esize,        /* number of neighbors of each vertex */
      &myPins.pinGID,     /* vertex GID for each neighbor */
      &zhg->Pin_Procs,    /* process owning each pinGID */
      &zhg->Ewgt);        /* weight of edge from vtx to neighbor */

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error");
      goto End;
    }

    /* Create a search structure to lookup vertex information by GID */

    lookup_myObjs = phg_create_GID_lookup_table(zhg->objGID, zhg->nObj, gid_size);
    if (!lookup_myObjs) MEMORY_ERROR;

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error");
      goto End;
    }

    /* Get the global number of each my pin vertices */

    msg_tag--;
    ierr = Zoltan_Comm_Create(&plan, zhg->nPins, zhg->Pin_Procs, comm, msg_tag, &nRequests);
  
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }
  
    if (nRequests > 0){
      gid_buf = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
      if (!gid_buf) MEMORY_ERROR;
      sendGnoBuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * nRequests);
      if (!sendGnoBuf) MEMORY_ERROR;
    }
  
    msg_tag--;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)myPins.pinGID, gid_chars, (char *)gid_buf);
  
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    phg_free_pins(&myPins);

    zhg->numHEdges = (int *)ZOLTAN_CALLOC(sizeof(int), zhg->nObj);
    if (zhg->nObj && !zhg->numHEdges) MEMORY_ERROR;
  
    for (i=0; i<nRequests; i++){
      j = phg_lookup_GID(lookup_myObjs, gid_buf + ( i * gid_size));
      if (j < 0) FATAL_ERROR("Unexpected vertex GID received");
      sendGnoBuf[i] = zhg->objGNO[j];
      zhg->numHEdges[j]++;  /* number of edges this vertex is in, in original graph */
    }

    ZOLTAN_FREE(&gid_buf);
    phg_free_GID_lookup_table(&lookup_myObjs);
    phg_free_objects(&myObjs);
  
    msg_tag--;

    zhg->pinGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_CALLOC(sizeof(ZOLTAN_GNO_TYPE) , zhg->nPins);
    if (zhg->nPins && !zhg->pinGNO) MEMORY_ERROR;
  
    ierr = Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *)sendGnoBuf, sizeof(ZOLTAN_GNO_TYPE),
                  NULL, (char *)zhg->pinGNO);
  
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    ZOLTAN_FREE(&sendGnoBuf);
    Zoltan_Comm_Destroy(&plan);

    /*
     * There may be two edge weights provided for a single edge, depending on 
     * how the application supplied the graph in the query functions.  
     * Find the one or two weights supplied for each edge.
     * If there are two weights, combine them as specified by the 
     * PHG_EDGE_WEIGHT_OPERATION.
     */

    /* Create a lookup object for all vertex pairs in graph */

    map2 = Zoltan_Map_Create(zz, 0, 2 * sizeof(ZOLTAN_GNO_TYPE), 1, zhg->nPins);
    if (map2 == NULL) goto End;

    for (i=0, iptr=0; i < zhg->nObj; i++){
      for (j=0; j < zhg->Esize[i]; j++, iptr++){
         gnos[0] = zhg->objGNO[i];
         gnos[1] = zhg->pinGNO[iptr];
         ierr = Zoltan_Map_Add(zz, map2, (char *)gnos, iptr + 1);
         if (ierr != ZOLTAN_OK) goto End;
      }
    }

    /* Search for a second weight for each edge, either on process or off process.
     * Most by far should be on process, so we'll look there first instead of
     * using collective communication for all vertex pairs.
     *
     * Even if we have no edge weights, we do this to find out about edges that
     * appear twice in the graph supplied by the application.
     */

    if (!use_all_neighbors){
      /* flag each edge that appears more than once */
      flag = (char *)ZOLTAN_CALLOC(sizeof(char), zhg->nPins);
      if (zhg->nPins && !flag) MEMORY_ERROR;
    }

    /* wgts will hold the other weight provided for the edge, if any */
    wgts = (float *)ZOLTAN_CALLOC(sizeof(float) , zhg->nPins * ew_dim);
    if (zhg->nPins && ew_dim && !wgts) MEMORY_ERROR;

    cnt = 0;

    for (i=0, k=0; i < zhg->nObj; i++){
      gnos[1] = zhg->objGNO[i];
      for (j=0; j < zhg->Esize[i]; j++, k++){
         gnos[0] = zhg->pinGNO[k];
         if (zhg->Pin_Procs[k] != zz->Proc){
           /* if there's another weight, it's on a different proc */
           cnt++;
         }
         else{
           ierr = Zoltan_Map_Find(zz, map2, (char *)gnos, &iptr);
           if (ierr != ZOLTAN_OK) goto End;

           if (iptr != ZOLTAN_NOT_FOUND){
             if (ew_dim){
               /* this proc provided weights for [v0,v1] and [v1, v0] */
               index = iptr - 1;
               dest = wgts + k * ew_dim; 
               src = zhg->Ewgt + index * ew_dim;
               for (dim = 0; dim < ew_dim; dim++){
                 dest[dim] = src[dim];
               }
             }
             if (!use_all_neighbors){
               flag[k] = 1;   /* flag that edge appears twice in graph */
             }
           }
           else{ /* I don't have a duplicate edge on my proc */
           }
         }
      }
    }

    if (cnt > 0){
      sendGnoBuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * 2 * cnt);
      if (!sendGnoBuf) MEMORY_ERROR;
      procBuf = (int *)ZOLTAN_MALLOC(sizeof(int) * cnt);
      if (!procBuf) MEMORY_ERROR;
      cnt = 0;

      for (i=0, k=0; i < zhg->nObj; i++){
        for (j=0; j < zhg->Esize[i]; j++, k++){
           if (zhg->Pin_Procs[k] != zz->Proc){
             procBuf[cnt] = zhg->Pin_Procs[k];
             sendGnoBuf[2*cnt] = zhg->pinGNO[k];
             sendGnoBuf[2*cnt + 1] = zhg->objGNO[i];
             cnt++;
           }
        }
      }
    }

    msg_tag--;
    ierr = Zoltan_Comm_Create(&plan, cnt, procBuf, comm, msg_tag, &nRequests);
  
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    ZOLTAN_FREE(&procBuf);
  
    if (nRequests > 0){
      sendFloatBuf = (float *) ZOLTAN_CALLOC(sizeof(float) , nRequests * ew_dim);
      if (ew_dim && !sendFloatBuf) MEMORY_ERROR;
      recvGnoBuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * 2 * nRequests);
      if (!recvGnoBuf) MEMORY_ERROR;
    }

    msg_tag--;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *)sendGnoBuf, sizeof(ZOLTAN_GNO_TYPE) * 2, (char *)recvGnoBuf);
  
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      goto End;
    }

    ZOLTAN_FREE(&sendGnoBuf);
  
    nEdge = 0;
    for (i=0; i<nRequests; i++){
      gnos[0] = recvGnoBuf[i*2];
      gnos[1] = recvGnoBuf[i*2 + 1];

      ierr = Zoltan_Map_Find(zz, map2, (char *)gnos, &iptr);
      if (ierr != ZOLTAN_OK) goto End;

      index = iptr - 1;

      if (iptr != ZOLTAN_NOT_FOUND){
        if (ew_dim){
          dest = sendFloatBuf + i * ew_dim; 
          src = zhg->Ewgt + index * ew_dim;
          for (dim = 0; dim < ew_dim; dim++){
            dest[dim] = src[dim];
          }
        }
        nEdge = 1;   /* flag that I have weights to send back */

        if (!use_all_neighbors){
          flag[index] = 1;   /* flag that edge appears twice in graph */
        }
      }
      else{ /* I don't have an edge corresponding to the edge on the other proc */
      }
    }

    Zoltan_Map_Destroy(zz, &map2);
    ZOLTAN_FREE(&recvGnoBuf);

    if (ew_dim){

      /* Before doing global communication, determine whether there is any information to share.  */
      rc = MPI_Allreduce(&nEdge, &w, 1, MPI_INT, MPI_MAX, comm);
      CHECK_FOR_MPI_ERROR(rc);
  
      if (w > 0){
        if (cnt > 0){
          recvFloatBuf = (float *)ZOLTAN_MALLOC(sizeof(float) * cnt * ew_dim);
          if (ew_dim && !recvFloatBuf) MEMORY_ERROR;
        }
      
        msg_tag--;
  
        ierr = Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *)sendFloatBuf, sizeof(float) * ew_dim,
                      NULL, (char *)recvFloatBuf);
      
        if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
          goto End;
        }
        ZOLTAN_FREE(&sendFloatBuf);
    
        if (cnt > 0){
          cnt = 0;
          for (k=0; k < zhg->nPins; k++){
            if (zhg->Pin_Procs[k] != zz->Proc){
              src = recvFloatBuf + cnt * ew_dim;
              dest = wgts + k * ew_dim;
              for (dim = 0; dim < ew_dim; dim++){
                dest[dim] = src[dim];
              }
              cnt++;
            }
          }
        }
    
        ZOLTAN_FREE(&recvFloatBuf);
      }
      else{
        ZOLTAN_FREE(&sendFloatBuf);
      }
    }

    Zoltan_Comm_Destroy(&plan);

    /* For all edges that have two weights, combine the weights according to
     * the PHG_EDGE_WEIGHT_OPERATION.
     */

    if (ew_dim){
      ierr = edge_weight_operation(zz, zhg->Ewgt, wgts, ew_dim, ew_op, zhg->nPins);
      if (ierr != ZOLTAN_OK) goto End;
    }

    ZOLTAN_FREE(&wgts);

    /************************************************************************* 
     * Convert the graph edges to hyperedges.
     */

    zhg->edgeWeightDim = ((ew_dim > 0) ? ew_dim : 1);

    ew_dim = zhg->edgeWeightDim;

    if (use_all_neighbors){
      /* 
       * Create a hyperedge out of each vertex, containing that vertex
       * and all of its neighbors.  The hyperedge weight will be the
       * maximum (not sum) of the weight of each original graph edge. 
       */

      zhg->nHedges = zhg->nObj;
      zhg->nPins = zhg->nPins + zhg->nObj;

      zhg->edgeGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * zhg->nObj);
      if (zhg->nObj && !zhg->edgeGNO) MEMORY_ERROR;

      wgts = (float *)ZOLTAN_MALLOC(sizeof(float) * ew_dim * zhg->nHedges);
      if (zhg->nHedges && ew_dim && !wgts) MEMORY_ERROR;
      procBuf = (int *)ZOLTAN_MALLOC(sizeof(int) * zhg->nPins);
      if (zhg->nPins && !procBuf) MEMORY_ERROR;
      edgeBuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * zhg->nPins);
      if (zhg->nPins && !edgeBuf) MEMORY_ERROR;

      k = 0;  /* index into new pin arrays */
      w = 0;  /* index into old nbor arrays */

      for (i=0; i < zhg->nHedges; i++){

        if (zhg->Ewgt){
          for (dim=0; dim < ew_dim; dim++){
            weight_val = 0.0;
            src = zhg->Ewgt + (w * ew_dim) + dim;
            for (j=0; j < zhg->Esize[i]; j++){
              /* Compute maximum graph edge weight */
              if (src[j*ew_dim] > weight_val)
                weight_val = src[j*ew_dim];
            }
            wgts[i*ew_dim + dim] = weight_val;
          }
        }
        else{
          wgts[i] = 1.0;
        }

        edgeBuf[k] = zhg->objGNO[i];
        procBuf[k] = zz->Proc;
        k++;

        for (j=0; j < zhg->Esize[i]; j++, k++){
          edgeBuf[k] = zhg->pinGNO[w+j];
          procBuf[k] = zhg->Pin_Procs[w+j];
        }

        w += zhg->Esize[i];

        zhg->Esize[i] = zhg->Esize[i] + 1;
      }

      ZOLTAN_FREE(&zhg->pinGNO);
      ZOLTAN_FREE(&zhg->Pin_Procs);
      ZOLTAN_FREE(&zhg->Ewgt);

      zhg->pinGNO = edgeBuf;
      zhg->Pin_Procs= procBuf;
      zhg->Ewgt = wgts;
    }
    else{
      /* 
       * Create a hyperedge out of each pair of neighboring vertices.  
       * The hyperedge weight will be the weight of the original graph edge.
       *
       * The user may or may not have specified graph edges twice
       * [v0, v1]  and [v1, v0].  We will only include one instance of
       * each edge.  (Edge weights for both have already been combined.)
       */

      nEdge = zhg->nPins;   /* upper bound on number of edges */

      cnt = 0;       /* actual number of hyperedges */

      edgeBuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * nEdge * 2);    /* pin gno */
      procBuf = (int *)ZOLTAN_MALLOC(sizeof(int) * nEdge * 2);    /* pin proc */

      if (nEdge && (!edgeBuf || !procBuf)) MEMORY_ERROR;

      wgts = (float *)ZOLTAN_MALLOC(sizeof(float) * nEdge * ew_dim);  /* edge weight */
      if (nEdge && ew_dim && !wgts) MEMORY_ERROR;

      k = 0;
      for (i=0; i < zhg->nObj; i++){
        gnos[0] = zhg->objGNO[i];
        for (j = 0; j < zhg->Esize[i]; j++, k++){
          gnos[1] = zhg->pinGNO[k];

          /* include each edge only once */

          if ((gnos[0] < gnos[1]) || !flag[k]){

            edgeBuf[cnt * 2]     = gnos[0];
            edgeBuf[cnt * 2 + 1] = gnos[1];

            procBuf[cnt * 2]     = zz->Proc;
            procBuf[cnt * 2 + 1] = zhg->Pin_Procs[k];

            if (zhg->Ewgt){
              src = zhg->Ewgt + (k * ew_dim);
              dest = wgts + (cnt * ew_dim);

              for (w=0; w < ew_dim; w++){
                dest[w] = src[w];
              }
            }
            else{
              wgts[cnt] = 1.0 + flag[k];
            }

            cnt++;
          }
        }
      }

      ZOLTAN_FREE(&flag);

      zhg->nHedges = cnt;
      zhg->nPins = cnt * 2;

      ZOLTAN_FREE(&zhg->edgeGNO);
      ZOLTAN_FREE(&zhg->Esize);
      ZOLTAN_FREE(&zhg->Ewgt);
      ZOLTAN_FREE(&zhg->pinGNO);
      ZOLTAN_FREE(&zhg->Pin_Procs);

      if (cnt > 0){
  
        if (ew_dim > 0){
          zhg->Ewgt = (float *)ZOLTAN_MALLOC(sizeof(float) * cnt * ew_dim);
          if (!zhg->Ewgt) MEMORY_ERROR;
          memcpy(zhg->Ewgt, wgts, ew_dim * cnt * sizeof(float));
          ZOLTAN_FREE(&wgts);
        }
  
        zhg->edgeGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * cnt);
        if (!zhg->edgeGNO) MEMORY_ERROR;
        zhg->Esize = (int *)ZOLTAN_MALLOC(sizeof(int) * cnt);
        if (!zhg->Esize) MEMORY_ERROR;
        zhg->pinGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE) * cnt * 2);
        if (!zhg->pinGNO) MEMORY_ERROR;
        zhg->Pin_Procs = (int *)ZOLTAN_MALLOC(sizeof(int) * cnt * 2);
        if (!zhg->Pin_Procs) MEMORY_ERROR;
  
        memcpy(zhg->pinGNO, edgeBuf, 2 * cnt * sizeof(ZOLTAN_GNO_TYPE));
  
        memcpy(zhg->Pin_Procs, procBuf, 2 * cnt * sizeof(int));
  
        for (i=0; i < cnt; i++){
          zhg->Esize[i] = 2;
        }
      }
      ZOLTAN_FREE(&edgeBuf);
      ZOLTAN_FREE(&procBuf);
    }
  } else {
    /* Partition without edge information?  Or return an error? */
    ZOLTAN_PRINT_WARN(zz->Proc, yo,
       "No edge information provided, partitioning vertices.");
  }

  /******************************************************************************/
  /* Done with graph or hypergraph queries. Impose a global hyperedge numbering */
  /******************************************************************************/

  ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, zhg->edgeGNO, zhg->nHedges,
                     randomizeInitDist, &zhg->globalHedges);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error assigning global numbers to edges");
    goto End;
  }

  /******************************************************************************/
  /* Set global number of pins                                                  */
  /******************************************************************************/

  tmpgno = (ZOLTAN_GNO_TYPE)zhg->nPins;
  rc = MPI_Allreduce(&tmpgno, &zhg->globalPins, 1, zoltan_gno_mpi_type, MPI_SUM, comm);
  CHECK_FOR_MPI_ERROR(rc);

  /***********************************************************************/
  /* If user requested ADD_OBJ_WEIGHT, modify object weights now.        */
  /***********************************************************************/

  if (add_obj_weight != PHG_ADD_NO_WEIGHT){
    add_vweight = 1;    /* may change in future to allow more than one */
    calcVwgt = (float *)ZOLTAN_CALLOC(sizeof(float), zhg->nObj);
    if (zhg->nObj && !calcVwgt) MEMORY_ERROR;

    if (add_obj_weight == PHG_ADD_UNIT_WEIGHT){
      for (i=0; i<zhg->nObj; i++){
        calcVwgt[i] = 1.0;
      }
    }
    else if (add_obj_weight == PHG_ADD_PINS_WEIGHT){
      for (i=0; i<zhg->nObj; i++){
        calcVwgt[i] = zhg->numHEdges[i];
      }
    }
  }
  else{
    add_vweight = 0;
    calcVwgt = NULL;
  }

  zhg->objWeightDim =     /* 1 or greater */
    zz->Obj_Weight_Dim +   /* 0 or more application supplied weights*/
    add_vweight;           /* 0 or 1 additional calculated weights */

  if (zhg->objWeightDim > 1) {
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
        (float *)ZOLTAN_MALLOC(sizeof(float) * zhg->nObj * zhg->objWeightDim);
      if ((zhg->nObj * zhg->objWeightDim) && !wgts){
        MEMORY_ERROR;
      }

      towgt = wgts;
      fromwgt = zhg->objWeight;

      for (i=0; i<zhg->nObj; i++){
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

    ZOLTAN_FREE(&zhg->objWeight);

    zhg->objWeight = wgts;
  }

End:

  ZOLTAN_FREE(&sendGnoBuf);
  ZOLTAN_FREE(&recvGnoBuf);
  ZOLTAN_FREE(&sendIdBuf);
  ZOLTAN_FREE(&recvIdBuf);
  ZOLTAN_FREE(&pinIdx);

  ZOLTAN_FREE(&gid_weights);
  ZOLTAN_FREE(&calcVwgt);
  ZOLTAN_FREE(&sendFloatBuf);
  ZOLTAN_FREE(&recvFloatBuf);

  ZOLTAN_FREE(&flag);

  ZOLTAN_FREE(&gid_buf);
  ZOLTAN_FREE(&pin_gid_buf);
  ZOLTAN_FREE(&ew_lids);
  ZOLTAN_FREE(&gid_requests);

  phg_free_objects(&myObjs);
  phg_free_pins(&myPins);
  phg_free_temp_vertices(&myHshVtxs);
  phg_free_ews(&myEWs);
  phg_free_temp_edges(&myHshEdges);

  phg_free_GID_lookup_table(&lookup_myHshEdges);
  phg_free_GID_lookup_table(&lookup_myHshVtxs);
  phg_free_GID_lookup_table(&lookup_myObjs);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}

/*****************************************************************************/

int Zoltan_PHG_GIDs_to_global_numbers(ZZ *zz, ZOLTAN_GNO_TYPE *gnos, int len, int randomize,
          ZOLTAN_GNO_TYPE *numGlobalObjects)
{
  static char *yo = "Zoltan_PHG_GIDs_to_global_number";
  int ierr = ZOLTAN_OK;
  ZOLTAN_GNO_TYPE *gtotal = NULL;
  int nProc = zz->Num_Proc;
  ZOLTAN_GNO_TYPE *mycnt = NULL, *gcnt = NULL;
  ZOLTAN_GNO_TYPE tmp;
  int i, rc;
  MPI_Comm comm = zz->Communicator;
  MPI_Datatype zoltan_gno_mpi_type;

  ZOLTAN_TRACE_ENTER(zz, yo);

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  /* The application uses global IDs which are unique but arbitrary.  Zoltan_PHG will
   * use global numbers which are consecutive integers beginning with zero.  Convert
   * the application's global IDs to our global numbers.
   */
  
  if (randomize) {
    /* Randomize the global numbers */
    gtotal = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC(3*nProc+1, sizeof(ZOLTAN_GNO_TYPE));
    if (!gtotal) MEMORY_ERROR;

    mycnt  = gtotal + nProc + 1;
    gcnt   = mycnt + nProc;

    /* Compute random processor bin. */
    /* Count how many local vtxs selected processor bin */
    Zoltan_Srand(Zoltan_Rand(NULL)+zz->Proc, NULL);
    for (i = 0; i < len; i++) {
      gnos[i] = (ZOLTAN_GNO_TYPE)Zoltan_Rand_InRange(NULL, nProc);
      mycnt[gnos[i]]++;
    }
    /* Compute prefix of mycnt */
    rc = MPI_Scan(mycnt, gcnt, nProc, zoltan_gno_mpi_type, MPI_SUM, comm);
    CHECK_FOR_MPI_ERROR(rc);

    rc = MPI_Allreduce(mycnt, gtotal, nProc, zoltan_gno_mpi_type, MPI_SUM, comm);
    CHECK_FOR_MPI_ERROR(rc);

    /* Compute first gno for vertices going to each target bin */
    for (tmp = 0, i = 0; i < nProc; i++) {
      gcnt[i] -= mycnt[i];
      tmp += gtotal[i];
      gtotal[i] = tmp - gtotal[i];
    }

    *numGlobalObjects = gtotal[nProc] = tmp;

    /* Assign gnos sequential from gcnt[bin]. */
    for (i=0; i< len; i++) {
      tmp = gnos[i];
      gnos[i] = (ZOLTAN_GNO_TYPE)(gtotal[tmp] + gcnt[tmp]);
      gcnt[tmp]++;
    }
  }
  else {
    /* Linearly order the input vertices */
    gtotal = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC((nProc+1) * sizeof(ZOLTAN_GNO_TYPE));
    if (!gtotal) MEMORY_ERROR;

    /* Construct gtotal[i] = the number of vertices on all procs < i. */
    /* Scan to compute partial sums of the number of objs */

    tmp = (ZOLTAN_GNO_TYPE)len;
    rc = MPI_Scan(&tmp, gtotal, 1, zoltan_gno_mpi_type, MPI_SUM, comm);
    CHECK_FOR_MPI_ERROR(rc);

    /* Gather data from all procs */

    rc = MPI_Allgather (&(gtotal[0]), 1, zoltan_gno_mpi_type, &(gtotal[1]), 1, zoltan_gno_mpi_type, comm);
    CHECK_FOR_MPI_ERROR(rc);
    gtotal[0] = 0;
    *numGlobalObjects = gtotal[nProc];

    for (i=0; i< len; i++) {
      gnos[i] =  (ZOLTAN_GNO_TYPE)gtotal[zz->Proc]+i;
    }
  }

End:

  ZOLTAN_FREE(&gtotal);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
/*****************************************************************************/

int Zoltan_Hypergraph_Queries(ZZ *zz, 
   int *num_lists,         /* output: number of edges */
   int *num_pins,          /* output: total number of pins in edges */
   ZOLTAN_ID_PTR *edg_GID, /* output: list of edge global IDs */
   int **row_ptr,          /* output: loc in vtx_GID for start of each edge */
                           /*         plus num_pins in last element         */
   ZOLTAN_ID_PTR *vtx_GID) /* output: vertex global ID for each pin */
{
static char *yo = "Zoltan_Hypergraph_Queries";
int ierr = ZOLTAN_OK;
int nl, np, format, have_pins, row_storage;
ZOLTAN_ID_PTR vid, eid;
int *rptr=NULL, *cptr=NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);
  /*
   * Call the pin query functions.  Pins may be provided in
   * compressed row storage format or compressed column storage
   * format.  Return compressed rows, converting if necessary.
   */

  *edg_GID = NULL;
  *vtx_GID = NULL;
  *row_ptr = NULL;
  *num_lists = *num_pins = 0;

  if (!zz->Get_HG_Size_CS || !zz->Get_HG_CS){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Hypergraph query functions undefined");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }

  /* Get size and type of compressed pin storage */

  zz->Get_HG_Size_CS(zz->Get_HG_Size_CS_Data, &nl, &np, &format, &ierr);

  if ((format != ZOLTAN_COMPRESSED_EDGE)&&(format != ZOLTAN_COMPRESSED_VERTEX)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Invalid compression format returned in Get_HG_Size_CS");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }
  ZOLTAN_TRACE_DETAIL(zz, yo, "done with Get_HG_Size_CS");

  have_pins = ((nl > 0) && (np > 0));
  row_storage = (format == ZOLTAN_COMPRESSED_EDGE);

  /* Get the hypergraph pins in compressed storage format */

  if (have_pins){
    if (!row_storage){    /* compressed column storage */

      vid = ZOLTAN_MALLOC_GID_ARRAY(zz, nl);
      cptr = (int *)ZOLTAN_MALLOC((nl+1) * sizeof(int));
      eid = ZOLTAN_MALLOC_GID_ARRAY(zz, np);

      if (!vid|| !cptr || !eid){
        Zoltan_Multifree(__FILE__, __LINE__, 3, &vid, &cptr, &eid);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "memory allocation");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_FATAL;
      }
      zz->Get_HG_CS(zz->Get_HG_CS_Data, zz->Num_GID,
               nl, np, format, vid, cptr, eid, &ierr);
      cptr[nl] = np;

      ZOLTAN_TRACE_DETAIL(zz, yo, "done with Get_HG_CS");


      if ((ierr == ZOLTAN_OK) || (ierr == ZOLTAN_WARN)){
        ierr = Convert_To_CSR(zz,
                     np,    /* number of pins doesn't change */
                     cptr,
                     &nl,     /* replace with number of rows */
                     &vid,    /* replace with pins           */
                     &rptr,   /* index into start of each row in vid */
                     &eid);   /* replace with row (edge) GIDs */
      }
      ZOLTAN_FREE(&cptr);
    }
    else{               /* compressed row storage */

      eid = ZOLTAN_MALLOC_GID_ARRAY(zz, nl);
      rptr = (int *)ZOLTAN_MALLOC((nl+1) * sizeof(int));
      vid = ZOLTAN_MALLOC_GID_ARRAY(zz, np);

      if (!vid || !rptr || !eid){
        Zoltan_Multifree(__FILE__, __LINE__, 3, &vid, &rptr, &eid);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "memory allocation");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_FATAL;
      }

      zz->Get_HG_CS(zz->Get_HG_CS_Data, zz->Num_GID,
                 nl, np, format, eid, rptr, vid, &ierr);

      ZOLTAN_TRACE_DETAIL(zz, yo, "done with Get_HG_CS");
      rptr[nl] = np;
    }

    *edg_GID = eid;
    *vtx_GID = vid;
    *row_ptr = rptr;
    *num_lists = nl;
    *num_pins = np;
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}


/*****************************************************************************/
/* This is called "Convert_To_CSR" but it also converts CSR to CSC.  The
 * conversion is symmetric.
 */

static int Convert_To_CSR(
    ZZ *zz, int num_pins, int *col_ptr,   /* input */
    int *num_lists,                       /* rest are input/output */
    ZOLTAN_ID_PTR *vtx_GID,
    int **row_ptr,
    ZOLTAN_ID_PTR *edg_GID)
{
static char *yo = "Convert_To_CSR";
int numVerts = *num_lists;
int numEdges, ierr, ht_size ;
int edg_cnt ;
ZOLTAN_ID_PTR egid, vgid;
int v, e, idx, found, npins;
struct _hash_node {
  ZOLTAN_ID_PTR egid;
  int loc; /* Will be used as counter for #vertices in an edge and later as
            * as bucket pointer to write the vertices. The value is -(#vertices)
            * when used as a counter */
  struct _hash_node *next;
} *hn=NULL, *tmp;
struct _hash_node **hash_table=NULL;
ZOLTAN_ID_PTR edges=NULL, pins=NULL;
int *eIdx=NULL, *vIdx=NULL;
int numGID = zz->Num_GID;

  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = ZOLTAN_OK;

  if (num_pins == 0){
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ierr;
  }

  /*
   * Convert from CCS to CRS (compressed columns to compressed rows)
   * We have the lists of edges for each vertex.  Create lists of
   * vertices for each edge.
   */

  /* SRSR : Guess numVerts == numEdges for hash table size.  */
  ht_size = Zoltan_Recommended_Hash_Size(numVerts) ;

  hash_table =
    (struct _hash_node **)ZOLTAN_CALLOC(ht_size, sizeof(struct _hash_node *));

  if (!hash_table){
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }

  /* For each edge, count how many vertices (pins) it has */

  egid = *edg_GID;
  numEdges = 0;

  for (e=0; e<num_pins; e++){
     idx = Zoltan_Hash(egid, numGID, (unsigned int)ht_size);
     found = 0;
     hn = hash_table[idx];

     while (hn){
       if (ZOLTAN_EQ_GID(zz, hn->egid, egid)){
         hn->loc--;
         found = 1;
         break;
       }
       else{
         hn = hn->next;
       }
     }
     if (!found){
       /* O(numEdges) calls to malloc may be costlier in some platforms. An
        * array of size O(numVerts) allocated first, and the reallocated if
        * numVerts < numEdges could be useful then. SRSR : The performance 
        * improvement in octopi was small.
        */
       hn = (struct _hash_node *)ZOLTAN_MALLOC(sizeof(struct _hash_node));
       if (!hn){
         ierr = ZOLTAN_MEMERR;
         goto End;
       }
       hn->egid = egid;
       hn->loc = -1;
       hn->next = hash_table[idx];
       hash_table[idx] = hn;
       numEdges++;
     }
     egid += numGID;
  }

  /* Create array of indices into location in pin list.
   * Create the corresponding list of unique edge IDs.                          
   */

  vIdx = (int *)ZOLTAN_MALLOC((numEdges+1) * sizeof(int));
  edges = ZOLTAN_MALLOC_GID_ARRAY(zz, numEdges);

  if (!vIdx || !edges){
    ZOLTAN_FREE(&vIdx);
    ZOLTAN_FREE(&edges);
    ierr = ZOLTAN_MEMERR;
  }

  /* Write out pins */

  pins = ZOLTAN_MALLOC_GID_ARRAY(zz, num_pins);
  if (!pins){
    ZOLTAN_FREE(&vIdx);
    ZOLTAN_FREE(&edges);
    ierr = ZOLTAN_MEMERR;
  }

  vgid = *vtx_GID;
  egid = *edg_GID;
  eIdx = col_ptr;

  vIdx[0] = 0;
  edg_cnt = 0 ;
  for (v=0; v < numVerts; v++){
    npins = eIdx[v+1] - eIdx[v];

    for (e=0; e < npins; e++){
      idx = Zoltan_Hash(egid, numGID, (unsigned int)ht_size);
      hn = hash_table[idx];

      while (hn){
        if (ZOLTAN_EQ_GID(zz, hn->egid, egid)){
            if (hn->loc < 0) {
                /* Never seen this edge before */
                hn->loc = -(hn->loc) ;
                vIdx[edg_cnt+1] = vIdx[edg_cnt] + hn->loc ;
                ZOLTAN_SET_GID(zz, edges + edg_cnt*numGID, hn->egid);
                hn->loc = vIdx[edg_cnt]; /* Use loc as the bucket pointer now */
                edg_cnt++ ;
            }

          ZOLTAN_SET_GID(zz, pins + (numGID * hn->loc), vgid);
          hn->loc++;

          break;
        }
        else{
          hn = hn->next;
        }
      }

      egid += numGID;
    }
    vgid += numGID;
  }

End:
  for (idx=0; idx<ht_size; idx++){
    hn = hash_table[idx];
    while (hn){
      tmp = hn;
      hn = hn->next;
      ZOLTAN_FREE(&tmp);
    }
  }
  ZOLTAN_FREE(&hash_table);

  *num_lists = numEdges;
  ZOLTAN_FREE(vtx_GID);
  *vtx_GID = pins;
  ZOLTAN_FREE(edg_GID);
  *edg_GID = edges;
  ZOLTAN_FREE(row_ptr);
  *row_ptr = vIdx;

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}
/*****************************************************************************/

int Zoltan_Graph_Queries( ZZ *zz, 
  int numVertex, ZOLTAN_ID_PTR vgid, ZOLTAN_ID_PTR vlid,        /* IN */
  int *tot_nbors, int **num_nbors,             /* OUT */
  ZOLTAN_ID_PTR *nbor_GIDs, int **nbor_Procs,  /* OUT */
  float **edgeWeights)                         /* OUT */
{
static char *yo = "Graph_Queries";
int ierr = ZOLTAN_OK;
int i;
int gid_size = zz->Num_GID;
int lid_size = zz->Num_LID;
int ew_dim = zz->Edge_Weight_Dim;
ZOLTAN_ID_PTR nbor_gids = NULL, gid_ptr = NULL;
int *nbor_procs = NULL, *proc_ptr = NULL;
float *gewgts = NULL, *wgt_ptr = NULL;
int sumNumEntries, temp;
int *numEdges = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  *tot_nbors = 0;
  *num_nbors = NULL;
  *nbor_GIDs = NULL;
  *nbor_Procs = NULL;
  *edgeWeights = NULL;

  ierr = Zoltan_Get_Num_Edges_Per_Obj(zz, numVertex, vgid, vlid,
       &numEdges,       /* number of neighbors for each vgid, length numVertex */
       &temp,           /* max neighbors of any vgid, we don't care */
       &sumNumEntries); /* total number of neighbors, sum of entries in numEdges */

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Zoltan_Get_Num_Edges_Per_Obj");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }

  if (sumNumEntries == 0){
    goto End;
  }
  
  nbor_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, sumNumEntries);
  nbor_procs = (int *)ZOLTAN_MALLOC(sizeof(int) * sumNumEntries);
  gewgts = (float *)ZOLTAN_MALLOC(sizeof(float) * ew_dim * sumNumEntries);

  if (!nbor_gids || !nbor_procs || (ew_dim && !gewgts)){
    Zoltan_Multifree(__FILE__, __LINE__, 3, &nbor_gids, &nbor_procs, &gewgts);
    ZOLTAN_FREE(&num_nbors);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }

  if (zz->Get_Edge_List_Multi){

    zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data,
      gid_size, lid_size, numVertex, vgid, vlid, numEdges,
      nbor_gids, nbor_procs, ew_dim, gewgts, &ierr);

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in edge list query function");
      Zoltan_Multifree(__FILE__, __LINE__, 3, &nbor_gids, &nbor_procs, &gewgts);
      ZOLTAN_FREE(&numEdges);
      goto End;
    }
  }
  else{

   gid_ptr = nbor_gids;
   proc_ptr = nbor_procs;
   wgt_ptr = gewgts;
 
    for (i=0; i<numVertex; i++){

      zz->Get_Edge_List(zz->Get_Edge_List_Data,
        gid_size, lid_size, 
        vgid + (i * gid_size), vlid + (i * lid_size), 
        gid_ptr, proc_ptr, ew_dim, wgt_ptr, &ierr);

      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in edge list query function");
        Zoltan_Multifree(__FILE__,__LINE__,3,&nbor_gids,&nbor_procs,&gewgts);
        ZOLTAN_FREE(&num_nbors);
        goto End;
      }

      gid_ptr += (numEdges[i] * gid_size);
      proc_ptr += numEdges[i];
      if (wgt_ptr) wgt_ptr += (ew_dim * numEdges[i]);
    }
  }

End:
  *tot_nbors = sumNumEntries;
  *num_nbors = numEdges;
  *nbor_GIDs = nbor_gids;
  *nbor_Procs = nbor_procs;
  *edgeWeights = gewgts;

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
/*****************************************************************************/

static int edge_weight_operation(ZZ *zz, float *dest, float *src, int ew_dim, int ew_op, int len)
{
  static char *yo = "edge_weight_operation";
  int k, w;
  int ierr = ZOLTAN_OK;

  for (k=0; k < len; k++){
    if (ew_op == PHG_FLAG_ERROR_EDGE_WEIGHTS){
      /* If the src weights are zero, it just means no weight was provided */
      for (w=0; w<ew_dim; w++){
        if ((src[w] != 0.0) && (src[w] != dest[w])){
          FATAL_ERROR("Different edge weights were supplied for the same edge");
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
    src += ew_dim;
    dest += ew_dim;
  }
End:
  return ierr;
}
/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
