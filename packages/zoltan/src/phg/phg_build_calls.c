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
#include "phg_hypergraph.h"
#include "parmetis_jostle.h"
    
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

#if VERBOSE_EDGE_INFO
static void show_edges(char *s, ZZ *zz, int num_lists, int num_pins, 
                int *edg_GID, int *row_ptr, int *vtx_GID);
#endif

static int ignore_some_edges(ZZ *, ZHG *, int, float, int, int *,
  ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, float *, ZOLTAN_ID_PTR, int);
static int remove_empty_edges(ZZ *zz, int *num_lists, int num_pins, 
   ZOLTAN_ID_PTR edg_GID, int *row_ptr, ZOLTAN_ID_PTR vtx_GID);
static int convert_to_CRS( ZZ *zz, int num_pins, int *col_ptr,
    int *num_lists, ZOLTAN_ID_PTR *vtx_GID,
    int **row_ptr, ZOLTAN_ID_PTR *edg_GID);
static int divideInterval(int rank, int myL, int myR,
                   int *l0, int *l1, int *r0, int *r1);
static int distribute_edges(ZZ *zz, int *num_lists, int *num_pins,
    ZOLTAN_ID_PTR *edg_GID, ZOLTAN_ID_PTR *edg_LID,
    int **row_ptr, ZOLTAN_ID_PTR *vtx_GID,
    int need_weights, int max_need_weights, char *out_need_list,
    int ew_table_size, void *htptr, float **edg_weight);
static int exchange(ZZ *zz, int proc, int *change, void *htptr,
   ZOLTAN_ID_PTR egid, ZOLTAN_ID_PTR inbuf, int inbufsize,
   int exchange_weights, int ew_num_edges, void *ewhtptr, float *edg_weights,
   char *out_need_list, int out_need_list_size,
   char *in_need_list, int in_need_list_size,
   float *out_weights, int out_weights_size,
   float *in_weights, int in_weights_size);
static int do_transfers(ZZ *zz, int proc, void *hn, char *mine, int numMatches,
  int *myMatch, ZOLTAN_ID_PTR egid, int *yourMatch, ZOLTAN_ID_PTR inbuf);
static float *find_weights(ZZ *zz,
              void *ewhtptr, int ew_num_edges, ZOLTAN_ID_PTR eid);
static char mpi_err_str[MPI_MAX_ERROR_STRING];
static int mpi_err_len;
/*****************************************************************************/

int Zoltan_HG_Hypergraph_Pin_Callbacks(
  ZZ *zz,
  ZHG *zhg,                /* Input:   Pointer to Zoltan's structure with
                                       GIDs, LIDs.
                              Output:  Removed edge fields of zhg are changed
                                       if dense edges are removed. */
  int gnVtx,               /* Input:   Global number of vertices in hgraph */
  float esize_threshold,   /* Input:   %age of gnVtx considered a dense edge */
  int return_removed,      /* Input:   flag indicating whether to return
                                       removed edge info */
  int ew_num_edges,        /* Input: num edges from edge weight queries */
  void *htptr,             /* Input: hash table for edge weights from query */
                           /*        we free hash table when done           */
  int htsize,              /* Input: size of hash table */
  int *nedges,             /* Output:  Number of hyperedges on this processor */
  ZOLTAN_ID_PTR *egids,    /* Output:  GIDs of hyperedges on this processor */
  ZOLTAN_ID_PTR *elids,    /* Output:  LIDs (if provided by application) */
  int **esizes,            /* Output:  # of vertices for each hyperedge on
                                       this processor */
  float **ewgts,           /* Output: Weights for each hyperedge on proc */
  int *npins,              /* Output:  # of pins on this processor =
                                       sum esizes[i], i = 0..nedges-1. */
  ZOLTAN_ID_PTR *pins)     /* Output:  vertex GIDs of pins */
{
/* Function to call the Zoltan Hypergraph pin callback functions.      */
/*                                                                     */
/* Processes return a set of pins, in either compressed row or         */
/* compressed column format.  We do some work to redistribute pins so  */
/* each process has complete rows (edges) and no two processes have    */
/* the same row.                                                       */
/*                                                                     */
/* If the application had defined the edge weight callbacks, the       */
/* supplied edges and weights are in the htptr table.  If they are not */
/* for the edges which our pins are part of, we we need to share       */
/* those while redistributing the pins.                                */

static char *yo = "Zoltan_HG_Hypergraph_Pin_Callbacks";
int ierr = ZOLTAN_OK;
int i, j, w;
ZOLTAN_ID_PTR vtx_GID=NULL, edg_GID=NULL, edg_LID, egptr, elptr;
float *edg_weight, *wptr;
int num_lists, num_pins;
int *row_ptr=NULL, *size_ptr;
struct _ewht{
  ZOLTAN_ID_PTR egid;
  ZOLTAN_ID_PTR elid;
  float *weights;
  struct _ewht *next;
} *ewNodes=NULL, *en;
struct _ewht **ewht = NULL;
int max_need_weights, need_weights;
char *need_list=NULL;
int dim = zz->Edge_Weight_Dim;
int lenGID = zz->Num_GID;
int lenLID = zz->Num_LID;

  ZOLTAN_TRACE_ENTER(zz, yo);

  ewht = (struct _ewht **)htptr;
  if (ewht){
    ewNodes = ewht[htsize]; /* so we can free it when done */
  }
  else{
    ewNodes = NULL;
  }

  ierr = Zoltan_Call_Hypergraph_Pin_Query(zz, &num_lists, &num_pins,
              &edg_GID, &row_ptr, &vtx_GID);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    Zoltan_Multifree(__FILE__, __LINE__, 3, &vtx_GID, &row_ptr, &edg_GID);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
       "Error calling hypergraph compressed pin storage query");
    goto End;
  }

  /*
   * Remove edges with no pins.  This would happen later in 
   * ignore_some_edges, but doing it now makes distribute_edges() faster.
   */

  ierr = 
    remove_empty_edges(zz, &num_lists, num_pins, edg_GID, row_ptr, vtx_GID);

  /*
   * If we have edge weights for any of our edges, apply them now.
   * Also, save the edge local ID.  The only way we get edge local
   * IDs is if the application returned them in the edge weight
   * query function.
   */
  edg_weight = NULL;
  need_weights = max_need_weights = 0;
  edg_LID = ZOLTAN_MALLOC_LID_ARRAY(zz, num_lists);

  if (num_lists && !edg_LID) MEMORY_ERROR;

  if (num_lists){
    memset(edg_LID, 0, num_lists * lenLID * sizeof(int));
  }

  if (dim > 0){
    edg_weight = (float *)ZOLTAN_MALLOC(dim * num_lists * sizeof(float));

    if (num_lists && !edg_weight){
      ZOLTAN_FREE(&edg_LID);
      MEMORY_ERROR;
    }

    for (i=0; i< dim*num_lists; i++){
      edg_weight[i] = 1.0;             /* default edge weight */
    }

    MPI_Allreduce(&ew_num_edges, &j, 1, MPI_INT, MPI_MAX, zz->Communicator);

    if (j > 0){
      /* If at least one process supplied weights in the query function
       * we need to process them.  Otherwise we go with the unit default.
       */
      need_weights = num_lists;
      need_list = (char *)ZOLTAN_MALLOC(num_lists * sizeof(char));

      if (num_lists && !need_list){
        ZOLTAN_FREE(&edg_LID);
        MEMORY_ERROR;
      }
      if (num_lists){
        memset((void *)need_list, 1, num_lists);
      }

      if (ew_num_edges && num_lists){
        egptr = edg_GID;
        elptr = edg_LID;
        wptr = edg_weight;
        for (i=0; i < num_lists; i++){
          j = Zoltan_Hash(egptr, lenGID, (unsigned int)htsize);
          en = ewht[j];
          while (en){
            if (ZOLTAN_EQ_GID(zz, egptr, en->egid)){
              ZOLTAN_SET_LID(zz, elptr, en->elid);
              for (w=0; w<dim; w++){
                wptr[w] = en->weights[w];
              }
              need_weights--;
              need_list[i] = 0;
              break;
            }
            en = en->next;
          }

          egptr += lenGID;
          elptr += lenLID;
          wptr += dim;
        }
      }

      MPI_Allreduce(&need_weights, &max_need_weights, 1,
                    MPI_INT, MPI_MAX, zz->Communicator);

      if (max_need_weights == 0){
        /* Each process included in it's edge weight query function
         * all edges associated with it's pins.  We are done processing
         * edge weights.
         */
        need_weights = 0;
        ZOLTAN_FREE(&need_list);
      }
    }
  }

  /*
   * Each processor has in general a set of incomplete rows.
   * Redistribute the pins so each processor has a set of complete
   * rows, and no edge (row) is on more than one processor.
   * This is necessary so we can enumerate all edges in the matrix.
   */

#if VERBOSE_EDGE_INFO 
  for (i=0; i<zz->Num_Proc; i++){
    if (i == zz->Proc){
      show_edges("BEFORE", zz, num_lists, num_pins, edg_GID, row_ptr, vtx_GID);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  ierr =
    distribute_edges(zz, &num_lists, &num_pins,
                     &edg_GID, &edg_LID, &row_ptr, &vtx_GID,
                     need_weights, max_need_weights, need_list,
                     htsize, htptr, &edg_weight);

#if VERBOSE_EDGE_INFO
  for (i=0; i<zz->Num_Proc; i++){
    if (i == zz->Proc){
      show_edges("AFTER", zz, num_lists, num_pins, edg_GID, row_ptr, vtx_GID);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  fflush(stdout);
#endif

  ZOLTAN_FREE(&need_list);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    Zoltan_Multifree(__FILE__, __LINE__, 5,
       &vtx_GID, &row_ptr, &edg_GID, &edg_LID, &edg_weight);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error redistributing edges");
    goto End;
  }

  *nedges = num_lists;
  *egids = edg_GID;
  *elids = edg_LID;
  *ewgts = edg_weight;
  size_ptr = *esizes = NULL;

  if (num_lists > 0){
    size_ptr = *esizes = (int *) ZOLTAN_MALLOC(num_lists * sizeof(int));
    for (i=0; i<num_lists-1; i++){
      size_ptr[i] = row_ptr[i+1] - row_ptr[i];
    }
    size_ptr[num_lists-1] = num_pins - row_ptr[num_lists - 1];

    ZOLTAN_FREE(&row_ptr);
  
    /* Remove dense edges from input list */
    ierr = ignore_some_edges(zz, zhg, gnVtx, esize_threshold, return_removed,
                   nedges, edg_GID, edg_LID, *esizes, edg_weight, vtx_GID, 0);
  
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from ignore_some_edges.");
      goto End;
    }
  } else {
    ZOLTAN_FREE(egids);
    ZOLTAN_FREE(elids);
    ZOLTAN_FREE(esizes);
    ZOLTAN_FREE(ewgts);
    ZOLTAN_FREE(&vtx_GID);
  }

  num_pins = 0;
  for (i=0; i<*nedges; i++){
    num_pins += size_ptr[i];
  }

  *npins = num_pins;
  *pins = vtx_GID;

End:

  ZOLTAN_FREE(&ewNodes);
  ZOLTAN_FREE(&ewht);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
/*****************************************************************************/

static int distribute_edges(ZZ *zz, int *num_lists, int *num_pins,
    ZOLTAN_ID_PTR *edg_GID, ZOLTAN_ID_PTR *edg_LID,
    int **row_ptr, ZOLTAN_ID_PTR *vtx_GID,
    int need_weights, int max_need_weights, char *out_need_list,
    int ew_table_size, void *htptr, float **edg_weights)
{
static char *yo = "distribute_edges";
int nprocs = zz->Num_Proc;
int rank = zz->Proc;
int lenGID = zz->Num_GID;
int lenLID = zz->Num_LID;
int ierr = ZOLTAN_OK;
int ew_dim = zz->Edge_Weight_Dim;
int i, j, maxEdges, nedges, npins, update_lists, pin_change;
int *rptr=NULL;
ZOLTAN_ID_PTR egid=NULL, eptr=NULL, elidptr=NULL, vptr=NULL, inbuf=NULL;
struct _egidNode {
  ZOLTAN_ID_PTR egidBuf;
  int idx;
  ZOLTAN_ID_PTR pins;
  char pinsAreCopies;
  struct _egidNode *next;
} *egidNode=NULL, *en=NULL;
struct _egidNode **ht=NULL;  /* hash table for my pin edges */
int l0, l1, lLen, myL, left;
int r0, r1, rLen, myR;
int relativeRank, inbufsize, have_pins, esize, v, e;
int in_need_list_size=0, in_weights_size=0, out_weights_size=0;
char *in_need_list=NULL;
float *in_weights=NULL, *out_weights=NULL, *wptr=NULL;

/* If there is an edge such that more than one process has pins for  */
/* that edge, assign the edge to only one process.  Transfer all     */
/* pins to that process.                                             */ 
/*                                                                   */
/* Also obtain edg_weights for my edges from processes which         */
/* supplied them in the edge weight query functions.                 */

  ZOLTAN_TRACE_ENTER(zz, yo);

  nedges = *num_lists;
  npins = *num_pins;
  MPI_Allreduce(&nedges, &maxEdges, 1, MPI_INT, MPI_MAX, zz->Communicator);

  if (maxEdges <= 0){
    ierr = ZOLTAN_WARN;
    goto End;
  }

  have_pins = ((nedges > 0) && (npins > 0));
  inbufsize = 2 * maxEdges + 1;
  inbuf = ZOLTAN_MALLOC_GID_ARRAY(zz, inbufsize);

  if (!inbuf) MEMORY_ERROR;

  /* Create a buffer to send to other processes listing:
   * --how many edges I have
   * --each edge GID followed by the number of pins (0 or greater) I have
   *
   * Create a structure for each edge listing the pins in that edge.
   * Create a hash table based on edge GID which accesses
   * these structures.
   */

  egid = ZOLTAN_MALLOC_GID_ARRAY(zz, 2*nedges+1);
  if (!egid){
    ZOLTAN_FREE(&inbuf);
    MEMORY_ERROR;
  }
  egid[0] = nedges;

  if (have_pins){
    eptr = egid + lenGID;
    ht = (struct _egidNode **)ZOLTAN_CALLOC(nedges+1,sizeof(struct _egidNode *));
    egidNode = (struct _egidNode *)ZOLTAN_MALLOC(nedges*sizeof(struct _egidNode));

    if (!ht || !egidNode){
      Zoltan_Multifree(__FILE__, __LINE__, 2, &inbuf, &ht);
      MEMORY_ERROR;
    }

    ht[nedges] = egidNode;
    vptr = *vtx_GID;

    for (i=0; i<nedges; i++){
      ZOLTAN_SET_GID(zz, eptr, *edg_GID + (i*lenGID));
      esize = (((i == nedges-1) ? npins : (*row_ptr)[i+1]) - (*row_ptr)[i]);

      egidNode[i].egidBuf = eptr;
      egidNode[i].idx     = i;
      egidNode[i].pins    = vptr;
      egidNode[i].pinsAreCopies = 0;

      vptr += (esize * lenGID);

      j = Zoltan_Hash(eptr, lenGID, (unsigned int)nedges);

      egidNode[i].next = ht[j];
      ht[j] = egidNode + i;

      eptr += lenGID;
      *eptr = esize;
      eptr += lenGID;
    }
  }
  else{
    ht = NULL;
    egidNode = NULL;
  }

  if (max_need_weights > 0){ /* buffers to request & exchange edge weights */
    in_need_list_size = maxEdges;
    out_weights_size = (max_need_weights * (1 + zz->Edge_Weight_Dim)) + 1;
    in_weights_size = (need_weights * (1 + zz->Edge_Weight_Dim)) + 1;

    in_need_list = (char *)ZOLTAN_MALLOC(in_need_list_size * sizeof(char));
    in_weights = (float *)ZOLTAN_MALLOC(in_weights_size * sizeof(float));
    out_weights = (float *)ZOLTAN_MALLOC(out_weights_size * sizeof(float));

    if (!in_need_list || !in_weights || !out_weights){
      Zoltan_Multifree(__FILE__, __LINE__, 5, &inbuf, &ht,
                       &in_need_list, &in_weights, &out_weights);
      MEMORY_ERROR;
    }
  }

  /*
   **Problem: Edges may be spread across processes.  (More
   * than one process owns pins for the same edge, although
   * no two processes own the same pin.)
   **Goal: Assign each edge to only one process.  Each process
   * then acquires all pins for each of it's edges.
   **Other problem: We assume the combined list of every
   * process' edge list is in general too large to fit in the
   * memory allocated to a single process. So we don't do a
   * MPI_Allgatherv of edge lists and have the processes review
   * all the lists at once and assign edges to processes.
   *
   * In nprocs-1 exchanges, each process exchanges edge ID lists
   * with every other process.  Following one exchange, the
   * two processes involved decide which process should own
   * each edge they have in common.  The owning process gets
   * all pins for the edge it owns.  This decision and exchange
   * must be completed before beginning the exchange with the next
   * process.  And a process is never assigned an edge that it
   * did not have at the beginning of the exchange.
   *
   * At the end of nprocs-1 exchanges, each edge is owned by
   * only one process, and that process has all pins for the edge,
   */

  /*
   * Divide processes roughly in half.  All processes in the
   * left half will perform an exchange with all processes in
   * the right half.  Recurse, dividing left and right roughly
   * in half.
   */

  myL = 0;
  myR = nprocs-1;
  update_lists = 0;
  ierr = ZOLTAN_OK;

  while ((myR - myL > 0) && (ierr != ZOLTAN_FATAL)){

    left = divideInterval(rank, myL, myR, &l0, &l1, &r0, &r1);

    if (left){                /* I'm in the left half */
      rLen = r1 - r0 + 1;
      relativeRank = rank - l0;

      for (i=0; i<rLen; i++){
        j = r0 + ((relativeRank + i) % rLen);

        ierr = exchange(zz, j, &pin_change, (void *)ht, egid, inbuf, inbufsize,
               max_need_weights, ew_table_size, (void *)htptr, *edg_weights,
               out_need_list, nedges,
               in_need_list, in_need_list_size,
               out_weights, out_weights_size,
               in_weights, in_weights_size);

        if (ierr != ZOLTAN_OK){
          break; 
        }

        if (!update_lists && pin_change){
          update_lists = 1;
        }
      }
      myL = l0;
      myR = l1;
    }
    else{                    /* I'm in the right half */
      lLen = l1 - l0 + 1;
      relativeRank = rank - r0;

      for (i=0; i<lLen; i++){
        j = l0 + ((relativeRank + lLen - i) % lLen);

        ierr = exchange(zz, j, &pin_change, (void *)ht, egid, inbuf, inbufsize,
               max_need_weights, ew_table_size, (void *)htptr, *edg_weights,
               out_need_list, nedges,
               in_need_list, in_need_list_size,
               out_weights, out_weights_size,
               in_weights, in_weights_size);

        if (ierr != ZOLTAN_OK){
          break; 
        }

        if (!update_lists && pin_change){
          update_lists = 1;
        }
      }
      myL = r0;
      myR = r1;
    }
  }

  ZOLTAN_FREE(&inbuf);
  ZOLTAN_FREE(&ht);
  ZOLTAN_FREE(&in_need_list);
  ZOLTAN_FREE(&in_weights);
  ZOLTAN_FREE(&out_weights);

  if (update_lists && (ierr != ZOLTAN_FATAL)){ /* Must rewrite my edge lists */

    nedges = npins = 0;

    for (i=0; i< egid[0]; i++) {   /* former number of edges */
      en = egidNode + i;
      if (en->egidBuf[lenGID] > 0){
        npins += en->egidBuf[lenGID];
        nedges++;
      }
    }

    if (npins > 0){
      *num_lists = nedges;
      *num_pins = npins;

      rptr = (int *)ZOLTAN_MALLOC(sizeof(int) * nedges);
      eptr = ZOLTAN_MALLOC_GID_ARRAY(zz, nedges);
      elidptr = ZOLTAN_MALLOC_LID_ARRAY(zz, nedges);
      vptr = ZOLTAN_MALLOC_GID_ARRAY(zz, npins);
      wptr = (float *)ZOLTAN_MALLOC(sizeof(float) * nedges * ew_dim);

      if (!eptr || !rptr || !elidptr || !vptr || (ew_dim && !wptr)){
        Zoltan_Multifree(__FILE__, __LINE__, 5,
          &eptr, &rptr, &elidptr, &vptr, &wptr);
        MEMORY_ERROR;
      }

      for (i=0, v=0, e=0; i< egid[0]; i++) {
        en = egidNode + i;
        npins = en->egidBuf[lenGID];

        if (npins > 0){
          ZOLTAN_SET_GID(zz, eptr + (e*lenGID), en->egidBuf);
          ZOLTAN_SET_LID(zz, elidptr + (e*lenLID), *edg_LID + (i*lenLID));
          ZOLTAN_COPY_GID_ARRAY(vptr + (v*lenGID), en->pins, zz, npins);
          for (j=0; j<ew_dim; j++){
            wptr[e*ew_dim+j] = (*edg_weights)[i*ew_dim+j];
          }
          rptr[e] = v;
          e++;
          v += npins;
        }
      }

      ZOLTAN_FREE(edg_GID);
      *edg_GID = eptr;
      ZOLTAN_FREE(edg_LID);
      *edg_LID = elidptr;
      ZOLTAN_FREE(vtx_GID);
      *vtx_GID = vptr;
      ZOLTAN_FREE(row_ptr);
      *row_ptr = rptr;
      ZOLTAN_FREE(edg_weights);
      *edg_weights = wptr;
    }
    else{
      *num_lists = 0;
      *num_pins = 0;
      ZOLTAN_FREE(edg_GID);
      ZOLTAN_FREE(vtx_GID);
      ZOLTAN_FREE(row_ptr);
      ZOLTAN_FREE(edg_weights);
    }
  }

End:

  if (egid) {
    for (i=0; i< egid[0]; i++) {
      if (egidNode[i].pinsAreCopies){
        ZOLTAN_FREE(&egidNode[i].pins);
      }
    }
  }
  ZOLTAN_FREE(&egidNode);
  ZOLTAN_FREE(&egid);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
static int divideInterval(int rank, int myL, int myR,
                   int *l0, int *l1, int *r0, int *r1)
{
int inLeftHalf= 0;
int center;

/* Divide interval of process ranks in half, if odd number of      */
/* processes, one "half" will have one more process than the other.*/

  center = (myL + myR) / 2;

  *l0 = myL;
  *l1 = center;
  *r0 = center+1;
  *r1 = myR;

  if (rank <= *l1){
    inLeftHalf = 1;
  }

  return inLeftHalf;
}
static int exchange(ZZ *zz, int proc, int *change, void *htptr,
     ZOLTAN_ID_PTR egid, ZOLTAN_ID_PTR inbuf, int inbufsize,
     int exchange_weights, int htsize, void *ewhtptr, float *edg_weights,
     char *out_need_list, int out_need_list_size,
     char *in_need_list, int in_need_list_size,
     float *out_weights, int out_weights_size,
     float *in_weights, int in_weights_size)
{
static char *yo = "exchange";
struct _egidNode {
  ZOLTAN_ID_PTR egidBuf;
  int idx;
  ZOLTAN_ID_PTR pins;
  char pinsAreCopies;
  struct _egidNode *next;
} *en;                    /* hash table for edges I have pins for */
struct _egidNode **ht;
int tag=0xff, need_tag=0xf1, weight_tag=0xf2;
int lenGID = zz->Num_GID;
int numMyEdges, numYourEdges, numAEdges, numBEdges, numMyPins, numYourPins;
int *myMatch=NULL, *yourMatch=NULL;
int *matchA, *matchB, nMatches;
int rc, i, j, w, nweights_found=0;
int favorA, favorB, idxA, idxB;
MPI_Request req, need_list_req, weight_req;
MPI_Status status;
ZOLTAN_ID_PTR eptr, eid, gidsA, gidsB;
float *wptr=NULL, *wgts;
char *I_Win = NULL; 
int *gidsAMatchIdx=NULL;
int dim = zz->Edge_Weight_Dim;
int ierr=ZOLTAN_OK;

  /* Two processes exchange lists of edges IDs and sizes.  If  */
  /* a match of non-empty edges is found, the edge is assigned */
  /* to only one process, and the pins are transferred.        */
  /*                                                           */
  /* Also during this exchange, a process supplies edge weights*/
  /* it may have that are required by the other process.       */
  /*                                                           */
  /* If this exchange results in a change of pins for the      */
  /* process, the change flag is set, indicating that edge/pin */
  /* lists need to be re-written when all exchanges are        */
  /* completed.                                                */

  *change = 0;

  ht = (struct _egidNode **)htptr;

  /* post receive for other process' edge GIDs, and possibly edge weights
      needed array and edge weights response */

  rc = MPI_Irecv(inbuf, inbufsize * lenGID, MPI_INT, proc,
            tag, zz->Communicator, &req);

  CHECK_FOR_MPI_ERROR(rc);

  if (exchange_weights){
    rc = MPI_Irecv(in_need_list, in_need_list_size, MPI_CHAR, proc, need_tag,
                 zz->Communicator, &need_list_req);
    CHECK_FOR_MPI_ERROR(rc);

    rc = MPI_Irecv(in_weights, in_weights_size, MPI_FLOAT, proc, weight_tag,
                 zz->Communicator, &weight_req);
    CHECK_FOR_MPI_ERROR(rc);
  }

  /* send my edge GIDs to that process, along with possibly my edge
       weights needed flags */

  numMyEdges = egid[0];

  rc = MPI_Send(egid, (2*numMyEdges + 1) * lenGID, MPI_INT, proc,
            tag, zz->Communicator);
  CHECK_FOR_MPI_ERROR(rc);

  if (exchange_weights){
    rc = MPI_Send(out_need_list, out_need_list_size, MPI_CHAR, proc,
               need_tag, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc);
  }


  /* await GIDs, then process to see if we have any in common,
       and possibly if we have any edge weights they need */

  rc = MPI_Wait(&req, &status);
  CHECK_FOR_MPI_ERROR(rc);

  if (exchange_weights){
    rc = MPI_Wait(&need_list_req, &status);
    CHECK_FOR_MPI_ERROR(rc);
    wptr = out_weights + 1;
    nweights_found = 0;
  }

  eptr = inbuf;
  numYourEdges = eptr[0]; 

  if (numYourEdges > (inbufsize-1)/2){
    /* Sanity check: number of edges exceeds maximum number of edges */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid gid message from remote process");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  j = ((numYourEdges < numMyEdges) ? numYourEdges : numMyEdges);

  if ((j <= 0) && !exchange_weights){
    return ZOLTAN_OK; 
  }

  myMatch   = (int *)ZOLTAN_MALLOC(j * sizeof(int));
  yourMatch = (int *)ZOLTAN_MALLOC(j * sizeof(int));

  if (j && (!myMatch || !yourMatch)){
    Zoltan_Multifree(__FILE__, __LINE__, 2, &myMatch, &yourMatch);
    MEMORY_ERROR
  }

  eptr += lenGID;
  nMatches = 0;

  for (i=0; i<numYourEdges; i++){

    eid         = eptr;

    if (exchange_weights && in_need_list[i]){
      wgts = find_weights(zz, ewhtptr, htsize, eid);
      if (wgts){
        *wptr++ = (float)i;
        for (w=0; w<dim; w++){
          *wptr++ = wgts[w];
        }
        nweights_found++;
      }
    }
    numYourPins = eptr[lenGID];
    numMyPins   = 0;

    if (numMyEdges && (numYourPins > 0)){
      j = Zoltan_Hash(eid, lenGID, (unsigned int)numMyEdges);

      en = ht[j];

      while (en){
        if (ZOLTAN_EQ_GID(zz, eid, en->egidBuf)){
          numMyPins = en->egidBuf[lenGID];
          break;
        }
        en = en->next;
      }

      if (numMyPins > 0){    /* save to evaluate later */
        myMatch[nMatches]   = en->idx;
        yourMatch[nMatches] = i;
        nMatches++;
      }
    }

    eptr += (2*lenGID);
  }

  /* If we are processing edge weights, send a message back with
       any of their needed edge weights that we have, await their
       weights and add them to my list of edge weights */

  if (exchange_weights){
    out_weights[0] = nweights_found;

    rc = MPI_Send(out_weights, (nweights_found * (dim+1)) + 1, MPI_FLOAT, proc,
             weight_tag, zz->Communicator);
    CHECK_FOR_MPI_ERROR(rc);

    rc = MPI_Wait(&weight_req, &status);
    CHECK_FOR_MPI_ERROR(rc);

    if (in_weights[0] > 0){
      wgts = in_weights + 1;

      for (i=0; i<in_weights[0]; i++){
        j = (int)*wgts++;
        wptr = edg_weights + (j * dim);
        for (w=0; w<dim; w++){
          wptr[w] = *wgts++;
        }
        out_need_list[j] = 0;
      }
    }
  }

  /* If we share edges in common, proceed to decide who gets the edge. */

  if (nMatches == 0){
    ierr = ZOLTAN_OK;
    goto End;
  }

  *change = 1;

  /* Determine which process gets each shared edge.  Be careful to
   * do the calculation in exactly the same way on both processes.
   * Call the lower rank proc "A", the higher rank proc is "B".
   */

  if (proc < zz->Proc){
    matchA    = yourMatch;
    gidsA     = inbuf;
    numAEdges = numYourEdges;
    matchB    = myMatch;
    gidsB     = egid;
    numBEdges = numMyEdges;
  }
  else{
    matchB    = yourMatch;
    gidsB     = inbuf;
    numBEdges = numYourEdges;
    matchA    = myMatch;
    gidsA     = egid;
    numAEdges = numMyEdges;
  }

  I_Win= (char *)ZOLTAN_MALLOC(sizeof(char) * nMatches);
  if (!I_Win) MEMORY_ERROR

  favorA = favorB = 0;

  /* This array allows us to evaluate matches in the same order
   * on both processes.
   */
  gidsAMatchIdx = (int *)ZOLTAN_MALLOC(sizeof(int) * numAEdges);
  if (!gidsAMatchIdx) MEMORY_ERROR

  for (i=0; i<numAEdges; i++){
    gidsAMatchIdx[i] = -1;
  }
  for (i=0; i<nMatches; i++){
    gidsAMatchIdx[matchA[i]] = i;
  }

  w = numAEdges;

  for (j=0; j<w; j++){

    if (gidsAMatchIdx[j] < 0) continue;

    i = gidsAMatchIdx[j];

    idxA = (2*matchA[i] + 2) * lenGID; /* location of A's num pins */
    idxB = (2*matchB[i] + 2) * lenGID; /* location of B's num pins */


    /* Give edge to the process having the most pins.  In case
     * of a tie, give it to the process having the least edges.
     * In case of a tie again, give it to the process of lowest rank.
     *
     * If a large imbalance in number of edges occurs, correct this
     * by favoring one process until the imbalance is corrected.
     */

    if (favorA){
      I_Win[i] = (matchA == myMatch);
      numBEdges--;
    }
    else if (favorB){
      I_Win[i] = (matchB == myMatch);
      numAEdges--;
    }
    else if (gidsA[idxA] < gidsB[idxB]){
      I_Win[i] = (matchB == myMatch);
      numAEdges--;
    }
    else if (gidsA[idxA] > gidsB[idxB]){
      I_Win[i] = (matchA == myMatch);
      numBEdges--;
    }
    else if (numAEdges < numBEdges){
      I_Win[i] = (matchA == myMatch);
      numBEdges--;
    }
    else if (numAEdges > numBEdges){
      I_Win[i] = (matchB == myMatch);
      numAEdges--;
    }
    else{
      I_Win[i] = (matchA == myMatch);
      numBEdges--;
    }

    favorA = numAEdges * 2 < numBEdges;
    favorB = numBEdges * 2 < numAEdges;
  }

  /* exchange pins, update my edge and pin information */

  ierr = do_transfers(zz, proc, (void *)ht, I_Win,
            nMatches, myMatch, egid, yourMatch, inbuf);

End:

  ZOLTAN_FREE(&gidsAMatchIdx);
  ZOLTAN_FREE(&I_Win);
  ZOLTAN_FREE(&myMatch);
  ZOLTAN_FREE(&yourMatch);

  return ierr;
}
static float *find_weights(ZZ *zz,
              void *ewhtptr, int htsize, ZOLTAN_ID_PTR eid)
{
struct _ewht{
  ZOLTAN_ID_PTR egid;
  ZOLTAN_ID_PTR elid;
  float *weights;
  struct _ewht *next;
} *en;
struct _ewht **ewht = NULL;
float *weights = NULL;
int idx;

  /* Search the edge weights returned by the edge weight query */
  /* function to see if I have the edge weights for this edge. */

  if (htsize < 1) return NULL;

  ewht = (struct _ewht **)ewhtptr;

  idx = Zoltan_Hash(eid, zz->Num_GID, (unsigned int)htsize);

  en = ewht[idx];

  while (en){
    if (ZOLTAN_EQ_GID(zz, en->egid, eid)){
      weights = en->weights;
      break;
    }
    en = en->next;
  }

  return weights;
}
static int do_transfers(ZZ *zz, int proc, void *hn, char *mine, int numMatches,
  int *myMatch, ZOLTAN_ID_PTR egid, int *yourMatch, ZOLTAN_ID_PTR inbuf)
{
static char *yo = "do_transfers";
struct _egidNode {
  ZOLTAN_ID_PTR egidBuf;
  int idx;
  ZOLTAN_ID_PTR pins;
  char pinsAreCopies;
  struct _egidNode *next;
} *egidNodes, *en;
struct _egidNode **egidTable;
ZOLTAN_ID_PTR eptr;
int lenGID = zz->Num_GID;
int recvSize, sendSize, i, j, nrecvs, idx, npins, nedges, nMyPins, rc;
int tagPins = 0xf0f0;
ZOLTAN_ID_PTR recvBuf=NULL, sendBuf=NULL, pinBuf=NULL;
MPI_Request reqPins;
MPI_Status status;
int ierr = ZOLTAN_OK;

  /* Two processes have pins for the same edge.  Transfer all the */
  /* pins to one process.                                         */

  egidTable = (struct _egidNode **)hn;
  nedges = egid[0];
  egidNodes = egidTable[nedges];

  /* Calculate message sizes */

  recvSize = 0;
  sendSize = 0;
  nrecvs = 0;

  for (i=0; i<numMatches; i++){
    if (mine[i]){
      j = (2 * yourMatch[i] + 2) * lenGID;
      recvSize += (2 + inbuf[j]);
      nrecvs++;
    }
    else{
      j = (2 * myMatch[i] + 2) * lenGID;
      sendSize += (2 + egid[j]);
    }
  }

  /* Post receive for incoming edge pins */

  if (recvSize > 0){
    recvBuf = (unsigned int *)ZOLTAN_MALLOC_GID_ARRAY(zz, recvSize);
    if (!recvBuf) MEMORY_ERROR

    rc = MPI_Irecv(recvBuf, recvSize * lenGID, MPI_INT, proc,
              tagPins, zz->Communicator, &reqPins);

    CHECK_FOR_MPI_ERROR(rc)
  }

  /*
   * Create a ZOLTAN_ID_TYPE array containing pins that are
   * going to the other process.  For each edge send:
   *
   *  your edge ID index
   *  number of pins
   *  vertex GID for each pin I have
   */

  if (sendSize > 0){
    sendBuf = (unsigned int *)ZOLTAN_MALLOC_GID_ARRAY(zz, sendSize);
    if (!sendBuf) MEMORY_ERROR

    eptr = sendBuf;

    for (i=0; i<numMatches; i++){

      if (mine[i] == 0){       /* I give this edge to other proc */

        j = (2 * myMatch[i] + 2) * lenGID;
        npins = egid[j];

        egid[j] = 0;             /* now I have no pins in edge */

        *eptr = yourMatch[i];                   /* your edge ID index */
        eptr += lenGID;

        *eptr = npins;                    /* number of pins to follow */
        eptr += lenGID;

        en = egidNodes + myMatch[i];

        ZOLTAN_COPY_GID_ARRAY(eptr, en->pins, zz, npins); /* the pins */

        eptr += (npins * lenGID);

        if (en->pinsAreCopies){
          ZOLTAN_FREE(&en->pins);
          en->pinsAreCopies = 0;
        }
      }
    }

    rc = MPI_Send(sendBuf, sendSize*lenGID, MPI_INT, proc, tagPins,
             zz->Communicator);

    CHECK_FOR_MPI_ERROR(rc)

    ZOLTAN_FREE(&sendBuf);
  }

  /* Add pins received from other process to my list.  */

  if (recvSize > 0){
    rc = MPI_Wait(&reqPins, &status);
    CHECK_FOR_MPI_ERROR(rc)

    eptr = recvBuf;

    for (i=0; i<nrecvs; i++){
      idx = *eptr;
      eptr += lenGID;
      npins = *eptr;
      eptr += lenGID;

      if (npins > 0){
        en = egidNodes + idx;
        nMyPins = en->egidBuf[lenGID];
        pinBuf = (unsigned int *)ZOLTAN_MALLOC_GID_ARRAY(zz, nMyPins + npins);
        if (nMyPins > 0){
          ZOLTAN_COPY_GID_ARRAY(pinBuf, en->pins, zz, nMyPins);
        }
        ZOLTAN_COPY_GID_ARRAY(pinBuf+(nMyPins*lenGID), eptr, zz, npins);
        if (en->pinsAreCopies){
          ZOLTAN_FREE(&en->pins);
        }
        else{
          en->pinsAreCopies = 1;
        }
        en->pins = pinBuf;
        en->egidBuf[lenGID] = nMyPins + npins;

        eptr += (npins * lenGID);
      }
    }
  }

End:

  ZOLTAN_FREE(&recvBuf);

  return ierr;
}
/*****************************************************************************/

int Zoltan_Call_Hypergraph_Pin_Query(ZZ *zz, 
   int *num_lists,         /* output: number of edges */
   int *num_pins,          /* output: total number of pins in edges */
   ZOLTAN_ID_PTR *edg_GID, /* output: list of edge global IDs */
   int **row_ptr,          /* output: loc in vtx_GID for start of each edge */
                           /*         plus num_pins in last element         */
   ZOLTAN_ID_PTR *vtx_GID) /* output: vertex global ID for each pin */
{
static char *yo = "Zoltan_Call_Hypergraph_Pin_Query";
int ierr = ZOLTAN_OK;
int nl, np, format, have_pins, row_storage;
ZOLTAN_ID_PTR vid, eid;
int *rptr, *cptr;

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

  if ((format != ZOLTAN_COMPRESSED_ROWS)&&(format != ZOLTAN_COMPRESSED_COLS)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Invalid compression format returned in Get_HG_Size_CS");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }
  ZOLTAN_TRACE_DETAIL(zz, yo, "done with Get_HG_Size_CS");

  have_pins = ((nl > 0) && (np > 0));
  row_storage = (format == ZOLTAN_COMPRESSED_ROWS);

  /* Get the hypergraph pins in compressed storage format */

  if (have_pins){
    if (!row_storage){    /* compressed column storage */

      vid = ZOLTAN_MALLOC_GID_ARRAY(zz, nl);
      cptr = (int *)ZOLTAN_MALLOC(nl * sizeof(int));
      eid = ZOLTAN_MALLOC_GID_ARRAY(zz, np);

      if (!vid|| !cptr || !eid){
        Zoltan_Multifree(__FILE__, __LINE__, 3, &vid, &cptr, &eid);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "memory allocation");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_FATAL;
      }
      zz->Get_HG_CS(zz->Get_HG_CS_Data, zz->Num_GID,
               nl, np, format, vid, cptr, eid, &ierr);

      ZOLTAN_TRACE_DETAIL(zz, yo, "done with Get_HG_CS");

      if ((ierr == ZOLTAN_OK) || (ierr == ZOLTAN_WARN)){
        ierr = convert_to_CRS(zz,
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
static int convert_to_CRS(
    ZZ *zz, int num_pins, int *col_ptr,   /* input */
    int *num_lists,                       /* rest are input/output */
    ZOLTAN_ID_PTR *vtx_GID,
    int **row_ptr,
    ZOLTAN_ID_PTR *edg_GID)
{
static char *yo = "convert_to_CRS";
int numVerts = *num_lists;
int numEdges, ierr, ht_size;
ZOLTAN_ID_PTR egid, vgid;
int v, e, idx, found, npins;
struct _hash_node {
  ZOLTAN_ID_PTR egid;
  int numVerts;
  int firstVert;
  int nextVert;
  struct _hash_node *next;
} *hn=NULL, *tmp;
struct _hash_node **hash_table=NULL;
ZOLTAN_ID_PTR edges=NULL, pins=NULL;
int *eIdx=NULL, *vIdx=NULL;
int numGID = zz->Num_GID;

  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = ZOLTAN_OK;

  if (num_pins == 0){
    return ierr;
  }

  /*
   * Convert from CCS to CRS (compressed columns to compressed rows)
   * We have the lists of edges for each vertex.  Create lists of
   * vertices for each edge.
   */

  ht_size = (int)sqrt((double)num_pins);

  if (ht_size < 10) ht_size = num_pins;

  hash_table =
    (struct _hash_node **)ZOLTAN_CALLOC(ht_size, sizeof(struct _hash_node *));

  if (!hash_table){
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
         hn->numVerts++;
         found = 1;
         break;
       }
       else{
         hn = hn->next;
       }
     }
     if (!found){
       hn = (struct _hash_node *)ZOLTAN_MALLOC(sizeof(struct _hash_node));
       if (!hn){
         ierr = ZOLTAN_MEMERR;
         goto End;
       }
       hn->egid = egid;
       hn->numVerts = 1;
       hn->next = hash_table[idx];
       hash_table[idx] = hn;
       numEdges++;
     }
     egid += numGID;
  }

  /* Create array of indices into the start of each edge's pins,
   * and the list of unique edge IDs.                          
   */

  vIdx = (int *)ZOLTAN_MALLOC((numEdges+1) * sizeof(int));
  edges = ZOLTAN_MALLOC_GID_ARRAY(zz, numEdges);

  if (!vIdx || !edges){
    ZOLTAN_FREE(&vIdx);
    ZOLTAN_FREE(&edges);
    ierr = ZOLTAN_MEMERR;
  }
  vIdx[0] = 0;
  e = 0;

  for (idx=0; idx < ht_size; idx++){
    hn = hash_table[idx];
    while (hn){
      ZOLTAN_SET_GID(zz, edges + e*numGID, hn->egid);
      hn->firstVert = vIdx[e];
      hn->nextVert  = 0;
      vIdx[e+1] = vIdx[e] + hn->numVerts;
      hn = hn->next;
      e++;
    }
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

  for (v=0; v < numVerts; v++){
    npins = ((v == (numVerts - 1)) ? num_pins : eIdx[v+1]) - eIdx[v];

    for (e=0; e < npins; e++){
      idx = Zoltan_Hash(egid, numGID, (unsigned int)ht_size);
      hn = hash_table[idx];

      while (hn){
        if (ZOLTAN_EQ_GID(zz, hn->egid, egid)){

          ZOLTAN_SET_GID(zz,
             pins + numGID*(hn->firstVert + hn->nextVert),
             vgid);

          hn->nextVert++;
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
  *row_ptr = vIdx;

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}


/*****************************************************************************/
int Zoltan_HG_Graph_Callbacks(
  ZZ *zz,
  ZHG *zhg,                /* Input:   Pointer to Zoltan's structure with
                                       GIDs, LIDs.
                              Output:  Removed edge fields of zhg are changed
                                       if dense edges are removed. */
  int gnVtx,               /* Input:   Global number of vertices in hgraph */
  float esize_threshold,   /* Input:   %age of gnVtx considered a dense edge */
  int return_removed,      /* Input:   flag indicating whether to return
                                       removed edge info */
  int *nedges,             /* Output:  Number of hyperedges on this processor */
  ZOLTAN_ID_PTR *egids,    /* Output:  GIDs of hyperedges on this processor */
  ZOLTAN_ID_PTR *elids,    /* Output:  LIDs of hyperedges on this processor */
  int **esizes,            /* Output:  # of vertices for each hyperedge on
                                       this processor */
  float **ewgts,           /* Output:  edge weights for each hyperedge on
                                       this processor */
  int *npins,              /* Output:  # of pins on this processor = 
                                       sum esizes[i], i = 0..nedges-1. */
  ZOLTAN_ID_PTR *pins,     /* Output:  vertex GIDs of pins */
  int **pin_procs          /* Output:  processors owning pin vertices */
)
{
/* Function to return hypergraph info built from the Zoltan Graph callback 
 * functions. Can be called from serial or parallel hypergraph partitioners.
 * Each vertex has an associated hyperedge containing the vertex and its
 * graph neighbors.  Dense hyperedges are removed.
 */
static char *yo = "Zoltan_HG_Graph_Callbacks";
int ierr = ZOLTAN_OK;
int i, j, k, tmp;
int cnt, ncnt;
int ewgtdim = zz->Edge_Weight_Dim;
int nwgt;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int nvtx = zhg->nObj;                /* Vertex info */
ZOLTAN_ID_PTR vgids = zhg->GIDs;     /* Vertex info */
ZOLTAN_ID_PTR vlids = zhg->LIDs;     /* Vertex info */
ZOLTAN_ID_PTR lid;
float *gewgts = NULL;     /* Graph-edge weights */

  ZOLTAN_TRACE_ENTER(zz, yo);

  *nedges = nvtx;   /* One hyperedge per graph vertex */
  
  if (*nedges > 0) {

    /* Get info about the edges:  GIDs, LIDs, sizes, edge weights */
    *egids = ZOLTAN_MALLOC_GID_ARRAY(zz, *nedges);
    *elids = ZOLTAN_MALLOC_LID_ARRAY(zz, *nedges);
    nwgt = *nedges * ewgtdim;
    if (nwgt) 
      *ewgts = (float *) ZOLTAN_MALLOC(nwgt * sizeof(float));
    if (!*egids || (num_lid_entries && !*elids) ||
        (nwgt && !*ewgts)) MEMORY_ERROR;


    for (i = 0; i < *nedges; i++) {
      /* One hyperedge per graph vertex; use same GIDs/LIDs */
      ZOLTAN_SET_GID(zz, &(*egids)[i*num_gid_entries],
                     &(vgids[i*num_gid_entries]));
      ZOLTAN_SET_LID(zz, &(*elids)[i*num_lid_entries],
                     &(vlids[i*num_lid_entries]));
    }

    ierr = Zoltan_Get_Num_Edges_Per_Obj(zz, nvtx, vgids, vlids, esizes, 
                                        &tmp, npins);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from Zoltan_Get_Num_Edges_Per_Obj");
      goto End;
    }

    /* Remove dense edges from input list */

    ierr = ignore_some_edges(zz, zhg, gnVtx, esize_threshold, return_removed,
                           nedges, *egids, *elids, *esizes, *ewgts, NULL, 1);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from ignore_some_edges.");
      goto End;
    }

    if (*nedges == 0){
      ZOLTAN_FREE(egids);
      ZOLTAN_FREE(elids);
      ZOLTAN_FREE(esizes);
      ZOLTAN_FREE(ewgts);
    }

    /* Now get lists of vertices that are in hyperedges */

    *npins = 0;
    for (i = 0; i < *nedges; i++)
      *npins += (*esizes)[i];

    *pins = ZOLTAN_MALLOC_GID_ARRAY(zz, (*npins + *nedges));
    *pin_procs = (int *) ZOLTAN_MALLOC((*npins + *nedges) * sizeof(int));
    if (ewgtdim)
      gewgts = (float *) ZOLTAN_MALLOC(*npins * ewgtdim * sizeof(float));

    if (!*pins || !*pin_procs || (ewgtdim && !gewgts)) MEMORY_ERROR;

    if (zz->Get_Edge_List_Multi)
      zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data,
                              num_gid_entries, num_lid_entries, *nedges, 
                              *egids, *elids, *esizes, *pins, 
                              *pin_procs, ewgtdim, gewgts, &ierr);
    else {
      cnt = 0;
      for (i = 0; i < *nedges; i++) {
        lid = (num_lid_entries ? &((*elids)[i*num_lid_entries]) : NULL);
        zz->Get_Edge_List(zz->Get_Edge_List_Data,
                          num_gid_entries, num_lid_entries, 
                          &((*egids)[i*num_gid_entries]), lid, 
                          &((*pins)[cnt]), &((*pin_procs)[cnt]), 
                          ewgtdim, &(gewgts[cnt*ewgtdim]), 
                          &ierr);
        cnt += (*esizes)[i];
      }
    }

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from getting Edge_Lists");
      goto End;
    }

    /* Post-process edges to add vgid[i] to hedge i. */
    cnt = *npins;
    *npins += *nedges;  /* Will add vgid[i] to hedge i */
    ncnt = *npins;
    for (i = *nedges-1; i >= 0; i--) {
      /* Copy the existing pins for the edge */
      for (j = 0; j < (*esizes)[i]; j++) {
        cnt--;
        ncnt--;
        ZOLTAN_SET_GID(zz, &((*pins)[ncnt]), &((*pins)[cnt]));
        (*pin_procs)[ncnt] = (*pin_procs)[cnt];
        for (k = 0; k < ewgtdim; k++)   /* sum the graph-edge wgts? */
          (*ewgts)[i*ewgtdim + k] += gewgts[cnt*ewgtdim + k];
      }
      /* Add egid[i] */
      ncnt--;
      ZOLTAN_SET_GID(zz, &((*pins)[ncnt]), &((*egids)[i*num_gid_entries]));
      (*pin_procs)[ncnt] = zz->Proc;
      (*esizes)[i]++;
    }
    ZOLTAN_FREE(&gewgts);
  }
  
End:
  
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
/*****************************************************************************/
static int ignore_some_edges (
  ZZ *zz,
  ZHG *zhg,                /* Input:   Pointer to Zoltan's structure with
                                       GIDs, LIDs.
                              Output:  Removed edge fields of zhg are changed
                                       if dense edges are removed. */
  int gnVtx,               /* Input:   Global number of vertices in hgraph */
  float esize_threshold,   /* Input:   %age of gnVtx considered a dense edge */
  int return_removed,      /* Input:   flag indicating whether to return
                                       removed edge info */
  int *nedges,             /* Output:  Number of hyperedges kept on this proc */
  ZOLTAN_ID_PTR egids,     /* Input:   GIDs of all hyperedges on this proc;
                              Output:  GIDs of kept hyperedges on this proc */
  ZOLTAN_ID_PTR elids,     /* Input:   LIDs of all hyperedges on this proc;
                              Output:  LIDs of kept hyperedges on this proc */
  int *esizes,             /* Input:   # of vtx for each hyperedge on this proc;
                              Output:  # of vtx for each kept hyperedge on
                                       this proc */
  float *ewgts,            /* Input:   For hypergraph input, edge weights for
                                       each hyperedge on this proc;
                              Output:  For hypergraph input, edge weights for 
                                       each kept hyperedge on this proc.
                              Ignored for graph input. */
  ZOLTAN_ID_PTR pins,      /* Input: NULL, or pin GIDs for each edge
                              Output: If not NULL on input, pins remaining */
  int graph_input          /* Input:   Indicates graph input. */
                      
)
{
/* Function to remove dense edges (> esize_threshold vertices)
 * and zero-sized edges (zero vertices) from input data.
 */
char *yo = "ignore_some_edges";
int ierr = ZOLTAN_OK;
int i, j;
int ewgtdim = zz->Edge_Weight_Dim;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int nkeep;                /* Number of edges below gesize_threshold. */
int nremove;              /* Number of edges to be removed; i.e., number of
                             edges above gesize_threshold. */
int nremove_size;         /* Number of pins in removed edges */
ZOLTAN_ID_PTR remove_egids = NULL;  /* Edge GIDs for removed edges */
ZOLTAN_ID_PTR remove_elids = NULL;  /* Edge LIDs for removed edges */
int *remove_esizes = NULL;          /* Edge sizes (# pins) for removed edges */
float *remove_ewgts = NULL;         /* Edge weights for removed edges */
float gesize_threshold;   /* Edges with more vertices than gesize_threshold
                             are considered to be dense. */
ZOLTAN_ID_PTR keep_pins, remove_pins, in_pins;

  /* Remove dense edges and zero-sized edges from input list */
  gesize_threshold = esize_threshold * gnVtx;
  nremove = 0;
  nremove_size = 0;
  for (i = 0; i < *nedges; i++) 
    /* Compute esizes[i]+graph_input; graph_input will add one GID to hedge */
    if (((esizes[i]+graph_input) > gesize_threshold) || 
        ((esizes[i]+graph_input) == 0)){
      nremove++;
      if (pins) nremove_size += esizes[i];
    }

  if (nremove) {
    if (return_removed) {
      /* Keep a record of removed edges so we can get their edge lists
       * later if needed (e.g., to evaluate total partition quality) */
      zhg->nRemove = nremove;
      zhg->Remove_EGIDs = remove_egids = ZOLTAN_MALLOC_GID_ARRAY(zz, nremove);
      zhg->Remove_ELIDs = remove_elids = ZOLTAN_MALLOC_LID_ARRAY(zz, nremove);
      memset(remove_elids, 0, sizeof(int) * nremove * zz->Num_LID);

      zhg->Remove_Esize = remove_esizes 
                        = (int *) ZOLTAN_MALLOC(nremove * sizeof(int));
      if (ewgtdim)
        zhg->Remove_Ewgt = remove_ewgts 
                         = (float *) ZOLTAN_CALLOC(nremove * ewgtdim,
                                                   sizeof(float));

      zhg->Remove_Pin_GIDs = ZOLTAN_MALLOC_LID_ARRAY(zz, nremove_size);

      if (!remove_egids || (num_lid_entries && !remove_elids) 
          || (nremove_size && !zhg->Remove_Pin_GIDs)
          || !remove_esizes || (ewgtdim && !remove_ewgts)) MEMORY_ERROR;
    }
      
    nremove = nkeep = 0;
    keep_pins = pins;
    remove_pins = zhg->Remove_Pin_GIDs;
    in_pins = pins;

    for (i = 0; i < *nedges; i++){
      /* Compare esizes[i]+graph_input; graph_input adds one GID to hedge */
      if (((esizes[i]+graph_input) <= gesize_threshold) &&
          ((esizes[i]+graph_input) > 0)) {
        /* Keep the edge in egids/elids to obtain its pins. */
        if (nkeep != i) {
          ZOLTAN_SET_GID(zz, &(egids[nkeep*num_gid_entries]),
                             &(egids[i*num_gid_entries]));
          if (num_lid_entries)
            ZOLTAN_SET_LID(zz, &(elids[nkeep*num_lid_entries]),
                               &(elids[i*num_lid_entries]));
          esizes[nkeep] = esizes[i];
          if (!graph_input)
            for (j = 0; j < ewgtdim; j++)
              ewgts[nkeep*ewgtdim + j] = ewgts[i*ewgtdim + j];

          if (pins){
            ZOLTAN_COPY_GID_ARRAY(keep_pins, in_pins, zz, esizes[i]);
            keep_pins += (esizes[i] * num_gid_entries);
          }
        }
        nkeep++;
      }
      else if (return_removed) {
        /* Remove the edges from egids/elids; don't want to have to 
           allocate memory for its pins */
        ZOLTAN_SET_GID(zz, &(remove_egids[nremove*num_gid_entries]),
                           &(egids[i*num_gid_entries]));
        if (num_lid_entries)
          ZOLTAN_SET_LID(zz, &(remove_elids[nremove*num_lid_entries]),
                             &(elids[i*num_lid_entries]));
        remove_esizes[nremove] = esizes[i];

        if (pins){
          ZOLTAN_COPY_GID_ARRAY(remove_pins, in_pins, zz, esizes[i]);
          remove_pins += (esizes[i] * num_gid_entries);
        }
        if (!graph_input)
          for (j = 0; j < ewgtdim; j++)
            remove_ewgts[nremove*ewgtdim + j] = ewgts[i*ewgtdim + j];
        nremove++;
      }
      in_pins += (esizes[i] * num_gid_entries);
    }
    *nedges = nkeep;
  }
End:
  return ierr;
}
static int remove_empty_edges(ZZ *zz, int *num_lists, int num_pins, 
   ZOLTAN_ID_PTR edg_GID, int *row_ptr, ZOLTAN_ID_PTR vtx_GID)
{
int nedges = 0;
ZOLTAN_ID_PTR old_e, new_e;
int i, size, npins;

  /* remove edges with no pins from the edge list */

  old_e = new_e = edg_GID;
  npins = 0;

  for (i=0; i<*num_lists; i++){
    size = ((i < *num_lists-1) ? row_ptr[i+1] : num_pins) - row_ptr[i];

    if (size > 0){
      if (new_e < old_e){
        ZOLTAN_SET_GID(zz, new_e, old_e);
        row_ptr[nedges] = npins;
      }
      new_e += zz->Num_GID;
      nedges++;
      npins += size;
    }
    old_e += zz->Num_GID;
  }

  *num_lists = nedges;

  return ZOLTAN_OK;
}

#if VERBOSE_EDGE_INFO
static void show_edges(char *s, ZZ *zz, int num_lists, int num_pins, 
                int *edg_GID, int *row_ptr, int *vtx_GID)
{
int i, j, size, sumsize=0;
int *v = vtx_GID;

  /* helpful in debugging */
  printf("%s> Process %d, %d edges, %d pins\n",s, zz->Proc, num_lists, num_pins);
  for (i=0; i<num_lists; i++){
    size = (i < num_lists-1 ? row_ptr[i+1] : num_pins) - row_ptr[i];
    sumsize += size;
    printf("Edge %d, size %d\n  ", edg_GID[i], size);
    for (j=0; j<size; j++){
      printf("%d ",   *v++);
    }
    printf("\n");
  }
  printf("Sum of edge sizes: %d\n",sumsize);
}
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
