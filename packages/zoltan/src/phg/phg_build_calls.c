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

#include "zz_const.h"
#include "zz_util_const.h"
#include "phg_hypergraph.h"
#include "parmetis_jostle.h"
    
#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}

static int ignore_some_edges(ZZ *, ZHG *, int, float, int, int *,
  ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int *, float *, int);
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
static int exchange(ZZ *zz, int proc, void *htptr, 
   ZOLTAN_ID_PTR egid, ZOLTAN_ID_PTR inbuf, int inbufsize,
   int exchange_weights, int ew_num_edges, void *ewhtptr, float *edg_weights,
   char *out_need_list, int out_need_list_size, 
   char *in_need_list, int in_need_list_size,
   float *out_weights, int out_weights_size, 
   float *in_weights, int in_weights_size);
static int do_transfers(ZZ *zz, int proc, void *hn, int numMatches,
  int *myMatch, ZOLTAN_ID_PTR egid, int *yourMatch, ZOLTAN_ID_PTR inbuf);
static float *find_weights(ZZ *zz, 
              void *ewhtptr, int ew_num_edges, ZOLTAN_ID_PTR eid);

/*****************************************************************************/

int Zoltan_HG_Hypergraph_Edge_Callbacks(
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
/* Function to call the Zoltan Hypergraph callback functions.      */
/* Processes return disjoint sets of edges, and know all pins in   */
/* the edge along with process owning the pin vertex.              */
/* Can be called from serial or parallel hypergraph partitioners */

static char *yo = "Zoltan_HG_Hypergraph_Callbacks";
int ierr = ZOLTAN_OK;
int i;
int ewgtdim = zz->Edge_Weight_Dim;
int nwgt;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;

  *nedges = zz->Get_Num_HG_Edges(zz->Get_Num_HG_Edges_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_HG_Edges");
    goto End;
  }


  if (*nedges > 0) {

    /* Get info about the edges:  GIDs, LIDs, sizes, edge weights */
    *egids = ZOLTAN_MALLOC_GID_ARRAY(zz, *nedges);
    *elids = ZOLTAN_MALLOC_LID_ARRAY(zz, *nedges);
    *esizes = (int *) ZOLTAN_MALLOC(*nedges * sizeof(int));
    nwgt = *nedges * ewgtdim;
    if (nwgt) 
      *ewgts = (float *) ZOLTAN_MALLOC(nwgt * sizeof(float));
    if (!*esizes || !*egids || (num_lid_entries && !*elids) ||
        (nwgt && !*ewgts)) MEMORY_ERROR;

    ierr = zz->Get_HG_Edge_Info(zz->Get_HG_Edge_Info_Data,
                                     num_gid_entries, num_lid_entries,
                                     *nedges, ewgtdim, 
                                     *egids, *elids, *esizes, *ewgts); 
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_HG_Edge_Info");
      goto End;
    }
                     
    /* Remove dense edges from input list */
    ierr = ignore_some_edges(zz, zhg, gnVtx, esize_threshold, return_removed,
                             nedges, *egids, *elids, *esizes, *ewgts, 0);

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from ignore_some_edges.");
      goto End;
    }

    /* Now get lists of vertices that are in hyperedges */

    *npins = 0;
    for (i = 0; i < *nedges; i++)
      *npins += (*esizes)[i];

    if (*npins) {
      *pins = ZOLTAN_MALLOC_GID_ARRAY(zz, *npins);
      *pin_procs = (int *) ZOLTAN_MALLOC(*npins * sizeof(int));
      if (!*pins || !*pin_procs) MEMORY_ERROR;

      ierr = zz->Get_HG_Edge_List(zz->Get_HG_Edge_List_Data, 
                                  num_gid_entries, num_lid_entries, *nedges, 
                                  *egids, *elids,
                                  *esizes, *pins, 
                                  *pin_procs);

      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,"Error returned from Get_HG_Edge_List");
        goto End;
      }
    }
  }
  
  /* 
   * KDDKDD -- Assuming hyperedges are given to Zoltan by one processor only.
   * KDDKDD -- Eventually, will loosen that constraint and remove duplicates.
   * KDDKDD -- Or the code might work (although with extra communication)
   * KDDKDD -- with the duplicates.
   * KDDKDD -- Might be easier once we have edge GIDs.
   */

End:
  
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
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
/* for the edges which our pins are part of, so we we need to share    */
/* those while redistributing the pins.                                */

static char *yo = "Zoltan_HG_Hypergraph_Pin_Callbacks";
int ierr = ZOLTAN_OK;
int i, j, have_pins, row_storage, w;
int num_gid_entries = zz->Num_GID;
ZOLTAN_ID_PTR vtx_GID=NULL, edg_GID=NULL, edg_LID, egptr, elptr;
float *edg_weight, *wptr;
int num_lists, num_pins;
int *row_ptr=NULL, *col_ptr, *size_ptr;
struct _ewht{
  ZOLTAN_ID_PTR egid;
  ZOLTAN_ID_PTR elid;
  float *weights;
  struct _ewht *next;
} *ewNodes=NULL, *en;
struct _ewht **ewht = NULL;
int max_need_weights, need_weights, format;
char *need_list=NULL;
int dim = zz->Edge_Weight_Dim;
int lenGID = zz->Num_GID;
int lenLID = zz->Num_LID;

  ewht = (struct _ewht **)htptr;
  if (ewht){
    ewNodes = ewht[ew_num_edges];
  }
  else{
    ewNodes = NULL;
  }

  /* Get size and type of compressed pin storage */

  ierr = zz->Get_HG_Size_CS(zz->Get_HG_Size_CS_Data,
                            &num_lists, &num_pins, &format);

  if ((format != ZOLTAN_COMPRESSED_ROWS)&&(format != ZOLTAN_COMPRESSED_COLS)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Invalid compression format returned in Get_HG_Size_CS");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  have_pins = ((num_lists > 0) && (num_pins > 0));
  row_storage = (format == ZOLTAN_COMPRESSED_ROWS);

  /* Get the hypergraph pins in compressed storage format */

  if (have_pins){
    if (!row_storage){    /* compressed column storage */

      vtx_GID = ZOLTAN_MALLOC_GID_ARRAY(zz, num_lists); 
      col_ptr = (int *)ZOLTAN_MALLOC(num_lists * sizeof(int));
      edg_GID = ZOLTAN_MALLOC_GID_ARRAY(zz, num_pins);

      if (!vtx_GID || !col_ptr || !edg_GID){
        Zoltan_Multifree(__FILE__, __LINE__, 3, &vtx_GID, &col_ptr, &edg_GID);
        MEMORY_ERROR;
      }
  
      ierr = zz->Get_HG_CS(zz->Get_HG_CS_Data, num_gid_entries,
               num_lists, num_pins, format,
               vtx_GID, col_ptr, edg_GID);

      if ((ierr == ZOLTAN_OK) || (ierr == ZOLTAN_WARN)){
        ierr = convert_to_CRS(zz, 
                     num_pins,    /* number of pins doesn't change */
                     col_ptr,   
                     &num_lists,  /* replace with number of rows */
                     &vtx_GID,    /* replace with pins           */
                     &row_ptr,    /* index into start of each row in vtx_GID */
                     &edg_GID);   /* replace with row (edge) GIDs */
      }

      ZOLTAN_FREE(&col_ptr);
    }
    else{               /* compressed row storage */

      edg_GID = ZOLTAN_MALLOC_GID_ARRAY(zz, num_lists); 
      row_ptr = (int *)ZOLTAN_MALLOC(num_lists * sizeof(int));
      vtx_GID = ZOLTAN_MALLOC_GID_ARRAY(zz, num_pins);

      if (!vtx_GID || !row_ptr || !edg_GID){
        Zoltan_Multifree(__FILE__, __LINE__, 3, &vtx_GID, &row_ptr, &edg_GID);
        MEMORY_ERROR;
      }

      ierr = zz->Get_HG_CS(zz->Get_HG_CS_Data, num_gid_entries,
                 num_lists, num_pins, format,
                 edg_GID, row_ptr, vtx_GID);

    }
  }

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    Zoltan_Multifree(__FILE__, __LINE__, 3, &vtx_GID, &row_ptr, &edg_GID);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
       "Error calling hypergraph compressed pin storage query");
    goto End;
  }

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
       * we need to process them.
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
          j = Zoltan_Hash(egptr, lenGID, (unsigned int)ew_num_edges);
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

  ierr = 
    distribute_edges(zz, &num_lists, &num_pins, 
                     &edg_GID, &edg_LID, &row_ptr, &vtx_GID,
                     need_weights, max_need_weights, need_list,
                     ew_num_edges, htptr, &edg_weight);

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
  size_ptr = *esizes = (int *) ZOLTAN_MALLOC(num_lists * sizeof(int));

  for (i=0; i<num_lists-1; i++){
    size_ptr[i] = row_ptr[i+1] - row_ptr[i];
  }
  size_ptr[num_lists-1] = num_pins - row_ptr[num_lists - 1];

  ZOLTAN_FREE(&row_ptr);

  *npins = num_pins;
  *pins = vtx_GID;
                     
  /* Remove dense edges from input list */
  ierr = ignore_some_edges(zz, zhg, gnVtx, esize_threshold, return_removed,
                           nedges, edg_GID, edg_LID, *esizes, edg_weight, 0);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from ignore_some_edges.");
    goto End;
  }

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
int i, j, maxEdges, nedges, npins, change, rc;
int *rptr;
ZOLTAN_ID_PTR egid, eptr, elidptr, vptr, inbuf;
struct _egidNode {
  ZOLTAN_ID_PTR egidBuf;
  int idx;
  ZOLTAN_ID_PTR pins;
  char pinsAreCopies;
  struct _egidNode *next;
} *egidNode, *en;
struct _egidNode **ht;  /* hash table for my pin edges */
int l0, l1, lLen, myL, left;
int r0, r1, rLen, myR;
int relativeRank, inbufsize, have_pins, esize, v, e; 
int in_need_list_size=0, in_weights_size=0, out_weights_size=0;
char *in_need_list=NULL;
float *in_weights=NULL, *out_weights=NULL;

/* If there is an edge such that more than one process has pins for */
/* that edge, assign the edge to only one process.  Transfer all    */
/* pins to that process.                                            */

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
   * than one process owns pins for the same edge.  Pins may
   * or may not overlap.)
   **Goal: Assign each edge to only one process.  Each process
   * owns all pins for each of it's edges.
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
   * all pins for the edge it owns.  This exchange and decision
   * must be completed before beginning the exchange with the next
   * process.  And a process is never assigned an edge that it
   * did not have at the beginning of the exchange.  
   *
   * At the end of nprocs-1 exchanges, each edge is owned by
   * only one process.
   */

  /*
   * Divide processes roughly in half.  All processes in the
   * left half will perform an exchange with all processes in
   * the right half.  Recurse, dividing left and right roughly
   * in half.
   */

  myL = 0;
  myR = nprocs-1;
  change = 0;

  while ((myR - myL > 0) && (ierr != ZOLTAN_FATAL)){

    left = divideInterval(rank, myL, myR, &l0, &l1, &r0, &r1);

    if (left){                /* I'm in the left half */
      rLen = r1 - r0 + 1;
      relativeRank = rank - l0;

      for (i=0; i<rLen; i++){
        j = r0 + ((relativeRank + i) % rLen);
        rc = exchange(zz, j, (void *)ht, egid, inbuf, inbufsize,
               max_need_weights, ew_table_size, (void *)htptr, *edg_weights,
               out_need_list, nedges,
               in_need_list, in_need_list_size,
               out_weights, out_weights_size,
               in_weights, in_weights_size);

        if (rc > 0){
          change = 1;
        }
        else if (rc < 0){
          ierr = ZOLTAN_FATAL;
          break;
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
        rc = exchange(zz, j, (void *)ht, egid, inbuf, inbufsize,
               max_need_weights, ew_table_size, (void *)htptr, *edg_weights,
               out_need_list, nedges,
               in_need_list, in_need_list_size,
               out_weights, out_weights_size,
               in_weights, in_weights_size);

        if (rc > 0){
          change = 1;
        }
        else if (rc < 0){
          ierr = ZOLTAN_FATAL;
          break;
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

  if (change && (ierr != ZOLTAN_FATAL)){ /* Must rewrite my edge lists */

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

      eptr = *edg_GID = ZOLTAN_REALLOC_GID_ARRAY(zz, *edg_GID, nedges);
      rptr = *row_ptr = (int *)ZOLTAN_REALLOC(*row_ptr, nedges);

      elidptr = ZOLTAN_MALLOC_LID_ARRAY(zz, nedges);
      vptr = ZOLTAN_MALLOC_GID_ARRAY(zz, npins);

      if (!eptr || !rptr || !elidptr || !vptr){
        Zoltan_Multifree(__FILE__, __LINE__, 4, 
          &eptr, &rptr, &elidptr, &vptr);
        MEMORY_ERROR;
      }

      for (i=0, v=0, e=0; i< egid[0]; i++) {
        en = egidNode + i;
        npins = en->egidBuf[lenGID];

        if (npins > 0){
          ZOLTAN_SET_GID(zz, eptr + (e*lenGID), en->egidBuf);
          ZOLTAN_SET_LID(zz, elidptr + (e*lenLID), *edg_LID + (i*lenLID));
          ZOLTAN_COPY_GID_ARRAY(vptr + (v*lenGID), en->pins, zz, npins);
          rptr[e] = v;
          e++;
          v += npins;
        }
      }
      ZOLTAN_FREE(edg_LID);
      *edg_LID = elidptr;
      ZOLTAN_FREE(vtx_GID);
      *vtx_GID = vptr;
    }
    else{
      *num_lists = 0;
      *num_pins = 0;
      ZOLTAN_FREE(edg_GID);
      ZOLTAN_FREE(vtx_GID);
      ZOLTAN_FREE(row_ptr);
    }
  }

End:

  for (i=0; i< egid[0]; i++) {
    if (egidNode[i].pinsAreCopies){
      ZOLTAN_FREE(&egidNode[i].pins);
    }
  }
  ZOLTAN_FREE(&egidNode);
  ZOLTAN_FREE(&egid);

  return ierr;
}
static int divideInterval(int rank, int myL, int myR, 
                   int *l0, int *l1, int *r0, int *r1)
{
int inLeftHalf= 0;
int center;

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
static int exchange(ZZ *zz, int proc, void *htptr, 
     ZOLTAN_ID_PTR egid, ZOLTAN_ID_PTR inbuf, int inbufsize,
     int exchange_weights, int ew_num_edges, void *ewhtptr, float *edg_weights,
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
int *matchA, *matchB, *Bidx, nMatches;
int i, j, k, w, nweights_found=0, lengthA;
MPI_Request req, need_list_req, weight_req;
MPI_Status status;
ZOLTAN_ID_PTR eptr, eid, gidsA, gidsB;
float *wptr=NULL, *wgts;
int dim = zz->Edge_Weight_Dim;
int ierr=ZOLTAN_OK, ret_val = 0;

  /* TODO put comments at top of every function                */
  /* Returns 1 if my edges lists changed during this exchange, */
  /* 0 if my edges did not change, -1 on error.                */
  /* If only edge weights changed, still returns 0 because no  */
  /* arrays need to be re-written.                             */

  ht = (struct _egidNode **)htptr;

  /* post receive for other process' edge GIDs, and possibly edge weights
      needed array and edge weights response */

  MPI_Irecv(inbuf, inbufsize * lenGID, MPI_INT, proc, 
            tag, zz->Communicator, &req);

  if (exchange_weights){
    MPI_Irecv(in_need_list, in_need_list_size, MPI_CHAR, proc, need_tag,
                 zz->Communicator, &need_list_req);
    MPI_Irecv(in_weights, in_weights_size, MPI_FLOAT, proc, weight_tag,
                 zz->Communicator, &weight_req);
  }

  /* send my edge GIDs to that process, along with possibly my edge
       weights needed flags */

  numMyEdges = egid[0];

  MPI_Send(egid, (2*numMyEdges + 1) * lenGID, MPI_INT, proc, 
            tag, zz->Communicator);

  if (exchange_weights){
    MPI_Send(out_need_list, out_need_list_size, MPI_CHAR, proc,
               need_tag, zz->Communicator);
  }


  /* await GIDs, then process to see if we have any in common,
       and possibly if we have any edge weights they need */

  MPI_Wait(&req, &status);

  if (exchange_weights){
    MPI_Wait(&need_list_req, &status);
    wptr = out_weights + 1;
    nweights_found = 0;
  }

  eptr = inbuf;
  numYourEdges = eptr[0];  /* TODO sanity check here */

  j = ((numYourEdges < numMyEdges) ? numYourEdges : numMyEdges);

  if ((j <= 0) && !exchange_weights){
    return 0;   /* no change to my edge lists */
  }

  myMatch   = (int *)ZOLTAN_MALLOC(j * sizeof(int));
  yourMatch = (int *)ZOLTAN_MALLOC(j * sizeof(int));

  if (!myMatch || !yourMatch){
    Zoltan_Multifree(__FILE__, __LINE__, 2, &myMatch, &yourMatch);
    MEMORY_ERROR
  }

  eptr += lenGID;
  nMatches = 0;

  for (i=0; i<numYourEdges; i++){

    eid         = eptr;

    if (exchange_weights && in_need_list[i]){
      wgts = find_weights(zz, ewhtptr, ew_num_edges, eid);
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

    if (numYourPins > 0){
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

    MPI_Send(out_weights, (nweights_found * (dim+1)) + 1, MPI_FLOAT, proc,
             weight_tag, zz->Communicator);

    MPI_Wait(&weight_req, &status);

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
    ZOLTAN_FREE(&myMatch);
    ZOLTAN_FREE(&yourMatch);
    return 0;  /* no change to my edge lists */
  }

  ret_val = 1;
  
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

  Bidx = (int *)ZOLTAN_MALLOC(numAEdges * sizeof(int)); 

  if (!Bidx) MEMORY_ERROR

  for (i=0; i<nMatches; i++){
    j = (2*matchA[i] + 2) * lenGID;
    gidsA[j] = -gidsA[j];          /* make pins negative to flag it */
    Bidx[matchA[i]] = matchB[i];   /* location of B's matching edge */
  }

  lengthA = numAEdges;
  k = 2 * lenGID;

  for (i=0; i<lengthA; i++){

    if (gidsA[k] < 0){

      if (numAEdges < numBEdges){
        j = (2 * Bidx[i] + 2) * lenGID;  /* matching edge in other array */
        gidsB[j] = -gidsB[j];            /* flag that B sends to A */
        gidsA[k] = -gidsA[k];            /* rather than A sends to B */
        numBEdges--;
      }
      else{
        numAEdges--;
      }
    }
    k += (2 * lenGID);
  }

  ZOLTAN_FREE(&Bidx);

  /* exchange pins, update my edge and pin information */

  ierr = do_transfers(zz, proc, (void *)ht, 
            nMatches, myMatch, egid, yourMatch, inbuf);

End:

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    ret_val = -1;
  }

  ZOLTAN_FREE(&myMatch);
  ZOLTAN_FREE(&yourMatch);

  return ret_val; 
}
static float *find_weights(ZZ *zz, 
              void *ewhtptr, int ew_num_edges, ZOLTAN_ID_PTR eid)
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

  ewht = (struct _ewht **)ewhtptr;

  idx = Zoltan_Hash(eid, zz->Num_GID, (unsigned int)ew_num_edges);

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
static int do_transfers(ZZ *zz, int proc, void *hn, int numMatches,
  int *myMatch, ZOLTAN_ID_PTR egid, int *yourMatch, ZOLTAN_ID_PTR inbuf)
{
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
int recvSize, sendSize, i, j, nrecvs, idx, npins, nedges, nMyPins;
int tagPins = 0xf0f0;
ZOLTAN_ID_PTR recvBuf=NULL, sendBuf=NULL, pinBuf=NULL;
MPI_Request reqPins;
MPI_Status status;

  /* TODO check for duplicate pins and filter them out */
  /* This is expensive, so maybe we should only do it if they */
  /* request this service with a parameter.                   */

  egidTable = (struct _egidNode **)hn;
  nedges = egid[0];
  egidNodes = egidTable[nedges];

  /* Calculate message sizes */

  recvSize = 0;
  sendSize = 0;
  nrecvs = 0;

  for (i=0; i<numMatches; i++){
    j = (2 * yourMatch[i] + 2) * lenGID;

    if (inbuf[j] < 0){
      recvSize += (2 - inbuf[j]);
      nrecvs++;
    }
    else{
      j = (2 * myMatch[i] + 2) * lenGID;
      sendSize += (2 - egid[j]);
    }
  }

  /* Post receive for incoming edge pins */

  if (recvSize > 0){
    recvBuf = (int *)ZOLTAN_MALLOC_GID_ARRAY(zz, recvSize);
    MPI_Irecv(recvBuf, recvSize * lenGID, MPI_INT, proc, 
              tagPins, zz->Communicator, &reqPins);
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
    sendBuf = (int *)ZOLTAN_MALLOC_GID_ARRAY(zz, sendSize);
    eptr = sendBuf;

    for (i=0; i<numMatches; i++){
      j = (2 * myMatch[i] + 2) * lenGID;
      npins = egid[j];

      if (npins < 0){       /* I give this edge to other proc */

        egid[j] = 0;             /* now I have no pins in edge */
        npins = -npins;

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

    MPI_Send(sendBuf, sendSize*lenGID, MPI_INT, proc, tagPins,
             zz->Communicator);

    ZOLTAN_FREE(&sendBuf);
  }

  /* Add pins received from other process to my list.  */

  if (recvSize > 0){
    MPI_Wait(&reqPins, &status);

    eptr = recvBuf;

    for (i=0; i<nrecvs; i++){
      idx = *eptr;
      eptr += lenGID;
      npins = *eptr;
      eptr += lenGID;

      if (npins > 0){
        en = egidNodes + idx;  
        nMyPins = en->egidBuf[lenGID];
        pinBuf = (int *)ZOLTAN_MALLOC_GID_ARRAY(zz, nMyPins + npins);
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

  ZOLTAN_FREE(&recvBuf);

  return ZOLTAN_OK;
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
int maxEdges = num_pins;
int numVerts = *num_lists;
int numEdges, ierr;
ZOLTAN_ID_PTR egid, vgid;
int v, e, idx, found, pin;
struct _hash_node {
  ZOLTAN_ID_PTR egid;
  char *verts;
  struct _hash_node *next;
} *hn, *tmp;
struct _hash_node **hash_table; 
ZOLTAN_ID_PTR ccs_egids=NULL, ccs_vgids=NULL;
int *p;

  ierr = ZOLTAN_OK;

  /*
   * Convert from CCS to CRS (compressed columns to compressed rows)
   * We have the lists of edges for each vertex.  Create lists of
   * vertices for each edge.
   */

  hash_table = 
    (struct _hash_node **)ZOLTAN_CALLOC(maxEdges, sizeof(struct _hash_node *));

  if (!hash_table){
    return ZOLTAN_MEMERR;
  }

  egid = *edg_GID;

  for (v=0; v<numVerts; v++){
    numEdges = ((v == numVerts-1) ? 
              (num_pins - col_ptr[v]) : (col_ptr[v+1] - col_ptr[v])); 

    for (e=0; e<numEdges; e++){

       idx = Zoltan_Hash(egid, zz->Num_GID, (unsigned int)maxEdges);
       found = 0;
       hn = hash_table[idx];

       while (hn){
         if (ZOLTAN_EQ_GID(zz, hn->egid, egid)){
           hn->verts[v] = 1;
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
         hn->verts = (char *)ZOLTAN_CALLOC(numVerts, sizeof(char));
         if (!hn->verts){
           ierr = ZOLTAN_MEMERR;
           goto End;
         }
         hn->verts[v] = 1;
         hn->next = hash_table[idx];
         hash_table[idx] = hn;
       }

       egid += zz->Num_GID;
    }
  }

  /* Build compressed row arrays from tables just calculated */

  numEdges = 0;
  for (idx=0; idx<maxEdges; idx++){
    hn = hash_table[idx];
    while (hn){
      numEdges++;
      hn = hn->next;
      if (numEdges > maxEdges){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid row/column information");
        ierr = ZOLTAN_FATAL;
        goto End;
      }
    }
  }

  p = (int *)ZOLTAN_MALLOC(numEdges * sizeof(int));
  egid = ccs_egids = ZOLTAN_MALLOC_GID_ARRAY(zz, numEdges);
  vgid = ccs_vgids = ZOLTAN_MALLOC_GID_ARRAY(zz, num_pins);

  if (!p || !egid || !vgid){
    Zoltan_Multifree(__FILE__, __LINE__, 3, &p, &egid, &vgid);
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  for (idx=0, e=0, v=0; idx<maxEdges; idx++){
    hn = hash_table[idx];
    while (hn){

      p[e] = v;
      ZOLTAN_SET_GID(zz, egid, hn->egid);

      for (pin=0; pin<numVerts; pin++){
        if (hn->verts[pin]){
          ZOLTAN_SET_GID(zz, vgid, *vtx_GID + pin * zz->Num_GID);
          vgid += zz->Num_LID;
          v++;
        }
      }
      egid += zz->Num_GID;
      e++;
      hn = hn->next;
    }
  }

End:
  for (idx=0; idx<maxEdges; idx++){
    hn = hash_table[idx];
    while (hn){
      ZOLTAN_FREE(&hn->verts);
      tmp = hn;
      hn = hn->next;
      ZOLTAN_FREE(&tmp);
    }
  }
  ZOLTAN_FREE(&hash_table);

  *num_lists = numEdges;
  ZOLTAN_FREE(vtx_GID);
  *vtx_GID = ccs_vgids;
  ZOLTAN_FREE(edg_GID);
  *edg_GID = ccs_egids;
  *row_ptr = p;

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
                              nedges, *egids, *elids, *esizes, *ewgts, 1);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from ignore_some_edges.");
      goto End;
    }

    /* Now get lists of vertices that are in hyperedges */

    *npins = 0;
    for (i = 0; i < *nedges; i++)
      *npins += (*esizes)[i];

    if (*npins) {
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
ZOLTAN_ID_PTR remove_egids = NULL;  /* Edge GIDs for removed edges */
ZOLTAN_ID_PTR remove_elids = NULL;  /* Edge LIDs for removed edges */
int *remove_esizes = NULL;          /* Edge sizes (# pins) for removed edges */
float *remove_ewgts = NULL;         /* Edge weights for removed edges */
float gesize_threshold;   /* Edges with more vertices than gesize_threshold
                             are considered to be dense. */

  /* Remove dense edges and zero-sized edges from input list */
  gesize_threshold = esize_threshold * gnVtx;
  nremove = 0;
  for (i = 0; i < *nedges; i++) 
    /* Compute esizes[i]+graph_input; graph_input will add one GID to hedge */
    if (((esizes[i]+graph_input) > gesize_threshold) || 
        ((esizes[i]+graph_input) == 0))
      nremove++;

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

      if (!remove_egids || (num_lid_entries && !remove_elids) 
          || !remove_esizes || (ewgtdim && !remove_ewgts)) MEMORY_ERROR;
    }
      
    nremove = nkeep = 0;
    for (i = 0; i < *nedges; i++)
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
        if (!graph_input)
          for (j = 0; j < ewgtdim; j++)
            remove_ewgts[nremove*ewgtdim + j] = ewgts[i*ewgtdim + j];
        nremove++;
      }
    *nedges = nkeep;
  }
End:
  return ierr;
}
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
