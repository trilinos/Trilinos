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
#include "parmetis_jostle.h"

#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}

#if VERBOSE_EDGE_INFO
static void show_edges(char *s, ZZ *zz, int num_lists, int num_pins, 
                int *edg_GID, int *row_ptr, int *vtx_GID);
#endif

static int convert_to_CRS( ZZ *zz, int num_pins, int *col_ptr,
    int *num_lists, ZOLTAN_ID_PTR *vtx_GID,
    int **row_ptr, ZOLTAN_ID_PTR *edg_GID);

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
/* Store/Search code required to ensure that vertex pairs are
 * only counted once when building hyperedges from neighboring
 * vertices.
 */
typedef struct _gid_node{
  ZOLTAN_ID_PTR gid;
  struct _gid_node *next;
}gid_node;

typedef struct _gid_list{
  gid_node *top;
  gid_node **gn;
  int size;
  int next_slot;
  int lenGID;
}gid_list;

static gid_list seen_gids;

static void initialize_gid_list(ZZ *zz, int nvtx)
{
  seen_gids.top = (gid_node *)ZOLTAN_MALLOC(sizeof(gid_node) * nvtx);
  seen_gids.gn = (gid_node **)ZOLTAN_CALLOC(sizeof(gid_node *) , nvtx);
  seen_gids.size = nvtx;
  seen_gids.next_slot = 0;
  seen_gids.lenGID = zz->Num_GID;
}
static void free_gid_list()
{
  ZOLTAN_FREE(&seen_gids.top);
  ZOLTAN_FREE(&seen_gids.gn);
}
static int find_in_gid_list(ZOLTAN_ID_PTR gid)
{
int i, j, same;
int gidlen = seen_gids.lenGID;
gid_node *gl;

  j = Zoltan_Hash(gid, gidlen, seen_gids.size);

  gl = seen_gids.gn[j];

  while (gl){
    same = 1;
    for (i=0; i<gidlen; i++){
      if (gl->gid[i] != gid[i]){
        same = 0;
        break;
      }
    }
    if (same) return 1;
    gl = gl->next;
  }
  return 0;
}
static void add_to_gid_list(ZOLTAN_ID_PTR gid)
{
int i, j, same;
int gidlen = seen_gids.lenGID;
gid_node *gl;

  j = Zoltan_Hash(gid, gidlen, seen_gids.size);

  gl = seen_gids.gn[j];

  while (gl){
    same = 1;
    for (i=0; i<gidlen; i++){
      if (gl->gid[i] != gid[i]){
        same = 0;
        break;
      }
    }
    if (same) return;
    gl = gl->next;
  }

  gl = seen_gids.top + seen_gids.next_slot;

  gl->gid = gid;
  gl->next = seen_gids.gn[j];

  seen_gids.gn[j] = gl;
  seen_gids.next_slot++;
}

#ifdef DEBUG_GRAPH_TO_HG
static void debug_graph_to_hg(
  int nedges, ZOLTAN_ID_PTR egids, ZOLTAN_ID_PTR elids,
  int *esizes, float *ewgts, int npins,
  ZOLTAN_ID_PTR pins, int *pin_procs, int ewgtdim, int lenGID, int lenLID)
{
  int i,j,k;
  ZOLTAN_ID_PTR nextpin;
  int *nextproc;

  nextpin = pins;
  nextproc = pin_procs;

  printf("%d hyperedges, %d pins\n",nedges,npins);
  for (i=0; i<nedges; i++){
    printf("GID ");
    for (j=0; j<lenGID; j++) printf("%d ", egids[i*lenGID+ j]);
    printf(" LID ");
    for (j=0; j<lenLID; j++) printf("%d ", elids[i*lenLID+ j]);
    printf(" weights ");
    for (j=0; j<ewgtdim; j++) printf("%f ", ewgts[i*ewgtdim+ j]);
    printf(" size %d\n",esizes[i]);

    for (j=0; j < esizes[i]; j++){
      printf("  ");
      for (k=0; k<lenGID; k++) printf("%d ", *nextpin++);
      printf(" (%d), ",*nextproc++);
      if (j && (j%10==0)) printf("\n");
    }
    printf("\n");
  }
}
#endif

int Zoltan_HG_Graph_Callbacks(
  ZZ *zz,
  ZHG *zhg,                /* Input:   Pointer to Zoltan's structure with
                                       GIDs, LIDs.
                              Output:  Removed edge fields of zhg are changed
                                       if dense edges are removed. */
  PHGPartParams *hgp,      /* Input:   Parameters */
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
 *
 * If PHG_FROM_GRAPH_METHOD="neighbors":
 * Each vertex yields a hyperedge containing the vertex and all its
 * graph neighbors.  If a graph has no edges, each resulting hyperedge
 * will contain exactly one vertex.
 *
 * If PHG_FROM_GRAPH_METHOD="pairs":
 * Every pair of vertices that are neighbors in the graph form a hyperedge.
 * If the graph has no edges, the resulting hypergraph will also have
 * no edges.
 *
 * Dense hyperedges are removed.
 */
static char *yo = "Zoltan_HG_Graph_Callbacks";
int ierr = ZOLTAN_OK;
int i, j, k;
int cnt, ncnt;
int ewgtdim = zz->Edge_Weight_Dim;
int nwgt, num_gids;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int nvtx = zhg->nObj;                /* Vertex info */
ZOLTAN_ID_PTR vgids = zhg->GIDs;     /* Vertex info */
ZOLTAN_ID_PTR vlids = zhg->LIDs;     /* Vertex info */
ZOLTAN_ID_PTR vtx_gid, lid_ptr;
int *num_nbors=NULL;
int max_nbors, tot_nbors, use_all_neighbors;
ZOLTAN_ID_PTR nbor_gids, gid_ptr, he_pins;
int *nbor_procs, *proc_ptr, *he_procs;
float *gewgts, *wgt_ptr, *he_wgts;
int num_pins, num_hedges, prefix_sum_hedges;

  ZOLTAN_TRACE_ENTER(zz, yo);

  *nedges = *npins = 0;
  *egids = *elids = *pins = NULL;
  *esizes = *pin_procs = NULL;
  *ewgts = NULL;

  if ((hgp->convert_str[0] == 'n') ||   /* "neighbors" */
      (hgp->convert_str[0] == 'N')){
    use_all_neighbors = 1;
  }
  else{                                 /* "pairs"     */
    use_all_neighbors = 0;
  }

  ierr = Zoltan_Get_Num_Edges_Per_Obj(zz, nvtx, vgids, vlids, 
                 &num_nbors,   /* array of #neighbors for each vertex */
                 &max_nbors,   /* maximum in num_nbors array          */
                 &tot_nbors);  /* total of num_nbors                  */

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Zoltan_Get_Num_Edges_Per_Obj");
    goto End;
  }
  
  if ((tot_nbors == 0) && (!use_all_neighbors)){
    /* 
     * Since there are no connected vertices, there are no hyperedges. 
     */
    ZOLTAN_FREE(&num_nbors); 
    goto End;
  }

  if (use_all_neighbors){
    num_gids = tot_nbors + nvtx;
  }
  else{
    num_gids = tot_nbors;
  }
  gid_ptr = nbor_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, num_gids);
  proc_ptr = nbor_procs = (int *)ZOLTAN_MALLOC(sizeof(int) * num_gids);
  wgt_ptr = gewgts = 
    (float *)ZOLTAN_MALLOC(sizeof(float) * ewgtdim * tot_nbors);

  if ((num_gids && (!nbor_gids || !nbor_procs)) || 
      (tot_nbors && ewgtdim && !gewgts)){
    ierr = ZOLTAN_MEMERR;
    Zoltan_Multifree(__FILE__, __LINE__, 3, &nbor_gids, &nbor_procs, &gewgts);
    ZOLTAN_FREE(&num_nbors);
    goto End;
  }

  if (zz->Get_Edge_List_Multi){

    zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data,
      num_gid_entries, num_lid_entries, nvtx, vgids, vlids, num_nbors,
      nbor_gids, nbor_procs, ewgtdim, gewgts, &ierr);

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in edge list query function");
      Zoltan_Multifree(__FILE__, __LINE__, 3, &nbor_gids, &nbor_procs, &gewgts);
      ZOLTAN_FREE(&num_nbors);
      goto End;
    }
  }
  else{
 
    for (i=0; i<nvtx; i++){

      zz->Get_Edge_List(zz->Get_Edge_List_Data,
        num_gid_entries, num_lid_entries, 
        vgids + (i * num_gid_entries), vlids + (i * num_lid_entries), 
        gid_ptr, proc_ptr, ewgtdim, wgt_ptr,
        &ierr);

      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in edge list query function");
        Zoltan_Multifree(__FILE__,__LINE__,3,&nbor_gids,&nbor_procs,&gewgts);
        ZOLTAN_FREE(&num_nbors);
        goto End;
      }

      gid_ptr += (num_nbors[i] * num_gid_entries);
      proc_ptr += num_nbors[i];
      if (wgt_ptr) wgt_ptr += (ewgtdim * num_nbors[i]);
    }
  }

  if (use_all_neighbors){
    *nedges = nvtx;   /* One hyperedge per graph vertex */
    
    if (*nedges > 0) {
  
      *egids = ZOLTAN_MALLOC_GID_ARRAY(zz, *nedges);
      *elids = ZOLTAN_MALLOC_LID_ARRAY(zz, *nedges);
      nwgt = *nedges * ewgtdim;
      if (nwgt) 
        *ewgts = (float *) ZOLTAN_CALLOC(nwgt, sizeof(float));
      if (!*egids || (num_lid_entries && !*elids) ||
          (nwgt && !*ewgts)) MEMORY_ERROR;

      if (!*egids || !*elids || (nwgt && !*ewgts)){
        ierr = ZOLTAN_MEMERR;
        Zoltan_Multifree(__FILE__,__LINE__,3,&nbor_gids,&nbor_procs,&gewgts);
        Zoltan_Multifree(__FILE__, __LINE__, 3, egids, elids, ewgts);
        goto End;
      }
  
      for (i = 0; i < *nedges; i++) {
        /* One hyperedge per graph vertex; use same GIDs/LIDs */
        ZOLTAN_SET_GID(zz, &(*egids)[i*num_gid_entries],
                       &(vgids[i*num_gid_entries]));
        ZOLTAN_SET_LID(zz, &(*elids)[i*num_lid_entries],
                       &(vlids[i*num_lid_entries]));
      }
  
      /* Post-process edges to add vgid[i] to hedge i. */

      cnt = tot_nbors;
      ncnt = num_gids;

      for (i = *nedges-1; i >= 0; i--) {
        /* Copy the existing pins for the edge */
        for (j = 0; j < num_nbors[i]; j++) {
          cnt--;
          ncnt--;
          ZOLTAN_SET_GID(zz, nbor_gids + (ncnt*num_gid_entries), 
                             nbor_gids + (cnt*num_gid_entries));
          nbor_procs[ncnt] = nbor_procs[cnt];
          for (k = 0; k < ewgtdim; k++)   /* sum the graph-edge wgts? */
            (*ewgts)[i*ewgtdim + k] += gewgts[cnt*ewgtdim + k];
        }
        /* Add egid[i] */
        ncnt--;
        ZOLTAN_SET_GID(zz, nbor_gids + (ncnt*num_gid_entries), 
                           &((*egids)[i*num_gid_entries]));
        nbor_procs[ncnt] = zz->Proc;
        num_nbors[i]++;
      }
      *esizes = num_nbors; 
      *npins = num_gids;
      *pins = nbor_gids;
      *pin_procs = nbor_procs;
      ZOLTAN_FREE(&gewgts);
    }
  }
  else if (!use_all_neighbors){

    /* tot_nbors is an upper bound on number of my hyperedges */

    he_pins = ZOLTAN_MALLOC_GID_ARRAY(zz, tot_nbors*2);
    he_procs = (int *)ZOLTAN_MALLOC(sizeof(int) * tot_nbors*2);
    he_wgts = (float *)ZOLTAN_MALLOC(sizeof(float) * ewgtdim * tot_nbors);

    if (tot_nbors && (!he_pins || !he_procs || (ewgtdim && !he_wgts))){
      ierr = ZOLTAN_MEMERR;
      Zoltan_Multifree(__FILE__,__LINE__,3,&nbor_gids,&nbor_procs,&gewgts);
      Zoltan_Multifree(__FILE__,__LINE__,3,&he_pins,&he_procs,&he_wgts);
      ZOLTAN_FREE(&num_nbors);
      goto End;
    }

    num_hedges = 0;
    num_pins = 0;
    initialize_gid_list(zz, nvtx);
    vtx_gid = vgids;
    gid_ptr = nbor_gids;
    proc_ptr = nbor_procs;
    wgt_ptr = gewgts;

    for (i=0; i<nvtx; i++){
      for (j=0; j < num_nbors[i]; j++){

        /* We'll add this edge if:
         *   The other vertex is on my process and I haven't seen it yet.
         *     OR
         *   The other vertex is on a remote process and my process rank
         *   is lower.  This causes an imbalance, but it will be 
         *   corrected when the 2D decomposition is done.
         */

        if ( ((*proc_ptr == zz->Proc) && !find_in_gid_list(gid_ptr)) || 
              (*proc_ptr > zz->Proc) ){

          ZOLTAN_SET_GID(zz, he_pins + (num_pins * num_gid_entries),
                             vtx_gid);
          ZOLTAN_SET_GID(zz, he_pins + ((num_pins+1) * num_gid_entries),
                             gid_ptr);
 
          he_procs[num_pins]     = zz->Proc;
          he_procs[num_pins + 1] = *proc_ptr;

          if (wgt_ptr){
            for (k=0; k<ewgtdim; k++){
              he_wgts[num_hedges*ewgtdim + k] = wgt_ptr[k];
            }
          }
            
          num_pins += 2;
          num_hedges += 1;
        }

        gid_ptr += num_gid_entries;
        if (wgt_ptr) wgt_ptr += ewgtdim;
        proc_ptr++;
      }

      add_to_gid_list(vtx_gid);
      vtx_gid += num_gid_entries;
    }
    free_gid_list();

    /* We'll assign edge global IDs in order across the processes */
    MPI_Scan(&num_hedges, &prefix_sum_hedges, 1, MPI_INT, MPI_SUM, 
             zz->Communicator);
 
    *esizes = (int *)ZOLTAN_MALLOC(sizeof(int) * num_hedges);
    *egids = ZOLTAN_MALLOC_GID_ARRAY(zz, num_hedges);
    *elids = ZOLTAN_MALLOC_GID_ARRAY(zz, num_hedges);

    if (ewgtdim && *npins && !gewgts) MEMORY_ERROR;
    gid_ptr = lid_ptr = NULL;

    if (num_hedges){
      if (!*esizes|| !*egids || !*elids){
        ierr = ZOLTAN_MEMERR;
        Zoltan_Multifree(__FILE__,__LINE__,3,&nbor_gids,&nbor_procs,&gewgts);
        Zoltan_Multifree(__FILE__,__LINE__,3,&he_pins,&he_procs,&he_wgts);
        Zoltan_Multifree(__FILE__,__LINE__,3,esizes,egids,elids);
        ZOLTAN_FREE(&num_nbors);
        goto End;
      }

      gid_ptr = *egids + (num_gid_entries - 1);
      lid_ptr = *elids + (num_lid_entries - 1);
    }

    k = prefix_sum_hedges - num_hedges; /* sum on processes prior to me */

    for (i=0; i<num_hedges; i++){
      *gid_ptr = k + i;
       gid_ptr += num_gid_entries;

      *lid_ptr = i;
       lid_ptr += num_lid_entries;

      (*esizes)[i] = 2;
    }

    *nedges = num_hedges;
    *npins = num_pins;

    *ewgts = he_wgts;
    *pins = he_pins;
    *pin_procs = he_procs;

    ZOLTAN_FREE(&num_nbors);
    Zoltan_Multifree(__FILE__,__LINE__,3,&nbor_gids,&nbor_procs,&gewgts);
  }

  /* Remove dense edges */

  ierr = Zoltan_HG_ignore_some_edges(zz, zhg, gnVtx,
                esize_threshold, return_removed,
                nedges, *egids, *elids, *esizes, NULL,
                *ewgts, *pins, *pin_procs);

End:

#ifdef DEBUG_GRAPH_TO_HG
  for (i=0; i<zz->Num_Proc; i++){
    if (i == zz->Proc){
      printf("Process %d:\n",i);
      debug_graph_to_hg(*nedges, *egids, *elids,
         *esizes, *ewgts, *npins, *pins, *pin_procs, ewgtdim,
         num_gid_entries, num_lid_entries);
      fflush(stdout);
    }
    MPI_Barrier(zz->Communicator);
  }
#endif

  
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
/*****************************************************************************/
int Zoltan_HG_ignore_some_edges (
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
  int *esizes,             /* Input:   # of vtx on this proc for each hyperedge
                              Output:  # of vtx on this proc for each kept he*/
  int *global_esizes,  /* Input:  sum of vtx over all procs for each hyperedge,
                                  use esizes if NULL
                          Output:  sum of vtx over all procs for each kept he*/
  float *ewgts,            /* Input:   For hypergraph input, edge weights for
                                       each hyperedge on this proc;
                              Output:  For hypergraph input, edge weights for 
                                       each kept hyperedge on this proc.
                              Ignored for graph input. */
  ZOLTAN_ID_PTR pins,      /* Input: NULL, or pin GIDs for each edge
                          Output: If not NULL on input, local pins remaining */
  int *pinProcs       /* Input: NULL, or process owning each pin vertex
                Output: If not NULL on input, processes owning remaining pins*/
)
{
/* Function to remove dense edges (> esize_threshold vertices)
 * and zero-sized edges (zero vertices) from input data.
 */
char *yo = "Zoltan_HG_ignore_some_edges";
int ierr = ZOLTAN_OK;
int i, j;
int ewgtdim = zz->Edge_Weight_Dim;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int nkeep;                /* Number of edges below gesize_threshold. */
int nremove;              /* Number of edges to be removed; i.e., number of
                             edges above gesize_threshold. */
int nremove_size;         /* Number of local pins in removed edges on proc */
ZOLTAN_ID_PTR remove_egids = NULL;  /* Edge GIDs for removed edges */
ZOLTAN_ID_PTR remove_elids = NULL;  /* Edge LIDs for removed edges */
int *use_esizes=NULL; /* Do we have a separate array of global sizes */
int *remove_esizes = NULL;          /* Edge sizes on proc for removed edges */
int *remove_global_esizes = NULL; /* Sum Edge sizes all procs for removed edges */
float *remove_ewgts = NULL;         /* Edge weights for removed edges */
float gesize_threshold;   /* Edges with more vertices than gesize_threshold
                             are considered to be dense. */
ZOLTAN_ID_PTR keep_pins, remove_pins, in_pins;
int *keep_pin_procs, *remove_pin_procs, *in_pin_procs;

  /* Remove dense edges and zero-sized edges from input list */
  gesize_threshold = esize_threshold * gnVtx;
  nremove = 0;
  nremove_size = 0;

  /*
   * If global_esizes is NULL, the size of an edge is just the number
   * of pins on this process.  Otherwise the size across all processes
   * is found in global_esizes.
    */
  if (global_esizes){
    use_esizes=global_esizes;
  }
  else{
    use_esizes=esizes;
  }

  for (i = 0; i < *nedges; i++) {

    if ((use_esizes[i] > gesize_threshold) || (use_esizes[i] == 0)){
      nremove++;
      if (pins || pinProcs){
        nremove_size += esizes[i];
      }
    }
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
      zhg->Remove_GEsize = remove_global_esizes 
                        = (int *) ZOLTAN_MALLOC(nremove * sizeof(int));
      if (ewgtdim)
        zhg->Remove_Ewgt = remove_ewgts 
                         = (float *) ZOLTAN_CALLOC(nremove * ewgtdim,
                                                   sizeof(float));

      if (pins && nremove_size){
        zhg->Remove_Pin_GIDs = ZOLTAN_MALLOC_GID_ARRAY(zz, nremove_size);
        if (!zhg->Remove_Pin_GIDs) MEMORY_ERROR;
      }
      else{
        zhg->Remove_Pin_GIDs = NULL;
      }
      if (pinProcs && nremove_size){
        zhg->Remove_Pin_Procs = (int *)ZOLTAN_MALLOC(sizeof(int) * nremove_size);
        if (!zhg->Remove_Pin_Procs) MEMORY_ERROR;
      }
      else{
        zhg->Remove_Pin_Procs = NULL;
      }

      if (!remove_egids || (num_lid_entries && !remove_elids) 
          || !remove_esizes || !remove_global_esizes
          || (ewgtdim && !remove_ewgts)) MEMORY_ERROR;
    }
      
    nremove = nkeep = 0;
    remove_pins = zhg->Remove_Pin_GIDs;
    keep_pins = pins;
    in_pins = pins;
    remove_pin_procs = zhg->Remove_Pin_Procs;
    keep_pin_procs = pinProcs;
    in_pin_procs = pinProcs;

    for (i = 0; i < *nedges; i++){

      if ((use_esizes[i] <= gesize_threshold) && (use_esizes[i] > 0)) {
        /* Keep the edge in egids/elids to obtain its pins. */
        if (nkeep != i) {
          ZOLTAN_SET_GID(zz, &(egids[nkeep*num_gid_entries]),
                             &(egids[i*num_gid_entries]));
          if (num_lid_entries)
            ZOLTAN_SET_LID(zz, &(elids[nkeep*num_lid_entries]),
                               &(elids[i*num_lid_entries]));
          esizes[nkeep] = esizes[i];
          if (global_esizes){
            global_esizes[nkeep] = global_esizes[i];
          }
          for (j = 0; j < ewgtdim; j++)
            ewgts[nkeep*ewgtdim + j] = ewgts[i*ewgtdim + j];

          if (pins){
            if (keep_pins < in_pins)
              ZOLTAN_COPY_GID_ARRAY(keep_pins, in_pins, zz, esizes[i]);
          }
          if (pinProcs){
            if (keep_pin_procs < in_pin_procs)
              memcpy(keep_pin_procs, in_pin_procs, esizes[i] * sizeof(int));
          }
        }
        nkeep++;
        keep_pins += (esizes[i] * num_gid_entries);
        keep_pin_procs += esizes[i];
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
        remove_global_esizes[nremove] = use_esizes[i];

        if (pins){
          ZOLTAN_COPY_GID_ARRAY(remove_pins, in_pins, zz, esizes[i]);
          remove_pins += (esizes[i] * num_gid_entries);
        }
        if (pinProcs){
          memcpy(remove_pin_procs, in_pin_procs, esizes[i] * sizeof(int));
          remove_pin_procs += esizes[i];
        }
        for (j = 0; j < ewgtdim; j++)
          remove_ewgts[nremove*ewgtdim + j] = ewgts[i*ewgtdim + j];
        nremove++;
      }
      in_pins += (esizes[i] * num_gid_entries);
      in_pin_procs += esizes[i];
    }
    *nedges = nkeep;
  }
End:
  return ierr;
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
