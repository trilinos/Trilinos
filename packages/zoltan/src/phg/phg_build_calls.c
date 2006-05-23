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
      *ewgts = (float *) ZOLTAN_CALLOC(nwgt, sizeof(float));
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

    ierr = Zoltan_HG_ignore_some_edges(zz, zhg, gnVtx, esize_threshold, 
      return_removed, nedges, *egids, *elids, *esizes, NULL, *ewgts, 
      NULL, NULL, 1);

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from Zoltan_HG_ignore_some_edges.");
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
      
    if (*npins + *nedges) {
      *pins = ZOLTAN_MALLOC_GID_ARRAY(zz, (*npins + *nedges));
      *pin_procs = (int *) ZOLTAN_MALLOC((*npins + *nedges) * sizeof(int));
      if (!*pins || !*pin_procs)
        MEMORY_ERROR;
    }
    if (ewgtdim)
      gewgts = (float *) ZOLTAN_MALLOC(*npins * ewgtdim * sizeof(float));

    if (ewgtdim && !gewgts) MEMORY_ERROR;

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
  int *pinProcs,      /* Input: NULL, or process owning each pin vertex
                Output: If not NULL on input, processes owning remaining pins*/
  int graph_input          /* Input:   Indicates graph input. */
                      
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
    /* Compute esizes[i]+graph_input; graph_input will add one GID to hedge */
    if (((use_esizes[i]+graph_input) > gesize_threshold) || 
        ((use_esizes[i]+graph_input) == 0)){
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
      /* Compare esizes[i]+graph_input; graph_input adds one GID to hedge */
      if (((use_esizes[i]+graph_input) <= gesize_threshold) &&
          ((use_esizes[i]+graph_input) > 0)) {
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
          if (!graph_input)
            for (j = 0; j < ewgtdim; j++)
              ewgts[nkeep*ewgtdim + j] = ewgts[i*ewgtdim + j];

          if (pins){
            if (keep_pins < in_pins)
              ZOLTAN_COPY_GID_ARRAY(keep_pins, in_pins, zz, esizes[i]);
            keep_pins += (esizes[i] * num_gid_entries);
          }
          if (pinProcs){
            if (keep_pin_procs < in_pin_procs)
              memcpy(keep_pin_procs, in_pin_procs, esizes[i] * sizeof(int));
            keep_pin_procs += esizes[i];
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
        remove_global_esizes[nremove] = use_esizes[i];

        if (pins){
          ZOLTAN_COPY_GID_ARRAY(remove_pins, in_pins, zz, esizes[i]);
          remove_pins += (esizes[i] * num_gid_entries);
        }
        if (pinProcs){
          memcpy(remove_pin_procs, in_pin_procs, esizes[i] * sizeof(int));
          remove_pin_procs += esizes[i];
        }
        if (!graph_input)
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
