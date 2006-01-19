/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002,2003, Sandia National Laboratories.          *
 * For more info, see the README file in the top-level Zoltan directory.     *   *****************************************************************************/
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
#include "parmetis_jostle.h"
#include "params_const.h"
#include "phg_hypergraph.h"

#define ZOLTAN_PRINT_VTX_NUM  0  /* print vertex number at beginning of line? */


/* Temporary prototypes. These functions are HG routines
   currently not compiled into Zoltan, but duplicated in this file. */
static int Zoltan_HG_Get_Hedges(ZZ *zz, int **p_hindex, 
           ZOLTAN_ID_PTR *p_edge_verts, int **p_edge_procs, 
           float **p_edge_wgts, int *glob_hedges, int *glob_pins);
static int Zoltan_HG_Get_Pins(ZZ *zz, int *nEdges, int **edgeSize,
                   ZOLTAN_ID_PTR *edgeIds, ZOLTAN_ID_PTR *vtxIds, 
                   int *nEwgts, ZOLTAN_ID_PTR *eWgtIds, float **eWgts);
static int Zoltan_HG_Print_Hedges(ZZ *zz, FILE *fp, 
           int *hindex, ZOLTAN_ID_PTR hevtxs, float *hewgts);
static int turn_off_reduce_dimensions(ZZ *zz);
static int merge_gids(ZZ *zz, ZOLTAN_ID_PTR *merged_egids, int size_merged,
           ZOLTAN_ID_PTR idbuf, int numIds, void *htptr, int htSize);
static int augment_search_structure(ZZ *zz, void *htptr,
     int maxEdges, ZOLTAN_ID_PTR merged_egids, int size_merged, 
     int prev_size_merged);
static int fan_in_edge_global_ids(ZZ *zz, int numEdges, ZOLTAN_ID_PTR egids);

/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for generating output files
 *  from Zoltan that describe the data given to Zoltan. These 
 *  functions are all callable by the application, and can be
 *  invoked independent of load-balancing (Zoltan_LB_Balance).
 */
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Generate_Files(ZZ *zz, char *fname, int base_index,
int gen_geom, int gen_graph, int gen_hg)
{
/*
 *  Generate up to four output files:
 *   a) Current assignment of objects to partitions (and procs?)
 *   b) Geometry of the objects, if geometry query functions are registered.
 *   c) Graph if graph query functions are available.
 *   d) Hypergraph if hypergraph query functions are available.
 *
 *  Input parameters:
 *   zz,         pointer to Zoltan struct
 *   fname,      the basename for the output files
 *   base_index, first vertex (object) is labelled 0 or 1?
 *   gen_geom,   generate geometry file? 
 *   gen_graph,  generate graph file? 
 *   gen_hg,     generate hypergraph file? 
 *
 *  The output is serialized, such that each processor
 *  will open and close each output file in turn.
 */

  int error=ZOLTAN_OK;
  ZOLTAN_ID_PTR local_ids = NULL;
  ZOLTAN_ID_PTR global_ids = NULL;
  FILE *fp;
  char full_fname[256];
  int *vtxdist, *xadj, *adjncy, *part, *adjproc, *edgeSize;
  int *heprocs, *hindex;
  ZOLTAN_ID_PTR hevtxs;
  float *float_vwgt, *ewgts, *hewgts, *eWgts, *wptr;
  double *xyz;
  int i, j, k, num_obj, num_geom, num_edges, reduce;
  int glob_nvtxs, glob_edges, glob_hedges, glob_pins, glob_ewgts;
  int print_vtx_num = ZOLTAN_PRINT_VTX_NUM;
  int have_edge_callbacks;
  int have_pin_callbacks;
  int nEdges, nEwgts;
  ZOLTAN_ID_PTR edgeIds, vtxIds, eWgtIds, eptr, vptr;
  float encodeProc;

  char *yo = "Zoltan_Generate_Files";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  vtxdist = xadj = adjncy = part = edgeSize = NULL;
  adjproc = NULL;
  float_vwgt = ewgts = hewgts = eWgts = NULL;
  xyz = NULL;
  heprocs = hindex = NULL;
  hevtxs = NULL;
  edgeIds = vtxIds = eWgtIds = NULL;

  /* Assign default file name if none was given. */
  if (fname==NULL) fname = "noname";

  /* Zoltan_Get_Obj_List allocates memory for all return lists. */
  error = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids,
                              zz->Obj_Weight_Dim, &float_vwgt, &part);
  if (error != ZOLTAN_OK){
    /* Return error */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Get_Obj_List returned error.");
    error = ZOLTAN_FATAL;
    goto End;
  }

  if (gen_graph){
    /* Build (ParMetis) graph data structures. */
    error = Zoltan_Build_Graph(zz, 1, 1, num_obj,
           global_ids, local_ids, zz->Obj_Weight_Dim, zz->Edge_Weight_Dim,
           &vtxdist, &xadj, &adjncy, &ewgts, &adjproc);
    if (error != ZOLTAN_OK && error != ZOLTAN_WARN){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_Build_Graph returned error.");
      goto End;
    }
    glob_nvtxs = vtxdist[zz->Num_Proc];
  }
  else{
    /* Compute global number of vertices. */
    MPI_Reduce(&num_obj, &glob_nvtxs, 1, MPI_INT, MPI_SUM, 0, 
        zz->Communicator);  
  }

  /* Local number of edges. */
  if (xadj==NULL || xadj[num_obj]==0)
    num_edges = 0;
  else
    num_edges = xadj[num_obj];
  /* Compute global number of edges. */
  MPI_Reduce(&num_edges, &glob_edges, 1, MPI_INT, MPI_SUM, 0, 
      zz->Communicator);  
  /* Assume no self-edges! */
  glob_edges /= 2;

  /* Build hypergraph, or get hypergraph data. */
  if (gen_hg){
    have_edge_callbacks =
        zz->Get_Num_HG_Edges != NULL &&
        zz->Get_HG_Edge_Info != NULL &&
        zz->Get_HG_Edge_List != NULL;
    have_pin_callbacks =
        zz->Get_HG_Size_CS != NULL && zz->Get_HG_CS != NULL ;

    if (!have_edge_callbacks && !have_pin_callbacks){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
        "Hypergraph output requested, "
        "but no corresponding query function was found.\n");

      error = ZOLTAN_FATAL;
      goto End;
    }
    if (have_edge_callbacks){
      /* error = Zoltan_HG_Build_Hypergraph(zz, &zhg, NULL); */
      /* Get data in parallel. Zoltan_HG_Build_Hypergraph
         currently only works in serial. */
      error = Zoltan_HG_Get_Hedges(zz, &hindex, &hevtxs, &heprocs, &hewgts,
              &glob_hedges, &glob_pins);
    }
    else{

      /* Use pin callbacks to get pins, edge weight callbacks to
       * get edge weights, calculate global number of edges.
       * Edge weights are not necessarily for edges of my pins.
       */
      glob_hedges = Zoltan_HG_Get_Pins(zz, &nEdges, &edgeSize,
                   &edgeIds, &vtxIds, 
                   &nEwgts, &eWgtIds, &eWgts);

      if (glob_hedges < 0){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
          "Error calling hypergraph pin query functions.\n");
        error = ZOLTAN_FATAL;
        goto End;
      }

      /* get the global number of pins */
       
      for (i=0, j=0; i <nEdges; i++){
        j += edgeSize[i];
      }
    
      MPI_Reduce(&j, &glob_pins, 1, MPI_INT, MPI_SUM, 0, 
          zz->Communicator);  

      /* Get the global number of edges that weights were
       * provided for.  More than one process may supply
       * weights for a given edge.
       */

      MPI_Reduce(&nEwgts, &glob_ewgts, 1, MPI_INT, MPI_SUM, 0, 
          zz->Communicator);  
    }
  }

  /**********************************************************/
  /* Write to files, serialized.                            */
  /* Note: This will be slow (not scalable) for many procs. */
  /**********************************************************/

  if (gen_geom){
    if (zz->Get_Num_Geom == NULL ||
     (zz->Get_Geom == NULL && zz->Get_Geom_Multi == NULL)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
     "Geometry output requested, but no corresponding query function was found.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }
    reduce = turn_off_reduce_dimensions(zz);  /* don't transform coordinates */

    error = Zoltan_Get_Coordinates(zz, num_obj, global_ids, local_ids,
                                   &num_geom, &xyz);

    if (reduce){
      Zoltan_Set_Param(zz, "REDUCE_DIMENSIONS", "1");
    }

    if (error != ZOLTAN_OK && error != ZOLTAN_WARN) {
      goto End;
    }
   }

  Zoltan_Print_Sync_Start(zz->Communicator, 0); 

  /* Write object assignments to file. */
  /* For now, only write partition number. */
  if (num_obj > 0){
    sprintf(full_fname, "%s.assign", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }
    for (i=0; i<num_obj; i++)
      fprintf(fp, "%d\n", part[i]);
      /* fprintf(fp, "%d\n", zz->Proc); */
    fclose(fp);
  }

  /* Write geometry to file, if applicable. */
  if (gen_geom){
    sprintf(full_fname, "%s.coords", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }
    for (i=0; i<num_obj; i++){
      for (j=0; j<num_geom; j++)
        fprintf(fp, "%f ", xyz[i*num_geom + j]);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  /* Write graph to file, if applicable. */
  /* Also create a minimal .graph file with no edges for geometric methods. */
  if (gen_geom || gen_graph){
    sprintf(full_fname, "%s.graph", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }

    /* If proc 0, write first line. */
    if (zz->Proc == 0){
      fprintf(fp, "%% First line: #vertices #edges weight_flag\n");
      fprintf(fp, "%d %d %1d%1d%1d", glob_nvtxs, glob_edges, 
        print_vtx_num, (zz->Obj_Weight_Dim>0), (zz->Edge_Weight_Dim>0));
      if (zz->Obj_Weight_Dim>1 || zz->Edge_Weight_Dim>1)
        fprintf(fp, " %d %d", zz->Obj_Weight_Dim, zz->Edge_Weight_Dim);
      fprintf(fp, "\n");
    }


    /* Print edge list for each node (object). */
    for (i=0; i<num_obj; i++){
      /* Print vertex number at beginning of line? */
      if (print_vtx_num){
        fprintf(fp, "%d ", vtxdist[zz->Proc]+base_index+i);
      }
      /* First print object (vertex) weight, if any. */
      for (k=0; k<zz->Obj_Weight_Dim; k++)
        fprintf(fp, "%f ", float_vwgt[i*(zz->Obj_Weight_Dim)+k]);
      if (gen_graph){
        /* If graph, then print neighbor list */
        for (j=xadj[i]; j<xadj[i+1]; j++){
          fprintf(fp, "%d ", adjncy[j]+base_index);
          /* Also print edge weight, if any. */
          for (k=0; k<zz->Edge_Weight_Dim; k++)
            fprintf(fp, "%f ", ewgts[j*(zz->Edge_Weight_Dim)+k]);
        }
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  Zoltan_Print_Sync_End(zz->Communicator, 0); 

  /* Separate synchronization for hypergraphs; this could be merged
     into the previous synchronization. */


  /* Write hypergraph to file, if applicable. */

  if (gen_hg && have_edge_callbacks){
    Zoltan_Print_Sync_Start(zz->Communicator, 0); 
    if (zz->Get_Num_HG_Edges == NULL || zz->Get_HG_Edge_List == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Hypergraph output requested, but no corresponding query function was found.\n");
      error = ZOLTAN_FATAL;
      Zoltan_Print_Sync_End(zz->Communicator, 0); 
      goto End;
    }
    sprintf(full_fname, "%s.hg", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      Zoltan_Print_Sync_End(zz->Communicator, 0); 
      goto End;
    }

    /* If proc 0, write first line. */
    if (zz->Proc == 0){
      fprintf(fp, "%% First line: #vertices #hyperedges #pins weight_flag\n");
      fprintf(fp, "%d %d %d %1d%1d%1d", glob_nvtxs, glob_hedges, 
        glob_pins, ZOLTAN_PRINT_VTX_NUM, 
        (zz->Obj_Weight_Dim>0), (zz->Edge_Weight_Dim>0));
      if (zz->Obj_Weight_Dim>1 || zz->Edge_Weight_Dim>1)
        fprintf(fp, " %d %d", zz->Obj_Weight_Dim, zz->Edge_Weight_Dim);
      fprintf(fp, "\n");
    }

    /* Each proc prints its part of the hgraph. */
    Zoltan_HG_Print_Hedges(zz, fp, hindex, hevtxs, hewgts);

    fclose(fp);
    Zoltan_Print_Sync_End(zz->Communicator, 0); 
  }

  if (gen_hg && have_pin_callbacks){

    /* Each proc prints its pins. For global IDs we print 1 integer.*/

    Zoltan_Print_Sync_Start(zz->Communicator, 0); 
    sprintf(full_fname, "%s.mtxp", fname);  /* matrixmarket plus */
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");

    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      Zoltan_Print_Sync_End(zz->Communicator, 0); 
      goto End;
    }

    /* If proc 0, write first line. */
    if (zz->Proc == 0){
      fprintf(fp, 
        "%% #hyperedges #vertices #pins dim-vertex-weights "
        "#edge-weight-entries dim-edge-weights\n");
      fprintf(fp, "%d %d %d %d %d %d", glob_hedges, glob_nvtxs,
        glob_pins, zz->Obj_Weight_Dim, glob_ewgts, zz->Edge_Weight_Dim);
      fprintf(fp, "\n");
    }

    eptr = edgeIds;
    vptr = vtxIds;
    encodeProc = (float)zz->Proc + .1;

    fseek(fp, 0, SEEK_END);
    for (i=0; i<nEdges; i++){
      for (j=0; j<edgeSize[i]; j++){
        fprintf(fp, "%d %d %4.1f\n",*eptr, *vptr, encodeProc);
        vptr += zz->Num_GID;
      }
      eptr += zz->Num_GID;
    }
    fflush(fp);
    Zoltan_Print_Sync_End(zz->Communicator, 0); 

    /* Each proc prints its vertex weights. */

    if (zz->Obj_Weight_Dim > 0){
      Zoltan_Print_Sync_Start(zz->Communicator, 0); 

      vptr = global_ids;
      wptr = float_vwgt;

      fseek(fp, 0, SEEK_END);
      for (i=0; i<num_obj; i++){
        fprintf(fp, "%d %d ",  zz->Proc, *vptr);
        for (j=0; j<zz->Obj_Weight_Dim; j++){
          fprintf(fp, "%f ", *wptr++);
        }
        vptr += zz->Num_GID;
        fprintf(fp, "\n");
      }
      fflush(fp);
      Zoltan_Print_Sync_End(zz->Communicator, 0); 
    }

    /* Each proc prints its edge weights. */

    if (zz->Edge_Weight_Dim > 0){

      Zoltan_Print_Sync_Start(zz->Communicator, 0); 

      eptr = eWgtIds;
      wptr = eWgts;

      fseek(fp, 0, SEEK_END);
      for (i=0; i<nEwgts; i++){
        fprintf(fp, "%d %d ",  zz->Proc, *eptr);
        for (j=0; j<zz->Edge_Weight_Dim; j++){
          fprintf(fp, "%f ", *wptr++);
        }
        eptr += zz->Num_GID;
        fprintf(fp, "\n");
      }

      fflush(fp);
      Zoltan_Print_Sync_End(zz->Communicator, 0); 
    }

    fclose(fp);
  }

End:
  ZOLTAN_FREE(&xyz);
  ZOLTAN_FREE(&vtxdist);
  ZOLTAN_FREE(&xadj);
  ZOLTAN_FREE(&adjncy);
  ZOLTAN_FREE(&adjproc);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&float_vwgt);
  ZOLTAN_FREE(&ewgts);
  ZOLTAN_FREE(&part);
  if ( have_edge_callbacks){
    ZOLTAN_FREE(&hindex);
    ZOLTAN_FREE(&hevtxs);
    ZOLTAN_FREE(&heprocs);
    ZOLTAN_FREE(&hewgts);
  }
  if ( have_pin_callbacks){
    ZOLTAN_FREE(&edgeSize);
    ZOLTAN_FREE(&edgeIds);
    ZOLTAN_FREE(&vtxIds);
    ZOLTAN_FREE(&eWgtIds);
    ZOLTAN_FREE(&eWgts);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return error;
}

/************************************************************************ 
 * EBEB: The routines below were copied from the hg directory because   *
 * Zoltan is distributed without hg. Duplicate functions should be      *
 * removed later.                                                       *
 ************************************************************************/

/* Each proc gets hypergraph data from query functions.
 * Compute some global values. 
 * All procs must participate due to collective communication.
 */
static int Zoltan_HG_Get_Hedges(ZZ *zz, int **p_hindex, 
           ZOLTAN_ID_PTR *p_edge_verts, int **p_edge_procs, 
           float **p_edge_wgts, int *glob_hedges, int *glob_pins)
{
  int i, ierr, cnt, j, nEdge, npins, numwgts, minproc;
  int loc_hedges, loc_pins;
  int *hindex = NULL, *edge_procs = NULL;
  ZOLTAN_ID_PTR edge_verts = NULL;
  ZOLTAN_ID_PTR edge_gids = NULL, edge_lids = NULL;
  float *edge_wgts = NULL;
  static char *yo = "Zoltan_HG_Get_Hedges";

  ZOLTAN_TRACE_ENTER(zz, yo);
  ierr = ZOLTAN_OK;

  /* Get hyperedge information from application through query functions. */

  nEdge = zz->Get_Num_HG_Edges(zz->Get_Num_HG_Edges_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_HG_Edges");
    goto End;
  }
  if (nEdge > 0) {
    edge_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, nEdge);
    edge_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, nEdge);
    (*p_hindex) = (int *) ZOLTAN_MALLOC((nEdge+1) * sizeof(int));
    hindex = *p_hindex;
    numwgts = nEdge * zz->Edge_Weight_Dim;
    if (numwgts){
       (*p_edge_wgts) = (float *) ZOLTAN_MALLOC(numwgts * sizeof(float));
       edge_wgts = *p_edge_wgts;
    }
    if (!edge_gids || (zz->Num_LID && !edge_lids) || 
        !hindex || (numwgts && !edge_wgts)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    ierr = zz->Get_HG_Edge_Info(zz->Get_HG_Edge_Info_Data, 
                                zz->Num_GID, zz->Num_LID, nEdge,
                                zz->Edge_Weight_Dim,
                                edge_gids, edge_lids, 
                                hindex, edge_wgts);
    npins = 0;
    for (i = 0; i < nEdge; i++) npins += hindex[i];

    (*p_edge_verts) = ZOLTAN_MALLOC_GID_ARRAY(zz, npins);
    edge_verts = *p_edge_verts;
    (*p_edge_procs) = (int *) ZOLTAN_MALLOC(npins * sizeof(int));
    edge_procs = *p_edge_procs;
    if (npins && (!edge_verts || !edge_procs)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    ierr = zz->Get_HG_Edge_List(zz->Get_HG_Edge_List_Data, 
                                zz->Num_GID, zz->Num_LID, 
                                nEdge, edge_gids, edge_lids,
                                hindex, edge_verts, edge_procs);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,"Error returned from Get_HG_Edge_List");
      goto End;
    }
  }
  /* Build hindex from edge_sizes */
  cnt = 0;
  for (i = 0; i < nEdge; i++) {
    j = hindex[i];
    hindex[i] = cnt;
    cnt += j;
  }
  hindex[nEdge] = cnt;

  /* Compute local hg statistics. */
  /* Make hindex negative if hedge is owned by another proc. */
  loc_hedges = loc_pins = 0;
  for (i = 0; i < nEdge; i++) {
    minproc = zz->Num_Proc;
    for (j=hindex[i]; j<hindex[i+1]; j++)
      if (edge_procs[j]<minproc)
        minproc = edge_procs[j];
    if (minproc == zz->Proc){  /* my hyperedge */
      loc_hedges++;
      loc_pins += (hindex[i+1] - hindex[i]); /* edge_size[i] */
    }
    else  /* lowest proc owns hyperedge, not me */
      hindex[i] = -hindex[i];
  }

  /* Sanity check */
  if (cnt != npins) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Input error:  Number of pins != sum of edge sizes");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /* Compute global #pins and #hyperedges, no duplicates */
  MPI_Allreduce(&loc_hedges, glob_hedges, 1, MPI_INT, MPI_SUM,
      zz->Communicator);
  MPI_Allreduce(&loc_pins, glob_pins, 1, MPI_INT, MPI_SUM,
      zz->Communicator);

End:
  /* Memory will be freed in calling function. */

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}

#define ABS(x) ((x)<0 ? -(x) : (x))

/* Each processor prints its hyperedges to file. */
static int Zoltan_HG_Print_Hedges(ZZ *zz, FILE *fp, 
           int *hindex, ZOLTAN_ID_PTR hevtxs, float *hewgts)
{
  int i, j, ierr, num_edges;
  char *yo = "Zoltan_HG_Print_Hedges";

  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = ZOLTAN_OK;
  num_edges = zz->Get_Num_HG_Edges(zz->Get_Num_HG_Edges_Data, &ierr);

  for (i=0; i<num_edges; i++){
    if (hindex[i]>=0){
      /* Only print hyperedges owned by me (to avoid duplicate hedges) */
      for (j=hindex[i]; j<ABS(hindex[i+1]); j++){ 
        /* EBEB - Print global ids as integers. */
        fprintf(fp, "%d ", (int) hevtxs[j]);
      }
    }
    fprintf(fp, "\n");
  }

  /* Print weights. EBEB - Not yet impl. */

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
static int turn_off_reduce_dimensions(ZZ *zz)
{
  int reduce=0;
  PARAM_VARS param[2] = {
     {"REDUCE_DIMENSIONS", NULL, "INT", 0},
     {NULL, NULL, NULL, 0}};

  Zoltan_Bind_Param(param, "REDUCE_DIMENSIONS", (void *)&reduce);
  Zoltan_Assign_Param_Vals(zz->Params, param, zz->Debug_Level, zz->Proc,
    zz->Debug_Proc);

  if (reduce){
    Zoltan_Set_Param(zz, "REDUCE_DIMENSIONS", "0");
  }

  return reduce;
}

static int Zoltan_HG_Get_Pins(ZZ *zz, int *nEdges, int **edgeSize,
            ZOLTAN_ID_PTR *edgeIds, ZOLTAN_ID_PTR *vtxIds, 
            int *nEwgts, ZOLTAN_ID_PTR *eWgtIds, float **eWgts)
{
  char *yo = "Zoltan_HG_Get_Pins";
  ZOLTAN_ID_PTR ew_gids = NULL, egids=NULL, vgids=NULL;
  ZOLTAN_ID_PTR lids = NULL;
  float *ew_weights = NULL;
  int num_pins, size, i;
  int dim = zz->Edge_Weight_Dim;
  int ew_num_edges = 0;
  int ierr = ZOLTAN_OK;
  int globalNumEdges = -1;
  int numEdges = 0;
  int *esize=NULL;

  /* get edge weights */

  if (dim && zz->Get_HG_Size_Edge_Weights && zz->Get_HG_Edge_Weights){

    ierr = zz->Get_HG_Size_Edge_Weights(
                 zz->Get_HG_Size_Edge_Weights_Data, &ew_num_edges);

    if (((ierr==ZOLTAN_OK)||(ierr==ZOLTAN_WARN)) && (ew_num_edges > 0)){
      ew_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, ew_num_edges);
      lids = ZOLTAN_MALLOC_LID_ARRAY(zz, ew_num_edges);
      ew_weights =
        (float *)ZOLTAN_MALLOC(sizeof(float) * ew_num_edges * dim);

      if (!ew_gids || !lids || !ew_weights){
        Zoltan_Multifree(__FILE__,__LINE__,3,&ew_gids,&lids,&ew_weights);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "memory allocation");
        goto End;
      }

      ierr = zz->Get_HG_Edge_Weights(zz->Get_HG_Edge_Weights_Data,
                  zz->Num_GID, zz->Num_LID, ew_num_edges, dim,
                  ew_gids, lids, ew_weights);

      ZOLTAN_FREE(&lids);
    }

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "hypergraph edge weight query function");
      goto End;
    }
  }

  /* get pins */

  ierr = Zoltan_Call_Hypergraph_Pin_Query(zz, &numEdges, &num_pins,
                &egids, &esize, &vgids);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    Zoltan_Multifree(__FILE__,__LINE__,2,&ew_gids,&ew_weights);
    ew_num_edges = 0;
    goto End;
  }

  /* esize array is actually index into vgids, we need size of each edge */

  for (i=0; i<numEdges; i++){
    size = ((i == numEdges-1) ? num_pins : esize[i+1]) - esize[i];
    esize[i] = size;
  }

  /* figure out how many distinct global edges there are */

  globalNumEdges = fan_in_edge_global_ids(zz, numEdges, egids);

  if (globalNumEdges < 0){
    Zoltan_Multifree(__FILE__,__LINE__,5,
            &ew_gids,&ew_weights,
            &egids, &vgids, &esize);
  }

End:

  *nEdges = numEdges;
  *edgeSize = esize;
  *edgeIds = egids;
  *vtxIds = vgids;
  *nEwgts = ew_num_edges;
  *eWgtIds = ew_gids;
  *eWgts = ew_weights;

  return globalNumEdges;
}

static int fan_in_edge_global_ids(ZZ *zz, int numEdges, ZOLTAN_ID_PTR egids)
{
int maxEdges, size_merged;
int lenGID = zz->Num_GID;
struct _gidht{
  int gidIdx;
  struct _gidht *next;
} *gidNode, *gidNext;
struct _gidht **ht=NULL;
int proc, myrank, nprocs, nbits, mask, i, prev_size_merged, numIds;
ZOLTAN_ID_PTR merged_egids, idbuf;
int allEdges;
int idbufSize = 0;
int gidTag = 0x1000;  /* any reason tags should not be these values? */
int sizeTag = 0x1001;
MPI_Status stat;

  merged_egids = idbuf = NULL;

  if (zz->Num_Proc == 1){
    return numEdges;
  }

  /*
   * The processes have possibly overlapping lists of edge global IDs.
   * We will fan in the lists, with each process merging its
   * IDs with the IDs it received.  Normally we wouldn't want to 
   * allocate storage for all possible edges on one process, because
   * in general this won't fit.  But Zoltan_Generate_Files is used
   * in a debugging mode, so hopefully this won't be a problem.
   */

  MPI_Allreduce(&numEdges, &maxEdges, 1, MPI_INT, MPI_MAX, zz->Communicator);  

  if (maxEdges == 0){
    return 0;
  }

  /* Prepare to create a search structure for my edge global IDs,
   * and any more that are sent to me.  We'll build the search
   * information lazily, just before we need it.
   */

  ht = (struct _gidht **)ZOLTAN_CALLOC(sizeof(struct _gidht *), maxEdges);
  prev_size_merged = 0;

  merged_egids = ZOLTAN_MALLOC_GID_ARRAY(zz, numEdges);
  ZOLTAN_COPY_GID_ARRAY(merged_egids, egids, zz, numEdges);
  size_merged = numEdges;

  /* Do the fan in (bit switching logarithmic fan-in) */

  myrank = zz->Proc;
  nprocs = zz->Num_Proc; 

  for (nbits=0; nbits<(sizeof(int)*8); nbits++){
    if ((nprocs >> nbits) == 0) break;
  }

  mask = 1 << nbits;

  for (i=0; i<nbits; i++){

    mask >>= 1;

    proc = myrank ^ mask;

    if (proc < nprocs){
      if (proc < myrank){
        MPI_Send(&size_merged, 1, MPI_INT, proc, sizeTag, zz->Communicator);
        if (size_merged > 0){
          MPI_Send(merged_egids, size_merged * lenGID, MPI_INT,
                 proc, gidTag, zz->Communicator);
        }
      }
 
      else{
        MPI_Recv(&numIds, 1, MPI_INT, proc, sizeTag, zz->Communicator, &stat);

        if (numIds > 0){

          if (numIds > idbufSize){
            idbuf = ZOLTAN_REALLOC_GID_ARRAY(zz, idbuf, numIds);
            idbufSize = numIds;
          }
          MPI_Recv(idbuf, numIds * lenGID, MPI_INT, 
                   proc, gidTag, zz->Communicator, &stat);

          augment_search_structure(zz, ht, maxEdges,
            merged_egids, size_merged, prev_size_merged);

          prev_size_merged = size_merged;

          size_merged = merge_gids(zz, &merged_egids, size_merged,
                         idbuf, numIds, ht, maxEdges);
        }
      }
    }

    if (myrank >= mask) break;  /* I'm done */
  }

  for (i=0; i<maxEdges; i++){
    gidNode = ht[i];

    while (gidNode){
      gidNext = gidNode->next;
      ZOLTAN_FREE(&gidNode);
      gidNode = gidNext;
    }
  }

  ZOLTAN_FREE(&idbuf);
  ZOLTAN_FREE(&merged_egids);
  ZOLTAN_FREE(&ht);

  /* node zero broadcasts the final number of edge global ids */

  if (myrank == 0){
    allEdges = size_merged;
  }

  MPI_Bcast(&allEdges, 1, MPI_INT, 0, zz->Communicator);

  return allEdges;
}
static int augment_search_structure(ZZ *zz, void *htptr,
     int maxEdges, ZOLTAN_ID_PTR merged_egids, int size_merged, 
     int prev_size_merged)
{
struct _gidht{
  int gidIdx;
  struct _gidht *next;
} *gidNode;
struct _gidht **ht=NULL;
int lenGID = zz->Num_GID;
int i, j;
ZOLTAN_ID_PTR eptr;

  if (prev_size_merged == size_merged) return ZOLTAN_OK;

  ht = (struct _gidht **)htptr;

  eptr = merged_egids + (prev_size_merged * lenGID);

  for (i=prev_size_merged; i<size_merged; i++){

    j = Zoltan_Hash(eptr, lenGID, maxEdges);

    gidNode = (struct _gidht *)ZOLTAN_MALLOC(sizeof(struct _gidht));
    gidNode->gidIdx = i;

    gidNode->next = ht[j];
    ht[j] = gidNode;

    eptr += lenGID;
  }

  return ZOLTAN_OK;
}
static int merge_gids(ZZ *zz, ZOLTAN_ID_PTR *merged_egids, int size_merged,
           ZOLTAN_ID_PTR idbuf, int numIds, void *htptr, int htSize)
{
struct _gidht{
  int gidIdx;
  struct _gidht *next;
} *gidNode;
struct _gidht **ht=NULL;
int numMerged, i, j;
int lenGID = zz->Num_GID;
ZOLTAN_ID_PTR newIds, mergedPtr, inPtr;

  if (numIds < 1){
    return size_merged;
  }

  ht = (struct _gidht **)htptr;

  newIds = ZOLTAN_REALLOC_GID_ARRAY(zz, *merged_egids, size_merged + numIds);

  mergedPtr = newIds + size_merged;
  numMerged = size_merged;
  inPtr = idbuf;

  for (i=0; i<numIds; i++){

    j = Zoltan_Hash(inPtr, lenGID, htSize);

    gidNode = ht[j];
    while (gidNode){
      if (ZOLTAN_EQ_GID(zz, inPtr, newIds + (gidNode->gidIdx * lenGID))){
        break;                           
      }
      gidNode = gidNode->next;
    }

    if (!gidNode){
      ZOLTAN_SET_GID(zz, mergedPtr, inPtr);
      mergedPtr += lenGID;
      numMerged++;
    }
    inPtr += lenGID;
  }

  *merged_egids = newIds;
  return numMerged;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

