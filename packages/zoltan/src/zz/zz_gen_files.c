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
#include "hypergraph.h"
#include "hg.h"

#define ZOLTAN_PRINT_VTX_NUM           0
#define ZOLTAN_PRINT_HGRAPH_FROM_GRAPH 0


/* Temporary prototypes. These functions should perhaps
   move to the hg module?? */
static int Zoltan_HG_Get_Hedges();
static int Zoltan_HG_Print_Hedges();

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

int Zoltan_Generate_Files(ZZ *zz, char *fname, int base_index)
{
/*
 *  Generate up to four output files:
 *   a) Current assignment of objects to partitions (and procs?)
 *   b) Geometry of the objects.
 *   c) Graph if graph query functions are available.
 *   d) Hypergraph if hypergraph query functions are available.
 *
 *  Input parameters:
 *   zz,         pointer to Zoltan struct
 *   fname,      the basename for the output files
 *   base_index, first vertex (object) is labelled 0 or 1?
 *
 *  The output is serialized, such that each processor
 *  will open and close each output file in turn.
 */

  int error=ZOLTAN_OK;
  ZOLTAN_ID_PTR local_ids = NULL;
  ZOLTAN_ID_PTR global_ids = NULL;
  ZHG *zhg;
  FILE *fp;
  char full_fname[256];
  int *vtxdist, *xadj, *adjncy, *part;
  int *vwgt, *adjwgt;
  int *hevtxs, *heprocs, *hindex;
  float *float_vwgt, *ewgts, *hewgts;
  double *xyz;
  int i, j, k, num_obj, num_geom, num_edges, nhedges; 
  int glob_edges, glob_hedges, glob_pins;
  int print_vtx_num = ZOLTAN_PRINT_VTX_NUM;
  char *yo = "Zoltan_Generate_Files";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  vtxdist = xadj = adjncy = vwgt = adjwgt = part = NULL;
  float_vwgt = ewgts = NULL;
  xyz = NULL;

  /* Assign default file name if none was given. */
  if (fname==NULL) fname = "noname";

  /* Zoltan_Get_Obj_List allocates memory for all return lists. */
  error = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids,
                              zz->Obj_Weight_Dim, &float_vwgt, &part);
  if (error){
    /* Return error */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Get_Obj_List returned error.");
    error = ZOLTAN_FATAL;
    goto End;
  }

  /* Build (ParMetis) graph data structures, or just get vtxdist. */
  error = Zoltan_Build_Graph(zz, 1, 1, num_obj,
         global_ids, local_ids, zz->Obj_Weight_Dim, zz->Edge_Weight_Dim,
         &vtxdist, &xadj, &adjncy, &ewgts);
  if (error != ZOLTAN_OK && error != ZOLTAN_WARN){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_Build_Graph returned error.");
    error = ZOLTAN_FATAL;
    goto End;
  }
  /* Local number of edges. */
  if (xadj==NULL || xadj[num_obj]==0)
    num_edges = 0;
  else
    num_edges = xadj[num_obj];
  /* Compute global number of edges. */
  MPI_Reduce(&num_edges, &glob_edges, 1, MPI_INT, MPI_SUM, 0, 
      zz->Communicator);  
  glob_edges /= 2;

  /* Build hypergraph, or get hypergraph data. */
  if (zz->Get_Num_HG_Edges != NULL && zz->Get_HG_Edge_List != NULL) {
    /* error = Zoltan_HG_Build_Hypergraph(zz, &zhg, NULL); */
    /* Get data in parallel. Zoltan_HG_Build_Hypergraph
       currently only works in serial. */
    error = Zoltan_HG_Get_Hedges(zz, &hindex, &hevtxs, &heprocs, &hewgts,
            &glob_hedges, &glob_pins);
  }

  /**********************************************************/
  /* Write to files, serialized.                            */
  /* Note: This will be slow (not scalable) for many procs. */
  /**********************************************************/
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
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_WARN;
      goto End;
    }
    for (i=0; i<num_obj; i++)
      fprintf(fp, "%d\n", part[i]);
      /* fprintf(fp, "%d\n", zz->Proc); */
    fclose(fp);
  }

  /* Write geometry to file, if applicable. */
  if (zz->Get_Num_Geom != NULL && zz->Get_Geom != NULL) {
    num_geom = zz->Get_Num_Geom(zz->Get_Num_Geom_Data, &error);
    if (error != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
             "Error returned from user function Get_Num_Geom.");
      goto End;
    }
    if (num_geom>0)
      xyz = (double *) ZOLTAN_MALLOC (num_geom * sizeof(double));

    sprintf(full_fname, "%s.coords", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_WARN;
      goto End;
    }
    for (i=0; i<num_obj; i++){
      zz->Get_Geom(zz->Get_Geom_Data, zz->Num_GID, zz->Num_LID,
              &(global_ids[i*zz->Num_GID]), &(local_ids[i*zz->Num_LID]),
              xyz, &error);
      if (error == ZOLTAN_FATAL || error == ZOLTAN_MEMERR) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,
             "Error returned from user function Get_Geom.");
        goto End;
      }
      for (j=0; j<num_geom; j++)
        fprintf(fp, "%f ", xyz[j]);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  /* Write graph to file, if applicable. */
  if (zz->Get_Num_Edges != NULL && zz->Get_Edge_List != NULL) {
    sprintf(full_fname, "%s.graph", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_WARN;
      goto End;
    }

    /* If proc 0, write first line. */
    if (zz->Proc == 0){
      fprintf(fp, "% First line: #vertices #edges weight_flag\n");
      fprintf(fp, "%d %d %1d%1d%1d", vtxdist[zz->Num_Proc], glob_edges, 
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
      /* Then print neighbor list */
      for (j=xadj[i]; j<xadj[i+1]; j++){
        fprintf(fp, "%d ", adjncy[j]+base_index);
        /* Also print edge weight, if any. */
        for (k=0; k<zz->Edge_Weight_Dim; k++)
          fprintf(fp, "%f ", ewgts[j*(zz->Edge_Weight_Dim)+k]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  /* Write hypergraph to file, if applicable. */
  if (ZOLTAN_PRINT_HGRAPH_FROM_GRAPH ||
      (zz->Get_Num_HG_Edges != NULL && zz->Get_HG_Edge_List != NULL)) {
    sprintf(full_fname, "%s.hg", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_WARN;
      goto End;
    }

    /* If proc 0, write first line. */
    if (zz->Proc == 0){
      fprintf(fp, "% First line: #vertices #hyperedges #pins weight_flag\n");
      fprintf(fp, "%d %d %d %1d%1d%1d", vtxdist[zz->Num_Proc], glob_hedges, 
        glob_pins, ZOLTAN_PRINT_VTX_NUM, 
        (zz->Obj_Weight_Dim>0), (zz->Edge_Weight_Dim>0));
      if (zz->Obj_Weight_Dim>1 || zz->Edge_Weight_Dim>1)
        fprintf(fp, " %d %d", zz->Obj_Weight_Dim, zz->Edge_Weight_Dim);
      fprintf(fp, "\n");
    }

    /* Each proc prints its part of the hgraph. */
    Zoltan_HG_Print_Hedges(zz, fp, nhedges, hindex, hevtxs, heprocs, hewgts);

    fclose(fp);
  }
  
  Zoltan_Print_Sync_End(zz->Communicator, 0); 

End:
  ZOLTAN_FREE(&xyz);
  ZOLTAN_FREE(&vtxdist);
  ZOLTAN_FREE(&xadj);
  ZOLTAN_FREE(&adjncy);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&float_vwgt);
  ZOLTAN_FREE(&vwgt);
  ZOLTAN_FREE(&ewgts);
  ZOLTAN_FREE(&part);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return error;
}


/* Each proc gets hypergraph data from query functions.
   Compute some global values. */
static int Zoltan_HG_Get_Hedges(ZZ *zz, int **p_hindex, 
           ZOLTAN_ID_PTR *p_edge_verts, int **p_edge_procs, 
           float **p_edge_wgts, int *glob_hedges, int *glob_pins)
{
  int i, ierr, cnt, j, nEdge, nInput, npins, numwgts, minproc;
  int *hindex , *edge_procs;
  ZOLTAN_ID_PTR edge_verts;
  float *edge_wgts;
  static char *yo = "Zoltan_HG_Get_Hedges";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Get hyperedge information from application through query functions. */

  nEdge = zz->Get_Num_HG_Edges(zz->Get_Num_HG_Edges_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_HG_Edges");
    goto End;
  }
  nInput = npins = zz->Get_Num_HG_Pins(zz->Get_Num_HG_Pins_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,"Error returned from Get_Max_HG_Edge_Size");    goto End;
  }
  if (nEdge > 0) {
    numwgts = nEdge * zz->Edge_Weight_Dim;
    edge_verts = ZOLTAN_MALLOC_GID_ARRAY(zz, npins);
    hindex = (int *) ZOLTAN_MALLOC(nEdge * sizeof(int));
    edge_procs = (int *) ZOLTAN_MALLOC(npins * sizeof(int));
    if (numwgts) edge_wgts = (float *) ZOLTAN_MALLOC(numwgts * sizeof(float));
    if (!edge_verts || !hindex || !edge_procs || (numwgts && !edge_wgts)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    ierr = zz->Get_HG_Edge_List(zz->Get_HG_Edge_List_Data, zz->Num_GID,
                                zz->Edge_Weight_Dim, nEdge, npins,
                                hindex, edge_verts, edge_procs, edge_wgts);
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

  /* Make hindex negative if hedge is owned by another proc. */
  for (i = 0; i < nEdge; i++) {
    minproc = zz->Num_Proc;
    for (j=hindex[i]; j<hindex[i+1]; j++)
      if (edge_procs[j]<minproc)
        minproc = edge_procs[j];
    if (minproc != zz->Proc)  /* lowest proc owns hyperedge, not me */
      hindex[i] = -hindex[i];
  }

  /* Sanity check */
  if (cnt != npins) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Input error:  Number of pins != sum of edge sizes");
    goto End;
  }

End:
  /* Return pointers to arrays allocated in this routine. */
  p_hindex = &hindex;
  p_edge_procs = &edge_procs;
  p_edge_verts = &edge_verts;
  p_edge_wgts = &edge_wgts;

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ZOLTAN_OK;
}

#define ABS(x) ((x)<0 ? -(x) : (x))

/* Each processor prints its hyperedges to file. */
static int Zoltan_HG_Print_Hedges(ZZ *zz, FILE *fp, int nhedges,
           int *hindex, ZOLTAN_ID_PTR hevtxs, int heprocs, float *hewgts)
{
  int i,j;
  char *yo = "Zoltan_HG_Print_Hedges";

  ZOLTAN_TRACE_ENTER(zz, yo);

  for (i=0; i<nhedges; i++){
    if (hindex[i]>=0){
      /* EBEB Print all hyperedges, no test for uniqueness. */
      for (j=hindex[i]; j<ABS(hindex[i+1]); j++){ 
        /* Print global ids as integers. */
        fprintf(fp, "%d ", (int) hevtxs[j]);
      }
    }
  }

  /* Print weights. EBEB - Not yet impl. */

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

