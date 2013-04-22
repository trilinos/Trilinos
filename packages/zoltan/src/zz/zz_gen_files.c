/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "zz_util_const.h"
#include "params_const.h"
#include "phg.h"

#define ZOLTAN_PRINT_VTX_NUM  0  /* print vertex number at beginning of line? */


/* Temporary prototypes. These functions are HG routines
   currently not compiled into Zoltan, but duplicated in this file. */
static ZOLTAN_GNO_TYPE Zoltan_HG_Get_Pins(ZZ *zz, int *nEdges, int **edgeSize,
                   ZOLTAN_ID_PTR *edgeIds, ZOLTAN_ID_PTR *vtxIds, 
                   int *nEwgts, ZOLTAN_ID_PTR *eWgtIds, float **eWgts);
static int turn_off_reduce_dimensions(ZZ *zz);
static ZOLTAN_GNO_TYPE merge_gids(ZZ *zz, ZOLTAN_ID_PTR *merged_egids, ZOLTAN_GNO_TYPE size_merged,
           ZOLTAN_ID_PTR idbuf, ZOLTAN_GNO_TYPE numIds, void *htptr, ZOLTAN_GNO_TYPE htSize);
static int augment_search_structure(ZZ *zz, void *htptr,
     ZOLTAN_GNO_TYPE maxEdges, ZOLTAN_ID_PTR merged_egids, ZOLTAN_GNO_TYPE size_merged, 
     ZOLTAN_GNO_TYPE prev_size_merged);
static ZOLTAN_GNO_TYPE fan_in_edge_global_ids(ZZ *zz, int numEdges, ZOLTAN_ID_PTR egids);

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
 *                 ignored for hypergraphs, Matrix Market standard says "1".
 *   gen_geom,   generate geometry file? 
 *   gen_graph,  generate graph file? 
 *   gen_hg,     generate hypergraph file? 
 *
 *  The output is serialized, such that each processor
 *  will open and close each output file in turn.
 */

  int error=ZOLTAN_OK;

  FILE *fp;
  char full_fname[256];
  int *xadj, *part, *adjproc, *edgeSize;

  int i, j, k, num_obj, num_geom, num_edges, reduce;
  int numPins, edgeOffset=0, vtxOffset=0;
  int print_vtx_num = ZOLTAN_PRINT_VTX_NUM;
  int have_pin_callbacks=0;
  int nEdges=0, nEwgts=0;
  int lenGID = zz->Num_GID;

  ZOLTAN_ID_PTR local_ids=NULL, global_ids=NULL;
  ZOLTAN_ID_PTR edgeIds, vtxIds, eWgtIds, eptr, vptr;
  ZOLTAN_ID_TYPE minid, maxid, minEdgeId, maxEdgeId, minVtxId, maxVtxId;

  ZOLTAN_GNO_TYPE *vtxdist, *adjncy;
  ZOLTAN_GNO_TYPE glob_edges=0, glob_hedges=0;
  ZOLTAN_GNO_TYPE glob_nvtxs, glob_pins, glob_ewgts;
  ZOLTAN_GNO_TYPE gno_val;

  MPI_Datatype zoltan_gno_mpi_type;

  float *float_vwgt, *ewgts, *eWgts, *wptr;
  double *xyz;

  char *yo = "Zoltan_Generate_Files";

  ZOLTAN_TRACE_ENTER(zz, yo);

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  vtxdist = adjncy = NULL;
  xadj = NULL;
  part = edgeSize = NULL;
  adjproc = NULL;
  float_vwgt = ewgts = eWgts = NULL;
  xyz = NULL;
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
    int graph_type = 0;
    /* Build (ParMetis) graph data structures. */
    error = Zoltan_Build_Graph(zz, &graph_type, 1, num_obj,
           global_ids, local_ids, zz->Obj_Weight_Dim, &zz->Edge_Weight_Dim,
           &vtxdist, &xadj, &adjncy, &ewgts, &adjproc);
    if (error != ZOLTAN_OK && error != ZOLTAN_WARN){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_Build_Graph returned error.");
      goto End;
    }
    glob_nvtxs = vtxdist[zz->Num_Proc];
  }
  else{
    /* Compute global number of vertices. */
    gno_val = (ZOLTAN_GNO_TYPE)num_obj;
    MPI_Allreduce(&gno_val, &glob_nvtxs, 1, zoltan_gno_mpi_type, MPI_SUM, zz->Communicator);  
  }

  /* Local number of edges. */
  if (xadj==NULL || xadj[num_obj]==0)
    num_edges = 0;
  else
    num_edges = xadj[num_obj];

  /* Compute global number of edges. */
  gno_val = (ZOLTAN_GNO_TYPE)num_edges;
  MPI_Reduce(&gno_val, &glob_edges, 1, zoltan_gno_mpi_type, MPI_SUM, 0, zz->Communicator);  
  /* Assume no self-edges! */
  glob_edges /= 2;

  /* Build hypergraph, or get hypergraph data. */
  /* We assume edge IDs and vertex (object) IDs are consecutive integers.
   * We'll write out one-based IDs.
   */
  if (gen_hg){
    have_pin_callbacks = zz->Get_HG_Size_CS != NULL && zz->Get_HG_CS != NULL;

    if (!have_pin_callbacks){
      gen_hg = 0;
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

      if (nEdges > 0){
        numPins = edgeSize[nEdges];
      }
      else{
        numPins = 0;
      }

      /* Get the global number of pins for process 0. 
       */
      gno_val = (ZOLTAN_GNO_TYPE)numPins;
      MPI_Reduce(&gno_val, &glob_pins, 1, zoltan_gno_mpi_type, MPI_SUM, 0, zz->Communicator);

      /* Get the global number of edges that weights were
       * provided for.  More than one process may supply
       * weights for a given edge.
       */
      gno_val = (ZOLTAN_GNO_TYPE)nEwgts;
      MPI_Reduce(&nEwgts, &glob_ewgts, 1, zoltan_gno_mpi_type, MPI_SUM, 0, zz->Communicator);

      /* We assume the Edge IDs and Vertex IDs are integers and
       * are contiguous.  Figure out what the lowest ID is.
       * Matrix Market files are one-based so we may adjust IDs.
       */
      minid = glob_hedges;
      maxid = 0;
      for (i=0; i<nEdges; i++){
        if (edgeIds[i*lenGID] < minid) minid = edgeIds[i*lenGID];
        if (edgeIds[i*lenGID] > maxid) maxid = edgeIds[i*lenGID];
      }
      MPI_Allreduce(&minid, &minEdgeId, 1, ZOLTAN_ID_MPI_TYPE, MPI_MIN, zz->Communicator);
      MPI_Allreduce(&maxid, &maxEdgeId, 1, ZOLTAN_ID_MPI_TYPE, MPI_MAX, zz->Communicator);

      minid = glob_nvtxs;
      maxid = 0;
      for (i=0; i<num_obj; i++){
        if (global_ids[i*lenGID] < minid) minid = global_ids[i*lenGID];
        if (global_ids[i*lenGID] > maxid) maxid = global_ids[i*lenGID];
      }
      MPI_Allreduce(&minid, &minVtxId, 1, ZOLTAN_ID_MPI_TYPE, MPI_MIN, zz->Communicator);
      MPI_Allreduce(&maxid, &maxVtxId, 1, ZOLTAN_ID_MPI_TYPE, MPI_MAX, zz->Communicator);

      if ( ((maxVtxId - minVtxId + 1) != glob_nvtxs) ||
           ((maxEdgeId - minEdgeId + 1) != glob_hedges)){
    
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Edge or Vertex IDs not consecutive.");
        error = ZOLTAN_FATAL;
        goto End;
      }

      edgeOffset = 1 - (int)minEdgeId;
      vtxOffset = 1 - (int)minVtxId;
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

    if (num_geom == 0){
      gen_geom = 0;
    }
   }

  Zoltan_Print_Sync_Start(zz->Communicator, 0); 

  /* Do we need this if we have pin callbacks?  Vertex assignment
   * is in the hypergraph file.
   */
  /* Write object assignments to file. */
  /* For now, only write partition number. */
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
  fflush(fp);
  fclose(fp);

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
    fflush(fp);
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
      fprintf(fp, ZOLTAN_GNO_SPEC " " ZOLTAN_GNO_SPEC " %1d%1d%1d", glob_nvtxs, glob_edges, 
        print_vtx_num, (zz->Obj_Weight_Dim>0), (zz->Edge_Weight_Dim>0));
      if (zz->Obj_Weight_Dim>1 || zz->Edge_Weight_Dim>1)
        fprintf(fp, " %d %d", zz->Obj_Weight_Dim, zz->Edge_Weight_Dim);
      fprintf(fp, "\n");
    }


    /* Print edge list for each node (object). */
    for (i=0; i<num_obj; i++){
      /* Print vertex number at beginning of line? */
      if (print_vtx_num){
        fprintf(fp, ZOLTAN_GNO_SPEC " ", vtxdist[zz->Proc]+base_index+i);
      }
      /* First print object (vertex) weight, if any. */
      for (k=0; k<zz->Obj_Weight_Dim; k++)
        fprintf(fp, "%f ", float_vwgt[i*(zz->Obj_Weight_Dim)+k]);
      if (gen_graph){
        /* If graph, then print neighbor list */
        for (j=xadj[i]; j<xadj[i+1]; j++){
          fprintf(fp, ZOLTAN_GNO_SPEC " ", adjncy[j]+base_index);
          /* Also print edge weight, if any. */
          for (k=0; k<zz->Edge_Weight_Dim; k++)
            fprintf(fp, "%f ", ewgts[j*(zz->Edge_Weight_Dim)+k]);
        }
      }
      fprintf(fp, "\n");
    }
    fflush(fp);
    fclose(fp);
  }

  Zoltan_Print_Sync_End(zz->Communicator, 0); 

  /* Separate synchronization for hypergraphs; this could be merged
     into the previous synchronization. */


  /* Write hypergraph to file, if applicable. */

  if (gen_hg){

    /* PLEASE READ: If you change the format of this "matrixmarket plus"
     * file, please change dr_hg_io.c:process_mtxp_file(), which 
     * reads the file for zdrive.
     *
     * This format is our own extension of the NIST Matrix Market file
     * format.  We wished to store vertex and edge weights, and also
     * pin, vertex weight and edge weight ownership data in the file.
     * Here are some rules from their design document:
     *  1. lines are limited to 1024 characters
     *  2. blank lines may appear anywhere after the first line
     *  3. numeric data on a line is separated by one or more blanks
     *  4. real data is in floating-point decimal format, can use "e" notation
     *  5. all indices are 1-based
     *  6. character data may be upper or lower case. 
     */

    sprintf(full_fname, "%s.mtxp", fname);  /* matrixmarket plus */

    Zoltan_Print_Sync_Start(zz->Communicator, 0); 
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
        "%%%%MatrixMarket+ distributed-matrix coordinate real general\n%%\n");
      fprintf(fp, 
        "%%rows = hyperedges, columns = vertices\n%%\n");
      fprintf(fp, 
        "%%#rows #columns #pins #procs dim-vertex-weights "
        "#edge-weight-entries dim-edge-weights\n%%\n");
      fprintf(fp, ZOLTAN_GNO_SPEC " " ZOLTAN_GNO_SPEC " " ZOLTAN_GNO_SPEC 
                  " %d %d " ZOLTAN_GNO_SPEC " %d\n", 
        glob_hedges, glob_nvtxs, glob_pins, 
        zz->Num_Proc,
        zz->Obj_Weight_Dim, glob_ewgts, zz->Edge_Weight_Dim);
      fprintf(fp, 
        "%%\n%%========================================================\n");
      fprintf(fp, 
        "%% Edge ID, Vertex ID, 1.0, Process ID\n%%\n");
      fprintf(fp, 
        "%% Edge and Vertex IDs are 1-based, process IDs are zero-based.\n%%\n");
    }

    eptr = edgeIds;
    vptr = vtxIds;

    for (i=0; i<nEdges; i++){
      for (j=0; j<edgeSize[i]; j++){
        fprintf(fp, ZOLTAN_ID_SPEC " " ZOLTAN_ID_SPEC " 1.0 %d\n",
                eptr[0] + edgeOffset, vptr[0] + vtxOffset, zz->Proc);
        vptr += lenGID;
      }
      eptr += lenGID;
    }
    fflush(fp);
    fclose(fp);
    Zoltan_Print_Sync_End(zz->Communicator, 0); 
    MPI_Barrier(zz->Communicator);

    /* Each proc prints its vertices and vertex weights. */

    Zoltan_Print_Sync_Start(zz->Communicator, 0); 

    vptr = global_ids;
    wptr = float_vwgt;

    fp = fopen(full_fname, "a");

    if (zz->Proc == 0){
      fprintf(fp, 
        "%%\n%%========================================================\n");
      fprintf(fp, 
        "%% Vertex ID, Vertex weights, Process ID\n%%\n");
    }

    for (i=0; i<num_obj; i++){
      fprintf(fp, ZOLTAN_ID_SPEC " ",  vptr[0] + vtxOffset);
      for (j=0; j<zz->Obj_Weight_Dim; j++){
        fprintf(fp, "%f ", *wptr++);
      }
      fprintf(fp, "%d\n",  zz->Proc);
      vptr += lenGID;
    }
    fflush(fp);
    fclose(fp);
    Zoltan_Print_Sync_End(zz->Communicator, 0); 
    MPI_Barrier(zz->Communicator);

    /* Each proc prints its edge weights. */

    if (zz->Edge_Weight_Dim > 0){

      Zoltan_Print_Sync_Start(zz->Communicator, 0); 

      eptr = eWgtIds;
      wptr = eWgts;

      fp = fopen(full_fname, "a");

      if (zz->Proc == 0){
        fprintf(fp, 
          "%%\n%%========================================================\n");
        fprintf(fp, 
          "%% Edge ID, Edge weights, Process ID\n%%\n");
      }

      for (i=0; i<nEwgts; i++){
        fprintf(fp, ZOLTAN_ID_SPEC " ",  eptr[0] + edgeOffset);
        for (j=0; j<zz->Edge_Weight_Dim; j++){
          fprintf(fp, "%f ", *wptr++);
        }
        fprintf(fp, "%d\n",  zz->Proc);
        eptr += lenGID;
      }

      fflush(fp);
      fclose(fp);
      Zoltan_Print_Sync_End(zz->Communicator, 0); 
      MPI_Barrier(zz->Communicator);
    }
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
  
  if (have_pin_callbacks){
    ZOLTAN_FREE(&edgeSize);
    ZOLTAN_FREE(&edgeIds);
    ZOLTAN_FREE(&vtxIds);
    ZOLTAN_FREE(&eWgtIds);
    ZOLTAN_FREE(&eWgts);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return error;
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

static ZOLTAN_GNO_TYPE Zoltan_HG_Get_Pins(ZZ *zz, int *nEdges, int **edgeSize,
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
  ZOLTAN_GNO_TYPE globalNumEdges = -1;
  int numEdges = 0;
  int *esize=NULL;

  /* get edge weights */

  if (dim && zz->Get_HG_Size_Edge_Wts && zz->Get_HG_Edge_Wts){

     zz->Get_HG_Size_Edge_Wts(
                 zz->Get_HG_Size_Edge_Wts_Data, &ew_num_edges, &ierr);

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

      zz->Get_HG_Edge_Wts(zz->Get_HG_Edge_Wts_Data,
                  zz->Num_GID, zz->Num_LID, ew_num_edges, dim,
                  ew_gids, lids, ew_weights, &ierr);

      ZOLTAN_FREE(&lids);
    }

    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "hypergraph edge weight query function");
      goto End;
    }
  }

  /* get pins */

  ierr = Zoltan_Hypergraph_Queries(zz, &numEdges, &num_pins, &egids, &esize, &vgids);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    Zoltan_Multifree(__FILE__,__LINE__,2,&ew_gids,&ew_weights);
    ew_num_edges = 0;
    goto End;
  }

  /* esize array is actually index into vgids, we need size of each edge */

  for (i=0; i<numEdges; i++){
    size = esize[i+1] - esize[i];
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

static ZOLTAN_GNO_TYPE fan_in_edge_global_ids(ZZ *zz, int numEdges, ZOLTAN_ID_PTR egids)
{
ZOLTAN_GNO_TYPE maxEdges, size_merged, prev_size_merged;
int lenGID = zz->Num_GID;
struct _gidht{
  int gidIdx;
  struct _gidht *next;
} *gidNode, *gidNext;
struct _gidht **ht=NULL;
int proc, myrank, nprocs, mask, i; 
unsigned int ui;
unsigned int nbits;
ZOLTAN_ID_PTR merged_egids, idbuf;
ZOLTAN_GNO_TYPE allEdges, numIds, nEdges;
ZOLTAN_GNO_TYPE idbufSize = 0;
int gidTag = 0x1000;  /* any reason tags should not be these values? */
int sizeTag = 0x1001;
MPI_Status stat;
MPI_Datatype zoltan_gno_mpi_type;

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

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

  nEdges = (ZOLTAN_GNO_TYPE)numEdges;

  MPI_Allreduce(&nEdges, &maxEdges, 1, zoltan_gno_mpi_type, MPI_MAX, zz->Communicator);  

  if (maxEdges == 0){
    return 0;
  }

  /* Prepare to create a search structure for my edge global IDs,
   * and any more that are sent to me.  We'll build the search
   * information lazily, just before we need it.
   */

  ht = (struct _gidht **)ZOLTAN_CALLOC(sizeof(struct _gidht *), maxEdges);
  prev_size_merged = 0;

  merged_egids = ZOLTAN_MALLOC_GID_ARRAY(zz, nEdges);
  ZOLTAN_COPY_GID_ARRAY(merged_egids, egids, zz, nEdges);
  size_merged = nEdges;

  /* Do the fan in (bit switching logarithmic fan-in) */

  myrank = zz->Proc;
  nprocs = zz->Num_Proc; 

  for (nbits=0; nbits<(sizeof(int)*8); nbits++){
    if ((nprocs >> nbits) == 0) break;
  }

  mask = 1 << nbits;

  for (ui=0; ui<nbits; ui++){

    mask >>= 1;

    proc = myrank ^ mask;

    if (proc < nprocs){
      if (proc < myrank){
        MPI_Send(&size_merged, 1, zoltan_gno_mpi_type, proc, sizeTag, zz->Communicator);
        if (size_merged > 0){
          MPI_Send(merged_egids, size_merged * lenGID, ZOLTAN_ID_MPI_TYPE, proc, gidTag, zz->Communicator);
        }
      }
 
      else{
        MPI_Recv(&numIds, 1, zoltan_gno_mpi_type, proc, sizeTag, zz->Communicator, &stat);

        if (numIds > 0){

          if (numIds > idbufSize){
            idbuf = ZOLTAN_REALLOC_GID_ARRAY(zz, idbuf, numIds);
            idbufSize = numIds;
          }
          MPI_Recv(idbuf, numIds * lenGID, ZOLTAN_ID_MPI_TYPE, proc, gidTag, zz->Communicator, &stat);

          augment_search_structure(zz, ht, maxEdges, merged_egids, size_merged, prev_size_merged);

          prev_size_merged = size_merged;

          size_merged = merge_gids(zz, &merged_egids, size_merged, idbuf, numIds, ht, maxEdges);
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

  MPI_Bcast(&allEdges, 1, zoltan_gno_mpi_type, 0, zz->Communicator);

  return allEdges;
}
static int augment_search_structure(ZZ *zz, void *htptr,
     ZOLTAN_GNO_TYPE maxEdges, ZOLTAN_ID_PTR merged_egids, ZOLTAN_GNO_TYPE size_merged, 
     ZOLTAN_GNO_TYPE prev_size_merged)
{
struct _gidht{
  int gidIdx;
  struct _gidht *next;
} *gidNode;
struct _gidht **ht=NULL;
int lenGID = zz->Num_GID;
int j;
ZOLTAN_GNO_TYPE i;
ZOLTAN_ID_PTR eptr;

  if (prev_size_merged == size_merged) return ZOLTAN_OK;

  ht = (struct _gidht **)htptr;

  eptr = merged_egids + (prev_size_merged * lenGID);

  for (i=prev_size_merged; i<size_merged; i++){

    j = Zoltan_Hash(eptr, lenGID, (int)maxEdges);

    gidNode = (struct _gidht *)ZOLTAN_MALLOC(sizeof(struct _gidht));
    gidNode->gidIdx = i;

    gidNode->next = ht[j];
    ht[j] = gidNode;

    eptr += lenGID;
  }

  return ZOLTAN_OK;
}
static ZOLTAN_GNO_TYPE merge_gids(ZZ *zz, ZOLTAN_ID_PTR *merged_egids, ZOLTAN_GNO_TYPE size_merged,
           ZOLTAN_ID_PTR idbuf, ZOLTAN_GNO_TYPE numIds, void *htptr, ZOLTAN_GNO_TYPE htSize)
{
struct _gidht{
  int gidIdx;
  struct _gidht *next;
} *gidNode;
struct _gidht **ht=NULL;
ZOLTAN_GNO_TYPE numMerged, i;
int j;
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

    j = Zoltan_Hash(inPtr, lenGID, (int)htSize);

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

