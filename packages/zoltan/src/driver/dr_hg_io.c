// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>

#include "dr_const.h"
#include "dr_externs.h"
#include "dr_input_const.h"
#include "dr_util_const.h"
#include "dr_par_util_const.h"
#include "dr_err_const.h"
#include "dr_output_const.h"
#include "dr_elem_util_const.h"
#include "dr_maps_const.h"
#include "ch_input_const.h"
#include "ch_init_dist_const.h"
#include "dr_hg_readfile.h"
#include "dr_compress_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#ifndef MAX_STR_LENGTH
#define MAX_STR_LENGTH 80
#endif

#define IS_BLANK(c) ((c == '\t') || (c == ' '))

#define COMMENT_CHAR '%'

#define PRINT_DEBUG_INFO 0

#if PRINT_DEBUG_INFO
static void debug_elements(int Proc, int Num_Proc, int num, ELEM_INFO_PTR el);
static void debug_lists(int Proc, int Num_Proc, int nedge, int *index, ZOLTAN_ID_TYPE *vtx, 
  int *vtx_proc, ZOLTAN_ID_TYPE *egid);
static void debug_pins(int Proc, int Num_Proc,  
          ZOLTAN_ID_TYPE nGlobalEdges, ZOLTAN_ID_TYPE nGlobalVtxs, 
          int vtxWDim, int edgeWDim,
          int nMyPins, ZOLTAN_ID_TYPE *myPinI, ZOLTAN_ID_TYPE *myPinJ,
          int nMyVtx, ZOLTAN_ID_TYPE *myVtxNum, float *myVtxWgts,
          int nMyEdgeWgts, ZOLTAN_ID_TYPE *myEWGno, float *myEdgeWgts);
#endif

static int dist_hyperedges( MPI_Comm, PARIO_INFO_PTR, int, int, int, 
    int *, int *, int **, int **, int **, int **, int *, float  **, short  *);

static int create_edge_lists(int, ZOLTAN_ID_TYPE *, ZOLTAN_ID_TYPE *,
      int *, ZOLTAN_ID_TYPE **, int **, ZOLTAN_ID_TYPE **);

static int process_mtxp_file(PARIO_INFO_PTR, char *, size_t , int , int , int, 
    ZOLTAN_ID_TYPE *, ZOLTAN_ID_TYPE *, int *, int *, int *, 
    ZOLTAN_ID_TYPE **, ZOLTAN_ID_TYPE **, int *, ZOLTAN_ID_TYPE **, 
    float **, int *, ZOLTAN_ID_TYPE **, float **, int);

static int read_mtxp_lines(char *, int, int, char **, int *);
static char *get_nth_token(char *, int , size_t , int , char );
static char *first_char(char *, size_t);
static void make_string(char *, char *);
static void section_starts(char *, char **, char **, char **);
static char *next_line(char *, size_t);
static char *next_line_of_data(char *, size_t, char *);
static char *next_line_of_comment(char *, size_t, char *);

static int my_vtx(int , ZOLTAN_ID_TYPE , ZOLTAN_ID_TYPE , ZOLTAN_ID_TYPE ,
       int , int , PARIO_INFO_PTR );
static int my_pin(ZOLTAN_ID_TYPE , ZOLTAN_ID_TYPE , int ,
       ZOLTAN_ID_TYPE , ZOLTAN_ID_TYPE , ZOLTAN_ID_TYPE , ZOLTAN_ID_TYPE ,
       int , int , PARIO_INFO_PTR );

static int _zoltan_sortFunc(const void *, const void *);

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* Read in MatrixMarket file.
 * For now, just call the hypergraph routine.
 * In the future we may want to allow both graphs and hypergraphs
 * to be read in MatrixMarket format.
 *
 */
int read_mm_file(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh
)
{
  return read_hypergraph_file(Proc, Num_Proc, prob, pio_info, mesh);
}

/* Read from file and set up hypergraph. */
int read_hypergraph_file(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh
)
{
  /* Local declarations. */
  const char  *yo = "read_hypergraph_file";
  char   cmesg[FILENAME_MAX+256];

  int    i, distributed_pins = 0, vertex, nextEdge;
  int    nvtxs = 0, nhedges = 0, npins = 0;
  int    vwgt_dim=0, hewgt_dim=0, vtx, edgeSize;
  int   *hindex = NULL, *hvertex_proc = NULL;
  int   *hgid = NULL, *hvertex = NULL;
  int edge, global_npins, gnhedges;
  int   gnvtxs;
  float *hewgts = NULL, *vwgts = NULL;
  ZOLTAN_FILE* fp = NULL;
  int base = 0;   /* Smallest vertex number; usually zero or one. */
  char filename[FILENAME_MAX+9];

  /* Variables that allow graph-based functions to be reused. */
  /* If no chaco.graph or chaco.coords files exist, values are NULL or 0,
   * since graph is not being built. If chaco.graph and/or chaco.coords
   * exist, these arrays are filled and values stored in mesh.
   * Including these files allows for comparison of HG methods with other
   * methods, along with visualization of results and comparison of
   * LB_Eval results.
   */
  int    ch_nvtxs = 0;        /* Temporary values for chaco_read_graph.   */
#ifdef KDDKDD
  int    ch_vwgt_dim = 0;     /* Their values are ignored, as vertex      */
#endif
  float *ch_vwgts = NULL;     /* info is provided by hypergraph file.     */
  int   *ch_start = NULL, *ch_adj = NULL, ch_ewgt_dim = 0;
  short *ch_assignments = NULL;
  float *ch_ewgts = NULL;
  int    ch_ndim = 0;
  float *ch_x = NULL, *ch_y = NULL, *ch_z = NULL;
  int    ch_no_geom = TRUE;   /* Assume no geometry info is given; reset if
				 it is provided. */
  int    file_error = 0;

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (Proc == 0) {

    /* Open and read the hypergraph file. */
    if (pio_info->file_type == HYPERGRAPH_FILE)
      sprintf(filename, "%s.hg", pio_info->pexo_fname);
    else if (pio_info->file_type == MATRIXMARKET_FILE)
      sprintf(filename, "%s.mtx", pio_info->pexo_fname);
    else {
	sprintf(cmesg, "fatal:  invalid file type %d", pio_info->file_type);
	Gen_Error(0, cmesg);
	return 0;
    }

    fp = ZOLTAN_FILE_open(filename, "r", pio_info->file_comp);
    file_error = (fp == NULL);
  }



  MPI_Bcast(&file_error, 1, MPI_INT, 0, zoltan_get_global_comm());

  if (file_error) {
    sprintf(cmesg,
      "fatal:  Could not open hypergraph file %s",pio_info->pexo_fname);
    Gen_Error(0, cmesg);
    return 0;
  }

  if (pio_info->file_type == HYPERGRAPH_FILE) {
    /* read the array in on processor 0 */
    if (Proc == 0) {
      if (HG_readfile(Proc, fp, &nvtxs, &nhedges, &npins,
		    &hindex, &hvertex, &vwgt_dim, &vwgts,
		    &hewgt_dim, &hewgts, &base) != 0){
	Gen_Error(0, "fatal: Error returned from HG_readfile");
	return 0;
      }
    }
  }
  else if (pio_info->file_type == MATRIXMARKET_FILE) {
    /*
     * pio_info->chunk_reader == 0  (the usual case)
     *   process 0 will read entire file in MM_readfile,
     *   and will distribute vertices in chaco_dist_graph and pins in
     *   dist_hyperedges later.   (distributed_pins==0)
     *
     * pio_info->chunk_reader == 1  ("initial read = chunks" in zdrive.inp)
     *   process 0 will read the file in chunks, and will send vertices
     *   and pins to other processes before reading the next chunk, all
     *   in MM_readfile.  (distributed_pins==1)
     */

    if (MM_readfile(Proc, Num_Proc, fp, pio_info,
		    &nvtxs,     /* global number of vertices */
		    &nhedges,   /* global number of hyperedges */
		    &npins,     /* local number of pins */
		    &hindex, &hvertex, &vwgt_dim, &vwgts,
		    &hewgt_dim, &hewgts, &ch_start, &ch_adj,
		    &ch_ewgt_dim, &ch_ewgts, &base, &global_npins)) {
      Gen_Error(0, "fatal: Error returned from MM_readfile");
      return 0;
    }

    if (Proc == 0) ZOLTAN_FILE_close(fp);

    if ((Num_Proc > 1) && pio_info->chunk_reader && (global_npins > Num_Proc)){
      distributed_pins = 1;
    }
    else{
      distributed_pins = 0;
    }
  }


#ifdef KDDKDD
 {
   /* If CHACO graph file is available, read it. */

   sprintf(filename, "%s.graph", pio_info->pexo_fname);

   fp = ZOLTAN_FILE_open(filename, "r", pio_info->file_comp);
   file_error =
#ifndef ZOLTAN_COMPRESS
     (fp == NULL);
#else
   fp.error;
#endif


   if (!file_error) {
      /* CHACO graph file is available. */
      /* Assuming hypergraph vertices are same as chaco vertices. */
      /* Chaco vertices and their weights are ignored in rest of function. */
      if (chaco_input_graph(fp, filename, &ch_start, &ch_adj, &ch_nvtxs,
		      &ch_vwgt_dim, &ch_vwgts, &ch_ewgt_dim, &ch_ewgts) != 0) {
	Gen_Error(0, "fatal: Error returned from chaco_input_graph");
	return 0;
      }
    }
   else
     ch_nvtxs = nvtxs;


    /* If coordinate file is available, read it. */
   sprintf(filename, "%s.coords", pio_info->pexo_fname);

   fp = ZOLTAN_FILE_open(filename, "r", pio_info->file_comp);
   file_error =
#ifndef ZOLTAN_COMPRESS
     (fp == NULL);
#else
   fp.error;
#endif

    if (!file_error) {
      /* CHACO coordinates file is available. */
      ch_no_geom = FALSE;
      if (chaco_input_geom(fpkdd, filename, ch_nvtxs, &ch_ndim,
			   &ch_x, &ch_y, &ch_z) != 0) {
	Gen_Error(0, "fatal: Error returned from chaco_input_geom");
	return 0;
      }
    }
 }
#else /* KDDKDD */
  ch_nvtxs = nvtxs;
#endif /* KDDKDD */


  {
     /* Read Chaco assignment file, if requested */
   if (pio_info->init_dist_type == INITIAL_FILE) {
     sprintf(filename, "%s.assign", pio_info->pexo_fname);

   fp = ZOLTAN_FILE_open(filename, "r", pio_info->file_comp);

   if (fp == NULL) {
     sprintf(cmesg, "Error:  Could not open Chaco assignment file %s; "
	     "initial distribution cannot be read",
	     filename);
     Gen_Error(0, cmesg);
     return 0;
   }
   else {
     /* read the coordinates in on processor 0 */
     ch_assignments = (short *) malloc(nvtxs * sizeof(short));
     if (nvtxs && !ch_assignments) {
       Gen_Error(0, "fatal: memory error in read_hypergraph_file");
       return 0;
     }
     /* closes fpassign when done */
     if (chaco_input_assign(fp, filename, ch_nvtxs, ch_assignments) != 0){
       Gen_Error(0, "fatal: Error returned from chaco_input_assign");
       return 0;
     }
   }
   }
 }

  MPI_Bcast(&base, 1, MPI_INT, 0, zoltan_get_global_comm());

  if (distributed_pins){
    gnhedges = nhedges;
    nhedges = 0;
    hewgt_dim = 0;
    hewgts = NULL;
    for (edge=0; edge<gnhedges; edge++){
      edgeSize = hindex[edge+1] - hindex[edge];
      if (edgeSize > 0) nhedges++;
    }
    hgid = (int *)malloc(nhedges * sizeof(int));
    hvertex_proc = (int *)malloc(npins * sizeof(int));
    nextEdge=0;
    vtx=0;
    for (edge=0; edge<gnhedges; edge++){
      edgeSize = hindex[edge+1] - hindex[edge];
      if (edgeSize > 0){
	hgid[nextEdge] = edge+1;
	if (nextEdge < edge){
	  hindex[nextEdge+1] = hindex[nextEdge] + edgeSize;
	}
	for (vertex=0; vertex<edgeSize; vertex++,vtx++){
	  hvertex_proc[vtx] = ch_dist_proc(hvertex[vtx], NULL, 1);
	}
	nextEdge++;
      }
    }
    gnvtxs = nvtxs;
    nvtxs = ch_dist_num_vtx(Proc, NULL);
    if (ch_start){    /* need to include only vertices this process owns */
      for (i=0,vertex=0; i<gnvtxs; i++){
	if ((ch_start[i+1] > ch_start[vertex]) || /* vtx has adjacencies so it's mine */
	    (ch_dist_proc(i, NULL, 0) == Proc))   /* my vtx with no adjacencies */
	  {
	  if (i > vertex){
	    ch_start[vertex+1] = ch_start[i+1];
	  }
	  vertex++;
	}
      }
    }
#if PRINT_DEBUG_INFO
    debug_lists(Proc, Num_Proc, nhedges, hindex, hvertex, hvertex_proc, hgid);
#endif
  } else{

    /* Distribute hypergraph graph */
    /* Use hypergraph vertex information and chaco edge information. */

    if (!chaco_dist_graph(zoltan_get_global_comm(), pio_info, 0, &gnvtxs, &nvtxs,
	     &ch_start, &ch_adj, &vwgt_dim, &vwgts, &ch_ewgt_dim, &ch_ewgts,
	     &ch_ndim, &ch_x, &ch_y, &ch_z, &ch_assignments)) {
      Gen_Error(0, "fatal: Error returned from chaco_dist_graph");
      return 0;
    }

    if (!dist_hyperedges(zoltan_get_global_comm(), pio_info, 0, base, gnvtxs, &gnhedges,
		       &nhedges, &hgid, &hindex, &hvertex, &hvertex_proc,
		       &hewgt_dim, &hewgts, ch_assignments)) {
      Gen_Error(0, "fatal: Error returned from dist_hyperedges");
      return 0;
    }
  }

  /* Initialize mesh structure for Hypergraph. */
  mesh->data_type = ZOLTAN_HYPERGRAPH;
  mesh->num_elems = nvtxs;
  mesh->vwgt_dim = vwgt_dim;
  mesh->ewgt_dim = ch_ewgt_dim;
  mesh->elem_array_len = mesh->num_elems + 5;
  mesh->num_dims = ch_ndim;
  mesh->num_el_blks = 1;

  mesh->gnhedges = (ZOLTAN_ID_TYPE)gnhedges;
  mesh->nhedges = nhedges;
  mesh->hewgt_dim = hewgt_dim;

  mesh->hindex = hindex;
  mesh->hvertex_proc = hvertex_proc;
  mesh->heNumWgts = nhedges;
  mesh->hewgts = hewgts;

  mesh->hgid = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nhedges);
  mesh->hvertex = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * hindex[nhedges]);

  for (i=0; i < nhedges; i++){
    mesh->hgid[i] = (ZOLTAN_ID_TYPE)hgid[i];
  }
  free(hgid);

  for (i=0; i < hindex[nhedges]; i++){
    mesh->hvertex[i] = (ZOLTAN_ID_TYPE)hvertex[i];
  }
  free(hvertex);

  mesh->heWgtId = NULL;

  mesh->eb_etypes = (int *) malloc (4 * mesh->num_el_blks * sizeof(int));
  if (!mesh->eb_etypes) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  mesh->eb_ids = mesh->eb_etypes + mesh->num_el_blks;
  mesh->eb_nnodes = mesh->eb_ids + mesh->num_el_blks;
  mesh->eb_nattrs = mesh->eb_nnodes + mesh->num_el_blks;

  mesh->eb_cnts = (ZOLTAN_ID_TYPE *) malloc ( mesh->num_el_blks * sizeof(ZOLTAN_ID_TYPE));

  mesh->eb_names = (char **) malloc (mesh->num_el_blks * sizeof(char *));
  if (!mesh->eb_names) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  mesh->eb_etypes[0] = -1;
  mesh->eb_ids[0] = 1;
  mesh->eb_cnts[0] = (ZOLTAN_ID_TYPE)nvtxs;
  mesh->eb_nattrs[0] = 0;
  /*
   * Each element has one set of coordinates (i.e., node) if a coords file
   * was provided; zero otherwise.
   */
  MPI_Bcast( &ch_no_geom, 1, MPI_INT, 0, zoltan_get_global_comm());
  if (ch_no_geom)
    mesh->eb_nnodes[0] = 0;
  else
    mesh->eb_nnodes[0] = 1;

  /* allocate space for name */
  mesh->eb_names[0] = (char *) malloc((MAX_STR_LENGTH+1) * sizeof(char));
  if (!mesh->eb_names[0]) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  strcpy(mesh->eb_names[0], "hypergraph");

  /* allocate the element structure array */
  mesh->elements = (ELEM_INFO_PTR) malloc (mesh->elem_array_len * sizeof(ELEM_INFO));
  if (!(mesh->elements)) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /*
   * initialize all of the element structs as unused by
   * setting the globalID to ZOLTAN_ID_INVALID 
   */
  for (i = 0; i < mesh->elem_array_len; i++)
    initialize_element(&(mesh->elements[i]));

  /*
   * now fill the element structure array with the
   * information from the Chaco file
   * Use hypergraph vertex information and chaco edge information.
   */

  if (!chaco_fill_elements(Proc, Num_Proc, prob, mesh, pio_info, (ZOLTAN_ID_TYPE)gnvtxs, nvtxs,
		     ch_start, ch_adj, vwgt_dim, vwgts, ch_ewgt_dim, ch_ewgts,
		     ch_ndim, ch_x, ch_y, ch_z, ch_assignments, base)) {
    Gen_Error(0, "fatal: Error returned from chaco_fill_elements");
    return 0;
  }
#if PRINT_DEBUG_INFO
  debug_elements(Proc, Num_Proc, mesh->num_elems,mesh->elements);
#endif

  safe_free((void **)(void *) &vwgts);
  safe_free((void **)(void *) &ch_ewgts);
  safe_free((void **)(void *) &ch_vwgts);
  safe_free((void **)(void *) &ch_x);
  safe_free((void **)(void *) &ch_y);
  safe_free((void **)(void *) &ch_z);
  safe_free((void **)(void *) &ch_start);
  safe_free((void **)(void *) &ch_adj);
  safe_free((void **)(void *) &ch_assignments);

 if (Debug_Driver > 3)
   print_distributed_mesh(Proc, Num_Proc, mesh);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/* Read from mtxp file (such as that written out by Zoltan_Generate_Files)
 * and set up the hypergraph. */

int read_mtxplus_file(
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh
)
{
  /* Local declarations. */
  const char  *yo = "read_mtxplus_file";
  char filename[FILENAME_MAX+9], cmesg[FILENAME_MAX+256];
  struct stat statbuf;
  int rc, fsize, i, j;
  char *filebuf=NULL;
  FILE* fp;
  ZOLTAN_ID_TYPE nGlobalEdges, nGlobalVtxs; 
  ZOLTAN_ID_TYPE *myPinI, *myPinJ, *myVtxNum, *myEWGno;
  ZOLTAN_ID_TYPE *edgeGno, *pinGno;
  int vtxWDim, edgeWDim;
  int nMyPins, nMyVtx, nMyEdgeWgts;
  float *myVtxWgts, *myEdgeWgts;
  int status;
  int numHEdges;
  int *edgeIdx;
  int num_my_lines, preprocessed;

  DEBUG_TRACE_START(Proc, yo);

  if (pio_info->init_dist_procs < 0){
    pio_info->init_dist_procs = Num_Proc;
  }

  num_my_lines = -1;
  preprocessed = 0;
  sprintf(filename, "%s.mtxp", pio_info->pexo_fname);
  if (pio_info->file_comp == GZIP)
    sprintf(filename, "%s.gz", filename);      /* but we don't uncompress?? TODO */

  if (pio_info->chunk_reader == 1 &&              /* read large file in chunks */
      pio_info->init_dist_type == INITIAL_OWNER)  /* each process gets its own objects */
  {
    /* Each process reads in the mtxp file and keeps only the parts that it owns.
     * Buffer is null terminated.
     * It's possible a process only has mtxp header and owns no vertices, etc.
     */                                     

    fsize = read_mtxp_lines(filename, Proc, Num_Proc, &filebuf, &num_my_lines);

    if (fsize == 0){
      Gen_Error(0, "fatal: insufficient memory or invalid mtxp file");  /* TODO better message */
      return 0;
    }

    preprocessed = 1;
  }
  else if (pio_info->chunk_reader == 0)
  {
    /* Process 0 reads the file and broadcasts it */
    if (Proc == 0) {
      fsize = 0;
  
      rc = stat(filename, &statbuf);
  
      if (rc == 0){
        fsize = statbuf.st_size;
        fp = fopen(filename, "r");
  
        if (!fp){
  	  fsize = 0;
        }
        else{
  	filebuf = (char *)malloc(fsize+1);
  
  	rc = fread(filebuf, 1, fsize, fp);
  
  	if (rc != fsize){
  	  free(filebuf);
  	  fsize = 0;
  	  fp = NULL;
  	}
  	else{
  	  filebuf[fsize] = 0;
  	  fsize++;
  	}
  	fclose(fp);
        }
      }
    }
  
    MPI_Bcast(&fsize, 1, MPI_INT, 0, zoltan_get_global_comm());
  
    if (fsize == 0) {
      sprintf(cmesg, "fatal:  Could not open/read hypergraph file %s", filename);
      Gen_Error(0, cmesg);
      return 0;
    }
  
    if (Proc > 0){
      filebuf = (char *)malloc(fsize);
    }
  
    MPI_Bcast(filebuf, fsize, MPI_BYTE, 0, zoltan_get_global_comm());
  }
  else{
    /* ERROR - we don't handle the zdrive.inp file request */
    Gen_Error(0, "fatal: read_mtxplus_file can not handle this request");
    return 0;
  }

  rc = process_mtxp_file(pio_info, filebuf, fsize, num_my_lines, Num_Proc, Proc,
	  &nGlobalEdges, &nGlobalVtxs, &vtxWDim, &edgeWDim,
          &nMyPins, &myPinI, &myPinJ,
  	  &nMyVtx, &myVtxNum, &myVtxWgts,
  	  &nMyEdgeWgts, &myEWGno, &myEdgeWgts, preprocessed);

  free(filebuf);
  
  MPI_Allreduce(&rc, &status, 1, MPI_INT, MPI_SUM, zoltan_get_global_comm());
  
  if (status != Num_Proc){
    Gen_Error(0, "fatal: invalid mtxp file");  /* TODO better message */
    return 0;
  }

#if PRINT_DEBUG_INFO
  debug_pins(Proc, Num_Proc, nGlobalEdges, nGlobalVtxs, vtxWDim, edgeWDim,
          nMyPins, myPinI, myPinJ,
  	  nMyVtx, myVtxNum, myVtxWgts,
  	  nMyEdgeWgts, myEWGno, myEdgeWgts);

#endif

  /*
   * From the lists of pins, create edge lists.  (Unless
   * the initial pin distribution is by column, in which
   * case we will test the hypergraph query interface's
   * ability to accept pins by column rather than row.)
   *
   * initial pins = row     each process gets full rows (hedges) initially
   *                        and format is ZOLTAN_COMPRESSED_EDGE
   *
   * initial pins = col     each process gets full columns (vertices) initially 
   *                        and format is ZOLTAN_COMPRESSED_VERTEX
   */

  if (pio_info->init_dist_pins != INITIAL_COL){  
    rc = create_edge_lists(nMyPins, myPinI, myPinJ, &numHEdges, &edgeGno, &edgeIdx, &pinGno);
    mesh->format = ZOLTAN_COMPRESSED_EDGE;
  }
  else{                                       
    /* actually creating vertex lists, since we switched
     * the role of I and J in the argument list.
     */
    rc = create_edge_lists(nMyPins, myPinJ, myPinI,
	    &numHEdges, &edgeGno, &edgeIdx, &pinGno);
    mesh->format = ZOLTAN_COMPRESSED_VERTEX;
  }

  MPI_Allreduce(&rc, &status, 1, MPI_INT, MPI_SUM, zoltan_get_global_comm());

  if (status != Num_Proc){
    return 0;
  }

  safe_free((void **)(void *)&myPinI);
  safe_free((void **)(void *)&myPinJ);

  /* Initialize mesh structure for Hypergraph. */
  mesh->data_type = ZOLTAN_HYPERGRAPH;
  mesh->num_elems = nMyVtx;
  mesh->vwgt_dim = vtxWDim;
  mesh->ewgt_dim = 0;
  mesh->elem_array_len = mesh->num_elems + 5;
  mesh->num_dims = 0;
  mesh->num_el_blks = 1;

  mesh->gnhedges = (ZOLTAN_ID_TYPE)nGlobalEdges;
  mesh->nhedges = numHEdges;     /* (or num vertices if CCS) */
  mesh->hewgt_dim = edgeWDim;

  mesh->hvertex_proc = NULL;     /* don't know don't care */
  mesh->hindex = edgeIdx;        /* (or vertex index if CCS) */
  mesh->heNumWgts = nMyEdgeWgts;
  mesh->hewgts = myEdgeWgts;

  if (numHEdges){
    mesh->hgid = (ZOLTAN_ID_TYPE *)malloc(numHEdges * sizeof(ZOLTAN_ID_TYPE));
    for (i=0; i < numHEdges; i++){
      mesh->hgid[i] = (ZOLTAN_ID_TYPE)edgeGno[i];
    }
    free(edgeGno);

    if (edgeIdx[numHEdges]){
      mesh->hvertex = (ZOLTAN_ID_TYPE *)malloc(edgeIdx[numHEdges] * sizeof(ZOLTAN_ID_TYPE));
      for (i=0; i < edgeIdx[numHEdges]; i++){
        mesh->hvertex[i] = (ZOLTAN_ID_TYPE)pinGno[i];
      }
      free(pinGno);
    }
  }

  if (nMyEdgeWgts){
    mesh->heWgtId = (ZOLTAN_ID_TYPE *)malloc(nMyEdgeWgts * sizeof(ZOLTAN_ID_TYPE));
    for (i=0; i < nMyEdgeWgts; i++){
      mesh->heWgtId[i] = (ZOLTAN_ID_TYPE)myEWGno[i];
    }
    free(myEWGno);
  }

  mesh->eb_etypes = (int *) malloc (4 * mesh->num_el_blks * sizeof(int));
  if (!mesh->eb_etypes) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  mesh->eb_ids = mesh->eb_etypes + mesh->num_el_blks;
  mesh->eb_nnodes = mesh->eb_ids + mesh->num_el_blks;
  mesh->eb_nattrs = mesh->eb_nnodes + mesh->num_el_blks;

  mesh->eb_cnts = (ZOLTAN_ID_TYPE *) malloc (mesh->num_el_blks * sizeof(ZOLTAN_ID_TYPE));
  if (!mesh->eb_cnts) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  mesh->eb_names = (char **) malloc (mesh->num_el_blks * sizeof(char *));
  if (!mesh->eb_names) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  mesh->eb_etypes[0] = -1;
  mesh->eb_ids[0] = 1;
  mesh->eb_cnts[0] = nGlobalVtxs;
  mesh->eb_nattrs[0] = 0;
  mesh->eb_nnodes[0] = 0;

  /* allocate space for name */
  mesh->eb_names[0] = (char *) malloc((MAX_STR_LENGTH+1) * sizeof(char));
  if (!mesh->eb_names[0]) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }
  strcpy(mesh->eb_names[0], "hypergraph");

  /* allocate the element structure array */
  mesh->elements = (ELEM_INFO_PTR) malloc (mesh->elem_array_len
					 * sizeof(ELEM_INFO));
  if (!(mesh->elements)) {
    Gen_Error(0, "fatal: insufficient memory");
    return 0;
  }

  /*
   * Write the element structure with the vertices and weights
   */
  for (i = 0; i < mesh->elem_array_len; i++) {
    initialize_element(&(mesh->elements[i]));
    if (i < mesh->num_elems){
      mesh->elements[i].globalID = (ZOLTAN_ID_TYPE)myVtxNum[i];
      mesh->elements[i].my_part  = Proc;
      for (j=0; j<vtxWDim; j++){
	mesh->elements[i].cpu_wgt[j] = myVtxWgts[i*vtxWDim + j];
      }
    }
  }

  safe_free((void **)(void *) &myVtxWgts);
  safe_free((void **)(void *) &myVtxNum);

 if (Debug_Driver > 3)
   print_distributed_mesh(Proc, Num_Proc, mesh);

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mm_cleanup(MESH_INFO_PTR mesh)
{
/* Matrix-market files are one-based, but when we read them, we converted
 * them to zero-based for convenience of indexing arrays.
 * However, we should give Zoltan and the output routines IDs 
 * that have the same base as the input (not converted from one-based 
 * to zero-based).  So if the input was from Matrix-Market files, 
 * add one to all GIDs.
 */
  
  int i, j, sum;
  int nhedges = mesh->nhedges;
  int npins = mesh->hindex[nhedges];

  ELEM_INFO_PTR elements = mesh->elements;
 
  for (i = 0; i < mesh->num_elems; i++) {
    ELEM_INFO_PTR current_elem = &(elements[i]);
    current_elem->globalID++;
    for (j = 0; j < current_elem->adj_len; j++) {
      if (current_elem->adj[j] == ZOLTAN_ID_INVALID) continue;
      if (current_elem->adj_proc[j] != mesh->proc) {
        current_elem->adj[j]++;
      }
    }
  }
  for (i = 0; i < nhedges; i++) mesh->hgid[i]++;
  for (i = 0; i < npins; i++) mesh->hvertex[i]++;
  sum = 0;
  for (i = 0; i < mesh->necmap; i++) sum += mesh->ecmap_cnt[i];
  for (i = 0; i < sum; i++) mesh->ecmap_neighids[i]++;
}

/*****************************************************************************/
static int dist_hyperedges(
  MPI_Comm comm,		/* MPI Communicator */
  PARIO_INFO_PTR pio_info,      /* Parallel IO info */
  int     host_proc,		/* processor where all the data is initially */
  int     base,                 /* indicates whether input is 0-based or
				   1-based (i.e., is lowest vertex number
				   0 or 1?). */
  int     gnvtxs,               /* global number of vertices */
  int     *gnhedges,		/* global number of hyperedges */
  int     *nhedges,		/* local number of hyperedges */
  int     **hgid,		/* global hyperedge numbers */
  int     **hindex,		/* Starting hvertex entry for hyperedges */
  int     **hvertex,		/* Array of vertices in hyperedges; returned
				   values are global IDs for vertices */
  int     **hvertex_proc,	/* Array of processor assignments for
				   vertices in hvertex.  */
  int     *hewgt_dim,           /* number of weights per hyperedge */
  float   **hewgts,		/* hyperedge weight list data */
  short   *assignments          /* assignments from Chaco file; may be NULL */
)
{
/*
 * Distribute hyperedges from one processor to all processors.
 * Vertex distribution is assumed already done through chaco_dist_graph.
 * The memory for the hyperedges on the host node is freed
 * and fresh memory is allocated for the distributed hyperedges.
 */

const char *yo = "dist_hyperedges";
int nprocs, myproc, i, h, p = 0;
int *old_hindex = NULL, *old_hvertex = NULL, *old_hvertex_proc = NULL;
int *size = NULL, *num_send = NULL;
int *send = NULL;
int *send_hgid = NULL, *send_hindex = NULL;
int *send_hvertex = NULL, *send_hvertex_proc = NULL;
int max_size, max_num_send;
int hcnt[2];
int hecnt, hvcnt;
float *old_hewgts = NULL;
float *send_hewgts = NULL;
MPI_Status status;
int num_dist_procs;
int hedge_init_dist_type;

  hedge_init_dist_type = (pio_info->init_dist_type != INITIAL_FILE
			  ? pio_info->init_dist_type
			  : INITIAL_LINEAR);

  /* Determine number of processors and my rank. */
  MPI_Comm_size (comm, &nprocs );
  MPI_Comm_rank (comm, &myproc );

  DEBUG_TRACE_START(myproc, yo);

  if (pio_info->init_dist_procs > 0 && pio_info->init_dist_procs <= nprocs)
    num_dist_procs = pio_info->init_dist_procs;
  else
    /* Reset num_dist_procs if not set validly by input */
    num_dist_procs = nprocs;

  /* Broadcast to all procs */
  MPI_Bcast( hewgt_dim, 1, MPI_INT, host_proc, comm);
  MPI_Bcast( nhedges, 1, MPI_INT, host_proc, comm);
  *gnhedges = *nhedges;

  /* Initialize */
  if (*gnhedges == 0) {
    *hindex = (int *) malloc(sizeof(int));
    if (!(*hindex)) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
    (*hindex)[0] = 0;
    *hvertex = NULL;
    *hvertex_proc = NULL;
    return 1;
  }

  if (nprocs == 1) {
    *nhedges = *gnhedges;
    *hgid = (int *) malloc(*gnhedges * sizeof(int));
    *hvertex_proc = (int *) malloc((*hindex)[*gnhedges] * sizeof(int));
    if ((*gnhedges && !(*hgid)) || ((*hindex)[*gnhedges] && !(*hvertex_proc))) {
      Gen_Error(0, "fatal: insufficient memory");
      return 0;
    }
    for (h = 0; h < *gnhedges; h++)
      (*hgid)[h] = h + base;           /* Want the same numbering than for vertices */
    for (h = 0; h < (*hindex)[*gnhedges]; h++)
      (*hvertex_proc)[h] = 0;
    return 1;
  }
  if (myproc == host_proc) {
    /* Store pointers to original data */
    old_hindex  = *hindex;
    old_hvertex = *hvertex;
    old_hewgts = *hewgts;
    old_hvertex_proc = (int *) malloc(old_hindex[*gnhedges] * sizeof(int));

    /* Allocate space for size and send flags */
    size = (int *) calloc(2 * nprocs, sizeof(int));
    num_send = size + nprocs;
    send = (int *) malloc(*gnhedges * sizeof(int *));
    if ((old_hindex[*gnhedges] && !old_hvertex_proc) || !size ||
        (*gnhedges && !send)) {
      Gen_Error(0, "fatal: memory error");
      return 0;
    }

    /* Determine to which processors hyperedges should be sent */
    for (h = 0; h < *gnhedges; h++) {
      if (hedge_init_dist_type == INITIAL_CYCLIC)
	p = h % num_dist_procs;
      else if (hedge_init_dist_type == INITIAL_LINEAR)
	p = (int) ((float) h * (float) num_dist_procs / (float)(*gnhedges));
      else if (hedge_init_dist_type == INITIAL_OWNER)
	p = ch_dist_proc(old_hvertex[old_hindex[h]], assignments, base);
      size[p] += (old_hindex[h+1] - old_hindex[h]);
      num_send[p]++;
      send[h] = p;
      for (i = old_hindex[h]; i < old_hindex[h+1]; i++)
	old_hvertex_proc[i] = ch_dist_proc(old_hvertex[i], assignments, base);
    }

    /* Determine size of send buffers and allocate them. */
    max_size = 0;
    max_num_send = 0;
    for (p = 0; p < nprocs; p++) {
      if (size[p] > max_size) max_size = size[p];
      if (num_send[p] > max_num_send) max_num_send = num_send[p];
    }

    send_hgid = (int *) malloc((max_num_send) * sizeof(int));
    send_hindex = (int *) malloc((max_num_send+1) * sizeof(int));
    send_hvertex = (int *) malloc(max_size * sizeof(int));
    send_hvertex_proc = (int *) malloc(max_size * sizeof(int));
    if (*hewgt_dim)
      send_hewgts = (float *) malloc(max_num_send*(*hewgt_dim)*sizeof(float));
    if ((max_num_send && !send_hgid) || !send_hindex ||
	(max_size && (!send_hvertex || !send_hvertex_proc)) ||
	(max_num_send && *hewgt_dim && !send_hewgts)) {
      Gen_Error(0, "fatal: memory error in dist_hyperedges");
      return 0;
    }

    /* Load and send data */
    for (p = 0; p < nprocs; p++) {

      if (p == myproc) continue;

      hecnt = 0;
      hvcnt = 0;
      send_hindex[0] = 0;
      for (h = 0; h < *gnhedges; h++) {
	if (send[h]==p) {
	  send_hgid[hecnt] = h + base; /* Want the same numbering than for vertices */
	  send_hindex[hecnt+1] = send_hindex[hecnt]
			       + (old_hindex[h+1] - old_hindex[h]);
	  for (i = 0; i < *hewgt_dim; i++)
	    send_hewgts[hecnt*(*hewgt_dim)+i] = old_hewgts[h*(*hewgt_dim)+i];
	  hecnt++;
	  for (i = old_hindex[h]; i < old_hindex[h+1]; i++) {
	    send_hvertex[hvcnt] = old_hvertex[i];
	    send_hvertex_proc[hvcnt] = old_hvertex_proc[i];
	    hvcnt++;
	  }
	}
      }

      hcnt[0] = hecnt;
      hcnt[1] = hvcnt;

      MPI_Send(hcnt, 2, MPI_INT, p, 1, comm);
      MPI_Send(send_hgid, hecnt, MPI_INT, p, 2, comm);
      MPI_Send(send_hindex, hecnt+1, MPI_INT, p, 3, comm);
      MPI_Send(send_hvertex, hvcnt, MPI_INT, p, 4, comm);
      MPI_Send(send_hvertex_proc, hvcnt, MPI_INT, p, 5, comm);
      if (*hewgt_dim)
	MPI_Send(send_hewgts, hecnt*(*hewgt_dim), MPI_FLOAT, p, 6, comm);
    }

    safe_free((void **)(void *) &send_hgid);
    safe_free((void **)(void *) &send_hindex);
    safe_free((void **)(void *) &send_hvertex);
    safe_free((void **)(void *) &send_hvertex_proc);
    safe_free((void **)(void *) &send_hewgts);

    /* Copy data owned by myproc into new local storage */
    *nhedges = num_send[myproc];
    *hgid = (int *) malloc(*nhedges * sizeof(int));
    *hindex = (int *) malloc((*nhedges+1) * sizeof(int));
    *hvertex = (int *) malloc(size[myproc] * sizeof(int));
    *hvertex_proc = (int *) malloc(size[myproc] * sizeof(int));
    if (*hewgt_dim)
      *hewgts = (float *) malloc(*nhedges * *hewgt_dim * sizeof(float));
    if ((*nhedges && !(*hgid)) || !(*hindex) ||
	(size[myproc] && (!(*hvertex) || !(*hvertex_proc))) ||
	(*nhedges && *hewgt_dim && !(*hewgts))) {
      Gen_Error(0, "fatal: memory error in dist_hyperedges");
      return 0;
    }

    hecnt = 0;
    hvcnt = 0;
    (*hindex)[0] = 0;

    for (h = 0; h < *gnhedges; h++) {
      if (send[h]==myproc) {
	(*hgid)[hecnt] = h + base;  /* Want the same numbering than for vertices */
	(*hindex)[hecnt+1] = (*hindex)[hecnt]
			   + (old_hindex[h+1] - old_hindex[h]);
	for (i = 0; i < *hewgt_dim; i++)
	  (*hewgts)[hecnt*(*hewgt_dim)+i] = old_hewgts[h*(*hewgt_dim)+i];
	hecnt++;
	for (i = old_hindex[h]; i < old_hindex[h+1]; i++) {
	  (*hvertex_proc)[hvcnt] = old_hvertex_proc[i];
	  (*hvertex)[hvcnt] = old_hvertex[i];  /* Global index */
	  hvcnt++;
	}
      }
    }

    safe_free((void **)(void *) &send);
    safe_free((void **)(void *) &size);
    safe_free((void **)(void *) &old_hindex);
    safe_free((void **)(void *) &old_hvertex);
    safe_free((void **)(void *) &old_hvertex_proc);
    safe_free((void **)(void *) &old_hewgts);
  }
  else {
    /* host_proc != myproc; receive hedge data from host_proc */
    MPI_Recv(hcnt, 2, MPI_INT, host_proc, 1, comm, &status);
    *nhedges = hcnt[0];
    *hgid = (int *) malloc(hcnt[0] * sizeof(int));
    *hindex = (int *) malloc((hcnt[0]+1) * sizeof(int));
    *hvertex = (int *) malloc(hcnt[1] * sizeof(int));
    *hvertex_proc = (int *) malloc(hcnt[1] * sizeof(int));
    if (*hewgt_dim)
      *hewgts = (float *) malloc(hcnt[0] * *hewgt_dim * sizeof(float));
    if ((hcnt[0] && !(*hgid)) || !(*hindex) ||
	(hcnt[1] && (!(*hvertex) || !(*hvertex_proc))) ||
	(hcnt[0] && *hewgt_dim && !(*hewgts))) {
      Gen_Error(0, "fatal: memory error in dist_hyperedges");
      return 0;
    }
    MPI_Recv(*hgid, hcnt[0], MPI_INT, host_proc, 2, comm, &status);
    MPI_Recv(*hindex, hcnt[0]+1, MPI_INT, host_proc, 3, comm, &status);
    MPI_Recv(*hvertex, hcnt[1], MPI_INT, host_proc, 4, comm, &status);
    MPI_Recv(*hvertex_proc, hcnt[1], MPI_INT, host_proc, 5, comm, &status);
    if (*hewgt_dim)
      MPI_Recv(*hewgts, hcnt[0]* *hewgt_dim, MPI_FLOAT, host_proc, 6, comm,
	       &status);
  }

  DEBUG_TRACE_END(myproc, yo);
  return 1;
}
/*****************************************************************************/
static int create_edge_lists(int nMyPins, ZOLTAN_ID_TYPE *myPinI, ZOLTAN_ID_TYPE *myPinJ,
      int *numHEdges, ZOLTAN_ID_TYPE **edgeGno, int **edgeIdx, ZOLTAN_ID_TYPE **pinGno)
{
int nedges, i, lid;
int *count, *start, *idx;
ZOLTAN_ID_TYPE *eidList, *pins, *match;

  *numHEdges = 0;
  *edgeGno = *pinGno = NULL;
  *edgeIdx = NULL;

  if (nMyPins == 0){
    *edgeIdx = (int *)malloc(sizeof(int));
    **edgeIdx = 0; /* number of pins */
    return 1;
  }

  eidList = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nMyPins);
  pins = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nMyPins);
  idx = (int *)malloc(sizeof(int) * nMyPins);

  if (!eidList || !idx || !pins){
    safe_free((void **)(void *)&eidList);
    safe_free((void **)(void *)&idx);
    safe_free((void **)(void *)&pins);
    Gen_Error(0, "memory allocation");
    return 0;
  }

  /* create a sorted list of unique edges IDs */

  memcpy(eidList, myPinI, sizeof(ZOLTAN_ID_TYPE) * nMyPins);
  qsort((void *)eidList, nMyPins, sizeof(ZOLTAN_ID_TYPE), _zoltan_sortFunc);

  for (i=1, nedges=1; i < nMyPins; i++){
    if (eidList[i] == eidList[nedges-1]) continue;
    if (nedges < i){
      eidList[nedges] = eidList[i];
    }
    nedges++;
  }

  eidList = (ZOLTAN_ID_TYPE *)realloc(eidList, nedges * sizeof(ZOLTAN_ID_TYPE));

  /* Count pins in each edge, map pins to local edge index */

  count = (int *)calloc(sizeof(int), nedges);

  if (!count){
    safe_free((void **)(void *)&eidList);
    safe_free((void **)(void *)&idx);
    safe_free((void **)(void *)&pins);
    Gen_Error(0, "memory allocation");
    return 0;
  }

  for (i = 0; i<nMyPins; i++){
    match = (ZOLTAN_ID_TYPE *)bsearch((const void *)(myPinI + i), (const void *)eidList,
		    nedges, sizeof(ZOLTAN_ID_TYPE), _zoltan_sortFunc);

    if (!match){
      safe_free((void **)(void *)&eidList);
      safe_free((void **)(void *)&idx);
      safe_free((void **)(void *)&pins);
      safe_free((void **)(void *)&count);
      Gen_Error(0, "bsearch failure");
      return 0;
    }

    lid = (int)(match - eidList);
    count[lid]++;
    idx[i] = lid;
  }

  /* Create edges lists */

  start = (int *)malloc(sizeof(int) * (nedges+1));
  if (!start){
    safe_free((void **)(void *)&eidList);
    safe_free((void **)(void *)&idx);
    safe_free((void **)(void *)&pins);
    safe_free((void **)(void *)&count);
    Gen_Error(0, "memory allocation");
    return 0;
  }
  start[0] = 0;
  for (i=0; i<nedges; i++){
    start[i+1] = start[i] + count[i];
    count[i] = 0;
  }

  for (i=0; i<nMyPins; i++){
    lid = idx[i];
    pins[start[lid] + count[lid]] = myPinJ[i];
    count[lid]++;
  }

  free(count);
  free(idx);

  *numHEdges = nedges;
  *edgeGno   = eidList;
  *edgeIdx   = start;
  *pinGno    = pins;

  return 1;
}

/*****************************************************************************/

static int process_mtxp_file(PARIO_INFO_PTR pio_info,
  char *filebuf, size_t fsize, int num_my_lines,
  int nprocs, int myrank,
  ZOLTAN_ID_TYPE *nGlobalEdges, ZOLTAN_ID_TYPE *nGlobalVtxs, 
  int *vtxWDim, int *edgeWDim,
  int *nMyPins, 
  ZOLTAN_ID_TYPE **myPinI, ZOLTAN_ID_TYPE **myPinJ,
  int *nMyVtx, ZOLTAN_ID_TYPE **myVtxNum, float **myVtxWgts,
  int *nMyEdgeWgts, ZOLTAN_ID_TYPE **myEdgeNum, float **myEdgeWgts,
  int preprocessed)
{
ZOLTAN_ID_TYPE nedges, nvtxs, npins, numew;
ZOLTAN_ID_TYPE eid, vid, i;
ZOLTAN_ID_TYPE myminPin=0, mymaxPin=0, myminVtx=0, mymaxVtx=0, myshare, share;
ZOLTAN_ID_TYPE *myi, *myj, *myvno, *myeno;
int ok, not_ok;
int vdim, edim, nFileProcs;
int countMyVtxs, countMyEdges, countMyPins;
int proc, mine, nextpin, rc, counter;
int nexte, nextv, nDistProcs=0;
float pinVal;
char *line, *token, *pinBuf, *vwgtBuf, *ewgtBuf;
float *myvwgt, *myewgt;
char cmesg[256];
char linestr[MATRIX_MARKET_MAX_LINE+1];

  *nGlobalEdges = *nGlobalVtxs = 0;
  *myPinI = *myPinJ = *myVtxNum = *myEdgeNum = NULL;
  *nMyPins = *nMyVtx = *nMyEdgeWgts = 0;
  *vtxWDim = *edgeWDim = 0;
  *myVtxWgts = *myEdgeWgts = NULL;

  ok = 1;
  not_ok = 0;

  line = first_char(filebuf, fsize);
  if (line && (*line == COMMENT_CHAR)){
    line = next_line_of_data(line, fsize, linestr);  /* skip comments & blanks */
  }

  if (!line){
    sprintf(cmesg, "Truncated file\n");
    Gen_Error(0, cmesg);
    return not_ok;
  }

  rc = sscanf(linestr, ZOLTAN_ID_SPEC ZOLTAN_ID_SPEC ZOLTAN_ID_SPEC " %d %d " ZOLTAN_ID_SPEC "%d",
	    &nedges, &nvtxs, &npins, &nFileProcs, &vdim, &numew, &edim);

  if (rc != 7){
    sprintf(cmesg, "%s\nFirst line should have 7 values in it\n",linestr);
    Gen_Error(0, cmesg);
    return not_ok;
  }

  *nGlobalEdges = nedges;
  *nGlobalVtxs = nvtxs;
  *vtxWDim = vdim;
  *edgeWDim = edim;

  if (preprocessed && (num_my_lines == 0)){
    return ok;
  }

  if (!preprocessed){

    myminPin = mymaxPin = -1;
    myminVtx = mymaxVtx = -1;
    nDistProcs = nprocs;

    if ((pio_info->init_dist_procs > 0) &&
        (pio_info->init_dist_procs < nprocs)){
      nDistProcs = pio_info->init_dist_procs;
    }
    
    if (pio_info->init_dist_type == INITIAL_LINEAR){ /* vertex distribution */
      if (myrank < nDistProcs){
        share = nvtxs / nDistProcs;
        i = nvtxs - (nDistProcs * share);
        myshare = ((myrank < (int)i) ? share+1 : share);
        myminVtx = myrank * myshare;
        if (myrank >= (int)i) myminVtx += i;
        mymaxVtx = myminVtx + myshare - 1;
      }
    }
    if (pio_info->init_dist_pins == INITIAL_LINEAR){ /* pin distribution */
      share = npins / nprocs;
      i = npins - (nprocs * share);
      myshare = ((myrank < (int)i) ? share+1 : share);
      myminPin = myrank * myshare;
      if (myrank >= (int)i) myminPin += i;
      mymaxPin = myminPin + myshare - 1;
    }
  }

  myvno = myeno = myi = myj = NULL;
  myewgt = myvwgt = NULL;
  pinBuf = vwgtBuf = ewgtBuf = NULL;

  /* Read through the pins, vertex weights and edge weights.
   * Accumulate all vertex and edge IDs, and map these to
   * consecutive global numbers beginning with zero.
   *
   * Also count my pins, my vertices, and the number of edges
   * for which I provide weights.
   */

  countMyPins = 0;
  countMyVtxs = 0;
  countMyEdges = 0;
  nexte = nextv = 0;

  section_starts(line, &pinBuf, &vwgtBuf, &ewgtBuf);

  if (preprocessed){

    if (pinBuf){
      line = pinBuf;
      while (line){
        countMyPins++;
        line = next_line(line, fsize);
        if (line && line[0] == COMMENT_CHAR) break;
      }
    }

    if (vwgtBuf){
      line = vwgtBuf;
      while (line){
        countMyVtxs++;
        line = next_line(line, fsize);
        if (line && line[0] == COMMENT_CHAR) break;
      }
    }

    if (ewgtBuf){
      line = ewgtBuf;
      while (line){
        countMyEdges++;
        line = next_line(line, fsize);
        if (line && line[0] == COMMENT_CHAR) break;
      }
    }
  }
  else{             /* file is not preprocessed */

    line = pinBuf;
    counter = 0;
    /* Skip any additional comment lines before pins begin
     * Zoltan_Generate_Files adds an extra comment line here to mtxp files.
     */
    while (line) {
      if (line[0] != COMMENT_CHAR)
        break;
      line = next_line(line, fsize);
    }
     
    while (line){          /* PINS */

      make_string(line, linestr);
      rc = sscanf(linestr, ZOLTAN_ID_SPEC ZOLTAN_ID_SPEC "%f %d", &eid, &vid, &pinVal, &proc);

      if ((rc != 4) || (eid < 1) || (eid > nedges) ||
	  (vid < 1) || (vid > nvtxs) || (proc < 0) || (proc >= nFileProcs)){
        sprintf(cmesg,"%s\nlooking for \"edge vertex pin process\"\n",linestr);
        Gen_Error(0, cmesg);
        goto failure;
      }

      eid -= 1;
      vid -= 1;
      mine = my_pin(eid, vid, proc, counter++, npins,
         myminPin, mymaxPin, myrank, nprocs, pio_info);

      if (mine){
        countMyPins++;
      }

      line = next_line(line, fsize);
      if (line && line[0] == COMMENT_CHAR) break;
    }

    line = vwgtBuf;
    /* Skip any additional comment lines before vwgts begin
     * Zoltan_Generate_Files adds an extra comment line here to mtxp files.
     */
    while (line) {
      if (line[0] != COMMENT_CHAR)
        break;
      line = next_line(line, fsize);
    }

    while(line) {        /* VERTICES and possibly WEIGHTS */
  
      make_string(line, linestr);
      rc = sscanf(linestr, ZOLTAN_ID_SPEC, &vid);

      token = get_nth_token(linestr, vdim + 1, strlen(linestr), 1, (char)0);
      if (token) proc = atoi(token);

      if ((rc != 1) || !token ||
  	 (vid < 1) || (vid > nvtxs) || (proc < 0) || (proc >= nFileProcs)){
        sprintf(cmesg,
        "%s\nlooking for \"vertex {optional weights} process\"\n",linestr);
        Gen_Error(0, cmesg);
        goto failure;
      }
  
      vid -= 1;
      mine = my_vtx(proc, vid, myminVtx, mymaxVtx, myrank, nDistProcs, pio_info);
  
      if (mine){
        countMyVtxs++;
      }
      line = next_line(line, fsize);
      if (line && line[0] == COMMENT_CHAR) break;
    }

    if (numew > 0){                      /* HYPEREDGE WEIGHTS */
      line = ewgtBuf;
      /* Skip any additional comment lines before ewgts begin
       * Zoltan_Generate_Files adds an extra comment line here to mtxp files.
       */
      while (line) {
        if (line[0] != COMMENT_CHAR)
          break;
        line = next_line(line, fsize);
      }

      while(line) {   
        make_string(line, linestr);
        rc = sscanf(linestr, ZOLTAN_ID_SPEC, &eid);
        token = get_nth_token(linestr, edim + 1, strlen(linestr), 1, (char)0);

        if (token) proc = atoi(token);
        if ((rc != 1) || !token ||
  	  (eid < 1) || (eid > nedges) || (proc < 0) || (proc >= nFileProcs)){
  	  sprintf(cmesg,
  	  "%s\nlooking for \"edge {optional weights} process\"\n",linestr);
  	  Gen_Error(0, cmesg);
  	  goto failure;
        }
        proc = atoi(token);
  
        if ((proc % nprocs) == myrank){
  	  countMyEdges++;
        }

        line = next_line(line, fsize);
        if (line && line[0] == COMMENT_CHAR) break;
      }
    }
  } 

  /* Start over at beginning of file and save my pins, and weights */

  mine = 1;

  if (countMyPins > 0){
    myi = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * countMyPins);
    myj = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * countMyPins);
    if (!myi || !myj){
      sprintf(cmesg,"memory allocation\n");
      Gen_Error(0, cmesg);
      goto failure;
    }
    nextpin = 0;

    make_string(pinBuf, linestr);
    line = pinBuf;
    counter = 0;

    while(line){

      sscanf(linestr, ZOLTAN_ID_SPEC ZOLTAN_ID_SPEC "%f %d", &eid, &vid, &pinVal, &proc);

      eid -= 1;
      vid -= 1;

      if (!preprocessed)
        mine = my_pin(eid, vid, proc, counter++, npins,
		  myminPin, mymaxPin, myrank, nprocs, pio_info);

      if (mine){
	myi[nextpin] = eid;
	myj[nextpin] = vid;
	nextpin++;
      }

      if (nextpin == countMyPins) break;
      line = next_line_of_data(line, fsize, linestr);
    }
  }

  if (countMyVtxs){
    myvno = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * countMyVtxs);
    if (vdim > 0){
      myvwgt = (float *)malloc(sizeof(float) * countMyVtxs * vdim);
    }
    if (!myvno || (vdim && !myvwgt)){
      sprintf(cmesg,"memory allocation\n");
      Gen_Error(0, cmesg);
      goto failure;
    }
    nextv = 0;

    make_string(vwgtBuf, linestr);
    line = vwgtBuf;
    counter = 0;

    while(line){

      sscanf(linestr, ZOLTAN_ID_SPEC, &vid);
      vid -= 1;

      if (!preprocessed){
        token = get_nth_token(linestr, vdim + 1, strlen(linestr), 1, (char)0);
        proc = atoi(token);
        mine = my_vtx(proc, vid, myminVtx, mymaxVtx, myrank, nDistProcs, pio_info);
      }

      if (mine){
        int jj;
	myvno[nextv] = vid;
	for (jj=0; jj<vdim; jj++){
	  token = get_nth_token(linestr, 1 + jj, strlen(linestr), 1, (char)0);
	  if (!token){
	    sprintf(cmesg,"%s\nCan't find %d vertex weights\n",linestr,vdim);
	    Gen_Error(0, cmesg);
	    goto failure;
	  }
	  myvwgt[nextv*vdim + jj] = (float)atof(token);
	}
	nextv++;
      }
      if (nextv == countMyVtxs) break;

      line = next_line_of_data(line, fsize, linestr);
    }
  }

  if (countMyEdges > 0){
    myeno = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * countMyEdges);
    myewgt = (float *)malloc(sizeof(float) * countMyEdges * edim);
    if (!myeno || !myewgt){
      sprintf(cmesg,"memory allocation\n");
      Gen_Error(0, cmesg);
      goto failure;
    }
    nexte = 0;

    make_string(ewgtBuf, linestr);
    line = ewgtBuf;

    while (line){
      sscanf(linestr, ZOLTAN_ID_SPEC, &eid);
      eid -= 1;

      if (!preprocessed){
        mine = 0;
        token = get_nth_token(linestr, edim + 1, strlen(linestr), 1, (char)0);
        proc = atoi(token);
        if ((proc % nprocs) == myrank)
          mine = 1;
      }

      if (mine){
        int jj;
	myeno[nexte] = eid;
	for (jj=0; jj<edim; jj++){
	  token = get_nth_token(linestr, 1 + jj, strlen(linestr), 1, (char)0);
	  if (!token){
	    sprintf(cmesg,"%s\nCan't find %d edge weights\n",linestr,edim);
	    Gen_Error(0, cmesg);
	    goto failure;
	  }
	  myewgt[nexte*edim + jj] = (float)atof(token);
	}
	nexte++;
      }
      if (nexte == countMyEdges) break;
      line = next_line_of_data(line, fsize, linestr);
    }
  }

  rc = 1;   /* success */
  goto done;

failure:
  if (myvno) free(myvno);
  if (myvwgt) free(myvwgt);
  if (myeno) free(myeno);
  if (myewgt) free(myewgt);
  if (myi) free(myi);
  if (myj) free(myj);
  nedges = nvtxs = vdim = edim = 0;
  countMyPins = countMyVtxs = countMyEdges = 0;
  rc = 0;

done:

  *nGlobalEdges = nedges;
  *nGlobalVtxs = nvtxs;
  *vtxWDim = vdim;
  *edgeWDim = edim;

  *nMyPins = countMyPins;
  *myPinI  = myi;
  *myPinJ  = myj;

  *nMyVtx = countMyVtxs;
  *myVtxNum = myvno;
  *myVtxWgts = myvwgt;

  *nMyEdgeWgts = countMyEdges;
  *myEdgeNum = myeno;
  *myEdgeWgts = myewgt;

  return rc;
}

/* Read the mtxp file and return only comments plus lines I "own",
 * null-terminate the buffer, and set the number of objects I own. 
 * Return the size of the buffer, 0 on failure.
 */

int read_mtxp_lines(char *fname, int my_rank, int nprocs, char **mine, int *num)
{
off_t fsize=0, bytes_read;
struct stat statbuf;
size_t inbufsize, outbufsize, nbytes;
char *inbuf=NULL, *outbuf=NULL, *c;
char *c_in, *c_out, *c_end;
char *c_in_end;
char *buf;
long count_mine;
int rc, max_sanity, total_line_found, line_size=0;
ZOLTAN_ID_TYPE nedges, nvtxs, npins, numew; 
int numOwners, vdim, edim;
int i=0;
int owning_proc;
int status;
int *new_owner=NULL;
FILE *fp=NULL;

  *mine = NULL;
  *num = 0;

  rc = stat(fname, &statbuf);

  if (rc){
    status = 0;
    goto End;
  }
  
  fsize = statbuf.st_size;
  fp = fopen(fname, "r");
  if (!fp){
    status = 0;
    goto End;
  }

  outbufsize = fsize / (nprocs - 1);
  outbuf = (char *)malloc(outbufsize+1);

  if (!outbuf){
    status = 0;
    goto End;
  }

#if PRINT_DEBUG_INFO
  if (fsize > 500){
    inbufsize = 500;
    inbuf = malloc(inbufsize+1);

    if (!inbuf){
      inbufsize *= .75;
      inbuf = malloc(inbufsize+1);
  
      if (!inbuf){
        inbufsize *= .5;
        inbuf = malloc(inbufsize+1);
        if (!inbuf){
          status = 0;
          goto End;
        }
      }
    }
  }
#else
  if (fsize > 10*1024){

    /* TODO a Zoltan utility that uses hwloc info to decide a reasonable amount of
     *     memory for Zoltan to use.  Maybe more of a TO_CONSIDER.
     */
    
    inbufsize = (fsize > 100*1024 ? 100*1024 : fsize);
    inbuf = malloc(inbufsize+1);

    if (!inbuf){
      inbufsize *= .75;
      inbuf = malloc(inbufsize+1);
  
      if (!inbuf){
        inbufsize *= .5;
        inbuf = malloc(inbufsize+1);
        if (!inbuf){
          status = 0;
          goto End;
        }
      }
    }
  }
#endif
  else{
    inbufsize = fsize;
    inbuf = malloc(inbufsize);
    if (!inbuf){
      status = 0;
      goto End;
    }
  }

  c_out = outbuf;
  max_sanity = MATRIX_MARKET_MAX_LINE+1;
  total_line_found=0;   /* totals at the top of file */
  bytes_read = 0;
  count_mine = 0;
  nbytes = 0;

  while (1){
    /* read next chunk */
    size_t rrc;
    rrc = fread(inbuf, 1, inbufsize, fp);
    if (rrc < inbufsize){
      if ((rrc == 0) && feof(fp)) break;
      if (ferror(fp)){
        status = 0;
        goto End;
      }
    }

    bytes_read += rrc;

    if ((bytes_read == fsize) && (inbuf[rrc-1] != '\n')){
      /* we assume last character in the file is new line */
      inbuf[rrc++] = '\n';
    }

    c_in = inbuf;
    c_in_end = inbuf + rrc;

    while (c_in < c_in_end){
      c_in = first_char(c_in, max_sanity);  /* skip blanks */

      c_end = c_in;
      line_size = 0;

      while ((c_end < c_in_end) && (*c_end != '\n') && (line_size < max_sanity)){
        c_end++;
        line_size++;
      }

      if (line_size == max_sanity){
        status = 0;
        goto End;
      }

      if (c_end == c_in_end){         /* file ends in mid-line */
        line_size = c_in_end - c_in;
        break;
      }

      *c_end = 0;                /* null terminate the line */

      line_size++;               /* plus '\n' */

      if (line_size > 1){

        if (!total_line_found){
          if (isdigit(c_in[0])){
            total_line_found=1;
            rc = sscanf(c_in, ZOLTAN_ID_SPEC ZOLTAN_ID_SPEC ZOLTAN_ID_SPEC " %d %d " ZOLTAN_ID_SPEC "%d",
	      &nedges, &nvtxs, &npins, &numOwners, &vdim, &numew, &edim);
            if (rc != 7){
              status = 0;
              goto End;
            }
            if (numOwners > nprocs){
              new_owner = (int *)malloc(sizeof(int) * numOwners);
              for (i=0; i < numOwners; i++){
                new_owner[i] = i % nprocs;
              }
            }
          }
        }
        else if (c_in[0] != COMMENT_CHAR){
          if (!total_line_found){
            status = 0;
            goto End;
          }
          c = get_nth_token(c_in, 0, max_sanity, -1, (char)0);
          if (c){
            owning_proc = atoi(c);
          }
          else{
            status = 0;
            goto End;
          }

          if ((owning_proc == my_rank) ||
              (new_owner && new_owner[owning_proc] == my_rank)){
            count_mine++;
          }
          else{
            line_size = 0;
          }
        }

        if (line_size > 0){

          /* Line is comment, totals at top, or one of mine */

          if (nbytes + line_size > outbufsize){
              outbufsize += line_size;
              outbufsize *= 1.2;
  
              buf = (char *)malloc(outbufsize+1);
              if (!buf){
                status = 0;
                goto End;
              }
              memcpy(buf, outbuf, nbytes);
              free(outbuf);
              c_out = buf + nbytes;
              outbuf = buf;
          }
  
          sprintf(c_out, "%s\n", c_in);
          c_out += line_size;
          nbytes += line_size;
        }
      }

      c_in = c_end+1;
    }

    if ((c_in < c_in_end) && (line_size > 0)){
      /* rewind fp to start of incomplete line */
      fseek(fp, -line_size, SEEK_CUR);
    }
  }

  *c_out++ = 0;   /* null terminate buffer */

  status = c_out - outbuf;

  *mine = outbuf;
  *num = count_mine;

End:
  if (fp) fclose(fp);
  if (new_owner) free(new_owner);
  if (inbuf) free(inbuf);
  if (outbuf && (*mine != outbuf)) free(outbuf);
  return status;
}

/*****************************************************************************
   functions to help read the .mtxp file

******************************************************************************/

static char *get_nth_token(char *line, int nth, /* get nth token (0-based) */
                         size_t max, 
                         int direction,      /* 1: from the front, -1: from the back */
                         char terminator)    /* null or \n terminated */
{
char *c1, *c2;
char *loc[10];    /* assume nth is at most 9 */
size_t l;
int found;

  if (nth > 9){
    return NULL;
    /* TODO - error */
  }

  if ((c2 = c1 = first_char(line, max)) == NULL) return NULL;

  l = (size_t)(c1-line);

  while ((*c2 != terminator) && (l < max)){ c2++; l++;} 

  if ((c2 == c1) ||            /* no tokens */
      (*c2 != terminator))     /* invalid line */
    return NULL;

  if (c2 == c1 + 1){  /* only one token in line */
    if (nth == 0)
      return c1;
    else
      return NULL;
  }

  c2--;

  while ((c2 > c1) && IS_BLANK(*c2)) c2--; 
  *(c2+1) = 0;

  found = 0;     /* null terminate last token */

  while (*c1){
    loc[found++] = c1++;
    while (*c1 && !IS_BLANK(*c1)) c1++;
    while (*c1 && IS_BLANK(*c1)) c1++;
  }

  if (nth < found){
    if (direction == 1)
      c1 = loc[nth];
    else
      c1 = loc[found - nth - 1];
  }
  else
    c1 = NULL;
  
  return c1;
}

static char *first_char(char *buf, size_t max)
{
char *c = buf;
size_t sanity = 0;

  while (*c && isspace(*c) && (sanity++ < max)) c++;

  if ((sanity > max) || !*c) return NULL;

  return c;
}

static void make_string(char *buf, char *str)
{
  *str = 0;

  /* replace new line character with NULL character */

  while (*buf && (*buf != '\n')){
    *str++ = *buf++;
  }

  *str = 0;
}

static void section_starts(char *buf, char **pins, char **vwgts, char **ewgts)
{
char *c, line[MATRIX_MARKET_MAX_LINE+1];
char *pin_string="Edge and Vertex IDs";  /* These strings are literal             */
char *vwgt_string="Vertex weights";      /* text from mtxp file.  If file changes */
char *ewgt_string="Edge weights";        /* this text has to change too.          */

  *pins = NULL;
  *vwgts = NULL;
  *ewgts = NULL;

  if (!buf) return;

  /* Find pointers to start of pin section, vertex weight section, edge weight section.
   * It's not an error if one or more of these are missing.
   */

  c = buf;

  while (c){
    c = next_line_of_comment(c, MATRIX_MARKET_MAX_LINE, line);
    if (c){
      if (strstr(line, pin_string)){
        c = *pins = next_line_of_data(c, MATRIX_MARKET_MAX_LINE, NULL);
      }
      else if (strstr(line, vwgt_string)){
        c = *vwgts = next_line_of_data(c, MATRIX_MARKET_MAX_LINE, NULL);
      }
      else if (strstr(line, ewgt_string)){
        c = *ewgts = next_line_of_data(c, MATRIX_MARKET_MAX_LINE, NULL);
      }
    }
  }

  return;
}

static char *next_line(char *buf, size_t max)
{
char *c = buf;
size_t sanity = 0;

  /* Skip current line, and any subsequent blank lines.
   * Return pointer to the next line. 
   */

  if (!buf) return NULL;

  /* go to end of line */
  while (*c && (*c != '\n') && (sanity++ < max)) c++;

  if (*c==0 || sanity == max || *(c+1) == 0)
    return NULL;

  c = first_char(c+1, max);

  return c;
}

static char *next_line_of_data(char *buf, size_t max, char *str)
{
char *c; 

  /* Skip current line, and any following comments or blank lines.
   * Return pointer to the next line. Also copy
   * the line to a null-terminated C string because sscanf
   * may perform really poorly on very long strings.
   */

  if (!buf) return NULL;

  if (str) *str = 0;

  c = buf;

  while (c){
    c = next_line(c, max);
    if (c && (c[0] != COMMENT_CHAR)){
      if (str) make_string(c, str);
      return c;
    }
  }

  return NULL;
}

static char *next_line_of_comment(char *buf, size_t max, char *str)
{
char *c; 

  /* Skip current line, and any following data or blank lines.
   * Return pointer to the next line. Also copy
   * the line to a null-terminated C string because sscanf
   * may perform really poorly on very long strings.
   */

  if (!buf) return NULL;

  if (str) *str = 0;

  c = buf;

  while (c){
    c = next_line(c, max);
    if (c && (c[0] == COMMENT_CHAR)){
      if (str) make_string(c, str);
      return c;
    }
  }

  return NULL;
}

/*************************************************************************/

static int my_vtx(int proc, 
       ZOLTAN_ID_TYPE vtx, ZOLTAN_ID_TYPE mymin, ZOLTAN_ID_TYPE mymax,
       int myrank, int nprocs, PARIO_INFO_PTR pio_info)
{
  int mine = -1;

  if (myrank >= nprocs) return 0;
  if (nprocs == 1) return 1;

  if (pio_info->init_dist_type == INITIAL_FILE){
    /* The process ID of the vertex is in the file */
    mine = ((proc % nprocs) == myrank);
  }
  else if (pio_info->init_dist_type == INITIAL_CYCLIC){
    /* Deal out the vertices in a random fashion */
    mine = ((int)(vtx % nprocs) == myrank);
  }
  else if (pio_info->init_dist_type == INITIAL_LINEAR){
    /* First process gets first nvtxs/nprocs vertices, and so on */
    mine = ((vtx >= mymin) && (vtx <= mymax));
  }

  return mine;
}
static int my_pin(ZOLTAN_ID_TYPE eid, ZOLTAN_ID_TYPE vid, int proc,
       ZOLTAN_ID_TYPE pin, ZOLTAN_ID_TYPE npins, 
       ZOLTAN_ID_TYPE mymin, ZOLTAN_ID_TYPE mymax,
       int myrank, int nprocs, PARIO_INFO_PTR pio_info)
{
  int mine = -1;

  if (nprocs == 1) return 1;

  if (pio_info->init_dist_pins == INITIAL_ZERO){
    /* Node zero initially has all pins */
    mine = (myrank == 0);
  }
  else if (pio_info->init_dist_pins == INITIAL_FILE){
    /* The process ID of the pin owner is in the file */
    mine = ((proc % nprocs) == myrank);
  }
  else if (pio_info->init_dist_pins == INITIAL_CYCLIC){
    /* Deal out the pins in a random fashion */
    mine = ((int)(pin % nprocs) == myrank);
  }
  else if (pio_info->init_dist_pins == INITIAL_LINEAR){
    /* First process gets first npins/nprocs pins, and so on */
    mine = ((pin >= mymin) && (pin <= mymax));
  }
  else if (pio_info->init_dist_pins == INITIAL_ROW){
    /* Each process gets entire rows (hyperedges) of pins, no
       row is split across processes  */

    mine = ((int)(eid % nprocs) == myrank);
  }
  else if (pio_info->init_dist_pins == INITIAL_COL){
    /* Each process gets entire columns of pins, no column is split
       across processes  */
    mine = ((int)(vid % nprocs) == myrank);
  }

  return mine;
}
/*************************************************************************/
static int _zoltan_sortFunc(const void *a, const void *b)
{
  ZOLTAN_ID_TYPE ia, ib;

  ia = *(ZOLTAN_ID_TYPE *)a;
  ib = *(ZOLTAN_ID_TYPE *)b;

  if (ia < ib){
    return -1;
  }
  else if (ia > ib){
    return 1;
  }
  else{
    return 0;
  }
}

/*************************************************************************/
#if PRINT_DEBUG_INFO
static void debug_elements(int Proc, int Num_Proc, int num, ELEM_INFO_PTR el)
{
  int i,e;
  for (i=0; i<Num_Proc; i++){
    if (i == Proc){
      printf("Process %d (%d elements):\n",i,num);
      for (e=0; e<num; e++){
	if (e%20==0) printf("\n    ");
	printf("%" ZOLTAN_ID_SPEC,el[e].globalID);
      }
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(zoltan_get_global_comm());
  }
}
static void debug_lists(int Proc, int Num_Proc, int nedge, int *index, ZOLTAN_ID_TYPE *vtx, int *vtx_proc, ZOLTAN_ID_TYPE *egid)
{
  int i,e,v,nvtxs;
  for (i=0; i<Num_Proc; i++){
    if (i == Proc){
      printf("Process %d\n",i);
      for (e=0; e<nedge; e++){
	nvtxs = index[e+1]-index[e];
	printf("%" ZOLTAN_ID_SPEC " ) ",egid[e]);
	for (v=0; v<nvtxs; v++){
	  if (v && (v%10==0)) printf("\n    ");
	  printf("%" ZOLTAN_ID_SPEC " (%d) ",*vtx++,*vtx_proc++);
	}
	printf("\n");
      }
      fflush(stdout);
    }
    MPI_Barrier(zoltan_get_global_comm());
  }
}
static void debug_pins(int Proc, int Num_Proc,  
          ZOLTAN_ID_TYPE nGlobalEdges, ZOLTAN_ID_TYPE nGlobalVtxs, 
          int vtxWDim, int edgeWDim,
          int nMyPins, ZOLTAN_ID_TYPE *myPinI, ZOLTAN_ID_TYPE *myPinJ,
          int nMyVtx, ZOLTAN_ID_TYPE *myVtxNum, float *myVtxWgts,
          int nMyEdgeWgts, ZOLTAN_ID_TYPE *myEWGno, float *myEdgeWgts)
{
int p,i,j,k;

  MPI_Barrier(zoltan_get_global_comm());
  for (p=0; p < Num_Proc; p++){
    if (p == Proc){
      printf("Process: %d\n",p);
      printf(ZOLTAN_ID_SPEC " global edges, " ZOLTAN_ID_SPEC " global vertices, vwgtdim %d, ewgtdim %d\n",
              nGlobalEdges, nGlobalVtxs, vtxWDim, edgeWDim);
      printf("%d pins: ",nMyPins);
      for (i=0; i < nMyPins; i++) printf("%d,%d ",myPinI[i],myPinJ[i]);
      printf("\n");
      printf("%d vertices: ",nMyVtx);
      for (i=0, j=0; i < nMyVtx; i++){
        printf(ZOLTAN_ID_SPEC " ",myVtxNum[i]);
        for (k=0; k < vtxWDim; k++){
          printf("%f", myVtxWgts[j++]);
          if (k < vtxWDim-1) printf(",");
        }
        if (vtxWDim > 0) printf(" ");
      }
      printf("\n");
      printf("%d edge weights: ",nMyEdgeWgts);
      for (i=0, j=0; i < nMyEdgeWgts; i++){
        printf(ZOLTAN_ID_SPEC " ",myEWGno[i]);
        for (k=0; k < edgeWDim; k++){
          printf("%f", myEdgeWgts[j++]);
          if (k < edgeWDim-1) printf(",");
        }
        if (edgeWDim > 0) printf(" ");
      }
      printf("\n");
      fflush(stdout);
    } 
    MPI_Barrier(zoltan_get_global_comm());
  }
  MPI_Barrier(zoltan_get_global_comm());
}
#endif
