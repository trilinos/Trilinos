// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <math.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "phg.h"
#include "zoltan_matrix.h"

static int
matrix_get_edges(ZZ *zz, Zoltan_matrix *matrix, ZOLTAN_ID_PTR *yGID, ZOLTAN_ID_PTR *pinID,
		 int nX, ZOLTAN_ID_PTR *xGID, ZOLTAN_ID_PTR *xLID, ZOLTAN_GNO_TYPE ** xGNO, float** xwgt, int use_full_dd);

  /* In build_graph.c, may be moved here soon */
extern int Zoltan_Get_Num_Edges_Per_Obj(
  ZZ *zz,
  int num_obj,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int **edges_per_obj,
  int *max_edges,
  int *num_edges
  );

int
Zoltan_Matrix_Build (ZZ* zz, Zoltan_matrix_options *opt, Zoltan_matrix* matrix,
  int request_GNOs,                /* Input:  Flag indicating calling code 
                                              needs translation of extra GIDs
                                              to GNOs; partial 2D coloring
                                              needs this feature. */
  int num_requested,               /* Input:  Local # of GIDs needing 
                                              translation to GNOs. */
  ZOLTAN_ID_PTR requested_GIDs,    /* Input:  Calling code requests the 
                                              GNOs for these GIDs */
  ZOLTAN_GNO_TYPE *requested_GNOs  /* Output: Return GNOs of 
                                              the requested GIDs.  */
)  
{
  static char *yo = "Zoltan_Matrix_Build";
  int ierr = ZOLTAN_OK;
  int nX;
  ZOLTAN_GNO_TYPE tmp;
  ZOLTAN_GNO_TYPE *xGNO = NULL;
  ZOLTAN_ID_PTR xLID=NULL;
  ZOLTAN_ID_PTR xGID=NULL;
  ZOLTAN_ID_PTR yGID=NULL;
  ZOLTAN_ID_PTR pinID=NULL;
  float *xwgt = NULL;
  int * Input_Parts=NULL;
  struct Zoltan_DD_Struct *dd = NULL;
  int *proclist = NULL;
  int *xpid = NULL;
  int i;
  int gno_size_for_dd;
  MPI_Datatype zoltan_gno_mpi_type;
  int use_full_dd = (opt->speed == MATRIX_FULL_DD);
  int fast_build_base = opt->fast_build_base;
  matrix->opts.speed = opt->speed;  
  matrix->opts.fast_build_base = opt->fast_build_base;

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (num_requested && (!requested_GIDs || !requested_GNOs)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Error in requested input; needed arrays are NULL.\n");
  }

  /* ZOLTAN_GNO_TYPE is >= ZOLTAN_ID_TYPE */
  gno_size_for_dd = sizeof(ZOLTAN_GNO_TYPE) / sizeof(ZOLTAN_ID_TYPE);

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

  memset (matrix, 0, sizeof(Zoltan_matrix)); /* Set all fields to 0 */
  memcpy (&matrix->opts, opt, sizeof(Zoltan_matrix_options));

  /**************************************************/
  /* Obtain vertex information from the application */
  /**************************************************/

  ierr = Zoltan_Get_Obj_List(zz, &nX, &xGID, &xLID,
			     zz->Obj_Weight_Dim, &xwgt,
			     &Input_Parts);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
    goto End;
  }
  ZOLTAN_FREE(&Input_Parts);
  ZOLTAN_FREE(&xwgt);

  /*******************************************************************/
  /* Assign vertex consecutive numbers (gnos)                        */
  /*******************************************************************/

  if (use_full_dd) {
    /* Zoltan computes a translation */
    /* Have to use Data Directory if request_GNOs is true. */
    if (nX) {
      xGNO = (ZOLTAN_GNO_TYPE*) ZOLTAN_MALLOC(nX*sizeof(ZOLTAN_GNO_TYPE));
      if (xGNO == NULL)
	MEMORY_ERROR;
    }
    ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, xGNO, nX, matrix->opts.randomize, &matrix->globalX);

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error assigning global numbers to vertices");
      goto End;
    }

    ierr = Zoltan_DD_Create (&dd, zz->Communicator, zz->Num_GID, gno_size_for_dd, 0, nX, 0);
    CHECK_IERR;

    /* Make our new numbering public */
    Zoltan_DD_Update (dd, xGID, (ZOLTAN_ID_PTR) xGNO, NULL,  NULL, nX);
    if (request_GNOs) {
      Zoltan_DD_Find(dd, requested_GIDs, (ZOLTAN_ID_PTR) requested_GNOs,
                     NULL, NULL, num_requested, NULL);
    }
  }
  else { /* We don't want to use the DD */
     /*
     * KDDKDD 2/10/11  This code cannot work when NUM_GID_ENTRIES>1.
     * KDDKDD 2/10/11  The assumption is that, if a user sets the
     * KDDKDD 2/10/11  appropriate parameter to enable this code, the user
     * KDDKDD 2/10/11  knows that his GIDs are compatible with integers.
     */
    if (sizeof(ZOLTAN_GNO_TYPE) != sizeof(ZOLTAN_ID_TYPE)){
      xGNO = (ZOLTAN_GNO_TYPE*) ZOLTAN_MALLOC(nX*sizeof(ZOLTAN_GNO_TYPE));
      if (nX && xGNO == NULL)
        MEMORY_ERROR;
      for (i=0; i < nX; i++)
        xGNO[i] = (ZOLTAN_GNO_TYPE)xGID[i] - fast_build_base;
    }
    else {
      xGNO = (ZOLTAN_GNO_TYPE *)xGID;
      if (fast_build_base)
        for (i = 0; i < nX; i++)
          xGNO[i] -= fast_build_base;
    }

    for (i = 0; i < num_requested; i++)
      requested_GNOs[i] = (ZOLTAN_GNO_TYPE)requested_GIDs[i]
                        - fast_build_base;
     
    tmp = (ZOLTAN_GNO_TYPE)nX; 
    MPI_Allreduce(&tmp, &matrix->globalX, 1, zoltan_gno_mpi_type, MPI_SUM, zz->Communicator);
  }

  /* I store : xGNO, xGID, xpid,  */

  ierr = Zoltan_DD_Create (&matrix->ddX, zz->Communicator, gno_size_for_dd, zz->Num_GID,
			   sizeof(int), matrix->globalX/zz->Num_Proc, 0);
  CHECK_IERR;

  /* Hope a linear assignment will help a little */
  if (matrix->globalX/zz->Num_Proc)
    Zoltan_DD_Set_Neighbor_Hash_Fn1(matrix->ddX, matrix->globalX/zz->Num_Proc);
  /* Associate all the data with our xGNO */
  xpid = (int*)ZOLTAN_MALLOC(nX*sizeof(int));
  if (nX >0 && xpid == NULL) MEMORY_ERROR;
  for (i = 0 ; i < nX ; ++i)
    xpid[i] = zz->Proc;

  Zoltan_DD_Update (matrix->ddX, (ZOLTAN_ID_PTR)xGNO, xGID, (char *)xpid, NULL, nX);
  ZOLTAN_FREE(&xpid);

  if (matrix->opts.pinwgt)
    matrix->pinwgtdim = zz->Edge_Weight_Dim;
  else
    matrix->pinwgtdim = 0;

  ierr = matrix_get_edges(zz, matrix, &yGID, &pinID, nX, &xGID, &xLID, &xGNO, &xwgt, use_full_dd);
  CHECK_IERR;
  matrix->nY_ori = matrix->nY;

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  if (matrix->opts.enforceSquare && matrix->redist) {
    /* Convert yGID to yGNO using the same translation as x */
    /* Needed for graph : rowID = colID */
    /* y and x may have different distributions */
    matrix->yGNO = (ZOLTAN_GNO_TYPE*)ZOLTAN_MALLOC(matrix->nY * sizeof(ZOLTAN_GNO_TYPE));
    if (matrix->nY && matrix->yGNO == NULL) {
      ZOLTAN_FREE(&pinID);
      MEMORY_ERROR;
    }
    ierr = Zoltan_DD_Find (dd, yGID, (ZOLTAN_ID_PTR)(matrix->yGNO), NULL, NULL,
		    matrix->nY, NULL);
    if (ierr != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc,yo,"Hyperedge GIDs don't match.\n");
      ierr = ZOLTAN_FATAL;
      ZOLTAN_FREE(&pinID);
      goto End;
    }
  }

  if (matrix->opts.local) { /* keep only local edges */
    proclist = (int*) ZOLTAN_MALLOC(matrix->nPins*sizeof(int));
    if (matrix->nPins && proclist == NULL) { 
      ZOLTAN_FREE(&pinID);
      MEMORY_ERROR;
    }
  }
  else
    proclist = NULL;

  /* Convert pinID to pinGNO using the same translation as x */
  if (use_full_dd) {
    matrix->pinGNO = (ZOLTAN_GNO_TYPE*)ZOLTAN_MALLOC(matrix->nPins* sizeof(ZOLTAN_GNO_TYPE));
    if ((matrix->nPins > 0) && (matrix->pinGNO == NULL)) {
        ZOLTAN_FREE(&pinID);
        MEMORY_ERROR;
    }

    ierr = Zoltan_DD_Find (dd, pinID, (ZOLTAN_ID_PTR)(matrix->pinGNO), NULL, NULL,
			   matrix->nPins, proclist);
    if (ierr != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc,yo,"Undefined GID found.\n");
      ierr = ZOLTAN_FATAL;
      goto End;
    }
    ZOLTAN_FREE(&pinID);
    Zoltan_DD_Destroy(&dd);
    dd = NULL;
  }
  else {
    if (sizeof(ZOLTAN_GNO_TYPE) != sizeof(ZOLTAN_ID_TYPE)){
      matrix->pinGNO = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(matrix->nPins * sizeof(ZOLTAN_GNO_TYPE));
      if (matrix->nPins && !matrix->pinGNO){
        ZOLTAN_FREE(&pinID);
        MEMORY_ERROR;
      }
      for (i=0; i < matrix->nPins; i++)
        matrix->pinGNO[i] = (ZOLTAN_GNO_TYPE)pinID[i] - fast_build_base;
      
      ZOLTAN_FREE(&pinID);
    }
    else{
      matrix->pinGNO = (ZOLTAN_GNO_TYPE *) pinID;
      if (fast_build_base)
        for (i=0; i < matrix->nPins; i++)
          matrix->pinGNO[i] -= fast_build_base;
      pinID = NULL;
    }
  }

/*   if (matrix->opts.local) {  /\* keep only local edges *\/ */
/*     int *nnz_list; /\* nnz offset to delete *\/ */
/*     int nnz;       /\* number of nnz to delete *\/ */
/*     int i; */

/*     nnz_list = (int*) ZOLTAN_MALLOC(matrix->nPins*sizeof(int)); */
/*     if (matrix->nPins && nnz_list == NULL) MEMORY_ERROR; */
/*     for (i = 0, nnz=0 ; i < matrix->nPins ; ++i) { */
/*       if (proclist[i] == zz->Proc) continue; */
/*       nnz_list[nnz++] = i; */
/*     } */
/*     ZOLTAN_FREE(&proclist); */
/*     Zoltan_Matrix_Delete_nnz(zz, matrix, nnz, nnz_list); */
/*   } */

  if (!matrix->opts.enforceSquare) {
    /* Hyperedges name translation is different from the one of vertices */
    matrix->yGNO = (ZOLTAN_GNO_TYPE*)ZOLTAN_CALLOC(matrix->nY, sizeof(ZOLTAN_GNO_TYPE));
    if (matrix->nY && matrix->yGNO == NULL) MEMORY_ERROR;

    /*     int nGlobalEdges = 0; */
    ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, matrix->yGNO, matrix->nY,
					     matrix->opts.randomize, &matrix->globalY);
    CHECK_IERR;

/*     /\**************************************************************************************** */
/*      * If it is desired to remove dense edges, divide the list of edges into */
/*      * two lists.  The ZHG structure will contain the removed edges (if final_output is true), */
/*      * and the kept edges will be returned. */
/*      ****************************************************************************************\/ */
/*     totalNumEdges = zhg->globalHedges; */

/*     ierr = remove_dense_edges_matrix(zz, zhg, edgeSizeThreshold, final_output, */
/*				     &nLocalEdges, &nGlobalEdges, &nPins, */
/*				     &edgeGNO, &edgeSize, &edgeWeight, &pinGNO, &pinProcs); */

/*     if (nGlobalEdges < totalNumEdges){ */
/*       /\* re-assign edge global numbers if any edges were removed *\/ */
/*       ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, edgeGNO, nLocalEdges, */
/*					       randomizeInitDist, &totalNumEdges); */
/*       if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) { */
/*	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error reassigning global numbers to edges"); */
/*	goto End; */
/*       } */
/*     } */

      /* We have to define ddY : yGNO, yGID, ywgt */
      ierr = Zoltan_DD_Create (&matrix->ddY, zz->Communicator, gno_size_for_dd, zz->Num_GID,
			       0, matrix->globalY/zz->Num_Proc, 0);
      /* Hope a linear assignment will help a little */
      if (matrix->globalY/zz->Num_Proc)
        Zoltan_DD_Set_Neighbor_Hash_Fn1(matrix->ddY, matrix->globalY/zz->Num_Proc);
      /* Associate all the data with our yGNO */
      Zoltan_DD_Update (matrix->ddY, (ZOLTAN_ID_PTR)matrix->yGNO, yGID, NULL, NULL, matrix->nY);
  }

 End:
  ZOLTAN_FREE(&xpid);
  ZOLTAN_FREE(&xLID);
  ZOLTAN_FREE(&xGNO);
  ZOLTAN_FREE(&xGID);
  ZOLTAN_FREE(&xwgt);
  ZOLTAN_FREE(&Input_Parts);
  ZOLTAN_FREE(&proclist);
  if (dd != NULL)
    Zoltan_DD_Destroy(&dd);
  /* Already stored in the DD */
  ZOLTAN_FREE(&yGID);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}


int
Zoltan_Matrix_Vertex_Info(ZZ* zz, const Zoltan_matrix * const m,
			  ZOLTAN_ID_PTR lid,
			  float *wwgt, int *input_part)
{
  static char *yo = "Zoltan_Matrix_Vertex_Info";
  int ierr = ZOLTAN_OK;
  int nX;
  ZOLTAN_ID_PTR l_gid = NULL;
  ZOLTAN_ID_PTR l_lid = NULL;
  float * l_xwgt = NULL;
  int *l_input_part = NULL;
  struct Zoltan_DD_Struct *dd = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (m->completed == 0) {
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  ierr = Zoltan_Get_Obj_List(zz, &nX, &l_gid, &l_lid,
			     zz->Obj_Weight_Dim, &l_xwgt,
			     &l_input_part);

  ierr = Zoltan_DD_Create (&dd, zz->Communicator, zz->Num_GID, zz->Num_LID, 
                           sizeof(float) * zz->Obj_Weight_Dim, nX, 0);
  CHECK_IERR;

    /* Make our new numbering public */
  Zoltan_DD_Update (dd, l_gid, l_lid, (char *) l_xwgt,l_input_part, nX);
  ZOLTAN_FREE(&l_gid);
  ZOLTAN_FREE(&l_lid);
  ZOLTAN_FREE(&l_xwgt);
  ZOLTAN_FREE(&l_input_part);

  ierr = Zoltan_DD_Find (dd, m->yGID, lid, (char *)wwgt, input_part,
		    m->nY, NULL);

 End:
  if (dd != NULL)
    Zoltan_DD_Destroy(&dd);
  ZOLTAN_FREE(&l_gid);
  ZOLTAN_FREE(&l_lid);
  ZOLTAN_FREE(&l_xwgt);
  ZOLTAN_FREE(&l_input_part);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


  /*
   * Each processor:
   *   owns a set of pins (nonzeros)
   *   may provide some edge weights
   *
   * We assume that no two processes will supply the same pin.
   * But more than one process may supply pins for the same edge.
   */
static int
matrix_get_edges(ZZ *zz, Zoltan_matrix *matrix, ZOLTAN_ID_PTR *yGID, ZOLTAN_ID_PTR *pinID, int nX,
		 ZOLTAN_ID_PTR *xGID, ZOLTAN_ID_PTR *xLID, ZOLTAN_GNO_TYPE **xGNO, float **xwgt, int use_full_dd)
{
  static char *yo = "matrix_get_edges";
  int ierr = ZOLTAN_OK;
  int hypergraph_callbacks = 0, graph_callbacks = 0;
  int *nbors_proc = NULL; /* Pointers are global for the function to ensure proper free */
  int *edgeSize = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);
  if (zz->Get_HG_Size_CS && zz->Get_HG_CS) {
    hypergraph_callbacks = 1;
  }
  if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
           (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) {
    graph_callbacks = 1;
  }

/* TEMPORARY FIX */
if (!graph_callbacks){
  fprintf(stderr,"Bug #5470: matrix_get_edges fails for hypergraph queries\n");
  return ZOLTAN_FATAL; 
}
hypergraph_callbacks=0;
/* TEMPORARY FIX */

  if (graph_callbacks && hypergraph_callbacks){
/*     if (hgraph_model == GRAPH) */
/*       hypergraph_callbacks = 0; */
    graph_callbacks = 1; /* I prefer graph (allow to do "inplace") ! */
  }

  if (hypergraph_callbacks) {
    matrix->redist = 1;
    if (use_full_dd || ((ZOLTAN_ID_PTR) *xGNO != *xGID))
      ZOLTAN_FREE(xGID);
    else
      *xGID = NULL;
    ZOLTAN_FREE(xLID);
    ZOLTAN_FREE(xGNO);
    ierr = Zoltan_Hypergraph_Queries(zz, &matrix->nY,
				     &matrix->nPins, yGID, &matrix->ystart,
				     pinID);
    CHECK_IERR;
    matrix->yend = matrix->ystart + 1;
  }
  else if (graph_callbacks) {
    int max_edges = 0;
    int vertex;
    int numGID, numLID;

    matrix->opts.enforceSquare = 1;
    matrix->nY = nX; /* It is square ! */
    matrix->yGNO = *xGNO;
    *yGID = NULL;
    matrix->ywgtdim = zz->Obj_Weight_Dim;
    *xwgt = NULL;

    numGID = zz->Num_GID;
    numLID = zz->Num_LID;

    /* TODO : support local graphs */
    /* TODO : support weights ! */
    /* Get edge data */
    Zoltan_Get_Num_Edges_Per_Obj(zz, matrix->nY, *xGID, *xLID,
				 &edgeSize, &max_edges, &matrix->nPins);

    (*pinID) = ZOLTAN_MALLOC_GID_ARRAY(zz, matrix->nPins);
    nbors_proc = (int *)ZOLTAN_MALLOC(matrix->nPins * sizeof(int));

    if (matrix->nPins && ((*pinID) == NULL || nbors_proc == NULL))
      MEMORY_ERROR;

    matrix->pinwgt = (float*)ZOLTAN_MALLOC(matrix->nPins*matrix->pinwgtdim*sizeof(float));
    if (matrix->nPins && matrix->pinwgtdim && matrix->pinwgt == NULL)
      MEMORY_ERROR;

    if (zz->Get_Edge_List_Multi) {
      zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data,
			      numGID, numLID,
			      matrix->nY, *xGID, *xLID,
			      edgeSize,
			      (*pinID), nbors_proc, matrix->pinwgtdim,
			      matrix->pinwgt, &ierr);
    }
    else {
      int edge;
      for (vertex = 0, edge = 0 ; vertex < matrix->nY ; ++vertex) {
	zz->Get_Edge_List(zz->Get_Edge_List_Data, numGID, numLID,
                          (*xGID)+vertex*numGID, (*xLID)+vertex*numLID,
                          (*pinID)+edge*numGID, nbors_proc+edge, matrix->pinwgtdim,
                          matrix->pinwgt+edge*matrix->pinwgtdim, &ierr);
	edge += edgeSize[vertex];
      }
    }
    CHECK_IERR;

    /* Not Useful anymore */
    ZOLTAN_FREE(xLID);
    if (use_full_dd || ((ZOLTAN_ID_PTR) *xGNO != *xGID)) {
      ZOLTAN_FREE(xGID);
      *xGNO = NULL;
    }
    else {
      *xGID = NULL;
      *xGNO = NULL;
    }
    ZOLTAN_FREE(&nbors_proc);

    /* Now construct CSR indexing */
    matrix->ystart = (int*) ZOLTAN_MALLOC((matrix->nY+1)*sizeof(int));
    if (matrix->ystart == NULL)
      MEMORY_ERROR;

    matrix->ystart[0] = 0;
    matrix->yend = matrix->ystart + 1;
    for (vertex = 0 ; vertex < matrix->nY ; ++vertex)
      matrix->ystart[vertex+1] = matrix->ystart[vertex] + edgeSize[vertex];
  }
  else {
    FATAL_ERROR ("You have to define Hypergraph or Graph queries");
  }

  if (matrix->opts.enforceSquare) {
    matrix->globalY = matrix->globalX;
    matrix->ddY = matrix->ddX;
    matrix->ywgtdim = zz->Obj_Weight_Dim;
  }

 End:
  ZOLTAN_FREE(&edgeSize);
  ZOLTAN_FREE(&nbors_proc);
  ZOLTAN_FREE(xLID);
  ZOLTAN_FREE(xGID);
  ZOLTAN_FREE(xGNO);
  ZOLTAN_FREE(xwgt);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

#ifdef __cplusplus
}
#endif
