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
#include "matrix.h"

static int
matrix_get_edges(ZZ *zz, Zoltan_matrix *matrix, ZOLTAN_ID_PTR *yGID, ZOLTAN_ID_PTR *pinID,
		 int nX, ZOLTAN_ID_PTR *xGID, ZOLTAN_ID_PTR *xLID, int** xGNO, float** xwgt);

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
Zoltan_Matrix_Build (ZZ* zz, Zoltan_matrix_options *opt, Zoltan_matrix* matrix)
{
  static char *yo = "Zoltan_Matrix_Build";
  int ierr = ZOLTAN_OK;
  int nX;
  int  *xGNO = NULL;
  ZOLTAN_ID_PTR xLID=NULL;
  ZOLTAN_ID_PTR xGID=NULL;
  ZOLTAN_ID_PTR yGID=NULL;
  ZOLTAN_ID_PTR pinID=NULL;
  float *xwgt = NULL;
  int * Input_Parts=NULL;
  struct Zoltan_DD_Struct *dd = NULL;
  int *proclist = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

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


  if (nX) {
    xGNO = (int*) ZOLTAN_MALLOC(nX*sizeof(int));
    if (xGNO == NULL)
      MEMORY_ERROR;
  }
  /*******************************************************************/
  /* Assign vertex consecutive numbers (gnos)                        */
  /*******************************************************************/

  ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, xGNO, nX, matrix->opts.randomize, &matrix->globalX);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error assigning global numbers to vertices");
    goto End;
  }

  ierr = Zoltan_DD_Create (&dd, zz->Communicator, zz->Num_GID, 1,
			   sizeof(float)/sizeof(int)*zz->Obj_Weight_Dim, nX, 0);
  CHECK_IERR;

  /* Make our new numbering public */
  Zoltan_DD_Update (dd, xGID, (ZOLTAN_ID_PTR) xGNO, (ZOLTAN_ID_PTR)xwgt,  NULL, nX);

  /* I store : xGNO, xGID, xwgt, Input_Part */
  ierr = Zoltan_DD_Create (&matrix->ddX, zz->Communicator, 1, zz->Num_GID,
			   zz->Obj_Weight_Dim*sizeof(float)/sizeof(int), matrix->globalX/zz->Num_Proc, 0);
  CHECK_IERR;

  /* Hope a linear assignment will help a little */
  Zoltan_DD_Set_Neighbor_Hash_Fn1(matrix->ddX, matrix->globalX/zz->Num_Proc);
  /* Associate all the data with our xGNO */
  Zoltan_DD_Update (matrix->ddX, (ZOLTAN_ID_PTR)xGNO, xGID, (ZOLTAN_ID_PTR) xwgt, Input_Parts, nX);

  ZOLTAN_FREE(&Input_Parts);


  if (matrix->opts.pinwgt)
    matrix->pinwgtdim = zz->Edge_Weight_Dim;
  else
    matrix->pinwgtdim = 0;

  ierr = matrix_get_edges(zz, matrix, &yGID, &pinID, nX, &xGID, &xLID, &xGNO, &xwgt);
  CHECK_IERR;
  matrix->nY_ori = matrix->nY;

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  if (matrix->opts.enforceSquare && matrix->redist) {
    /* Convert yGID to yGNO using the same translation as x */
    /* Needed for graph : rowID = colID */
    /* y and x may have different distributions */
    matrix->ywgt = (float*)ZOLTAN_MALLOC(matrix->ywgtdim * sizeof(float));
    if (matrix->ywgtdim && matrix->ywgt == NULL)
      MEMORY_ERROR;
    Zoltan_DD_Find (dd, yGID, (ZOLTAN_ID_PTR)(matrix->yGNO), (ZOLTAN_ID_PTR)matrix->ywgt, NULL,
		    matrix->nY, NULL);
  }

  matrix->pinGNO = (int*)ZOLTAN_MALLOC(matrix->nPins* sizeof(int));
  if ((matrix->nPins > 0) && (matrix->pinGNO == NULL)) MEMORY_ERROR;

  if (matrix->opts.local) { /* keep only local edges */
    proclist = (int*) ZOLTAN_MALLOC(matrix->nPins*sizeof(int));
    if (matrix->nPins && proclist == NULL) MEMORY_ERROR;
  }
  else
    proclist = NULL;
  /* Convert pinID to pinGNO using the same translation as x */
  Zoltan_DD_Find (dd, pinID, (ZOLTAN_ID_PTR)(matrix->pinGNO), NULL, NULL,
		  matrix->nPins, proclist);
  ZOLTAN_FREE(&pinID);
  Zoltan_DD_Destroy(&dd);
  dd = NULL;

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
    matrix->yGNO = (int*)ZOLTAN_CALLOC(matrix->nY, sizeof(int));
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
      ierr = Zoltan_DD_Create (&matrix->ddY, zz->Communicator, 1, zz->Num_GID,
			       matrix->ywgtdim*sizeof(float)/sizeof(int), matrix->globalY/zz->Num_Proc, 0);
      /* Hope a linear assignment will help a little */
      Zoltan_DD_Set_Neighbor_Hash_Fn1(matrix->ddY, matrix->globalY/zz->Num_Proc);
      /* Associate all the data with our yGNO */
      Zoltan_DD_Update (matrix->ddY, (ZOLTAN_ID_PTR)matrix->yGNO,
			yGID, (ZOLTAN_ID_PTR) matrix->ywgt, NULL, matrix->nY);
  }

 End:
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
  ZOLTAN_FREE(&matrix->ywgt);

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
		 ZOLTAN_ID_PTR *xGID, ZOLTAN_ID_PTR *xLID, int **xGNO, float **xwgt)
{
  static char *yo = "Zoltan_Matrix_Build";
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

  if (graph_callbacks && hypergraph_callbacks){
/*     if (hgraph_model == GRAPH) */
/*       hypergraph_callbacks = 0; */
    graph_callbacks = 1; /* I prefer graph (allow to do "inplace") ! */
  }

  if (hypergraph_callbacks) {
    matrix->redist = 1;
    ZOLTAN_FREE(xGID);
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
    *xGNO = NULL;
    *yGID = NULL;
    matrix->ywgtdim = zz->Obj_Weight_Dim;
    matrix->ywgt = *xwgt;
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
    ZOLTAN_FREE(xGID);
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
