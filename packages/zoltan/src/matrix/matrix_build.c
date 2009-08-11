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

int
Zoltan_Matrix_Build (ZZ* zz, Zoltan_matrix* matrix)
{
  static char *yo = "Zoltan_Matrix_Build";
  int ierr = ZOLTAN_OK;
  int randomizeInitDist = 0;
  int enforceSquare = 0;
  int nX;
  int  *xGNO = NULL;
  ZOLTAN_ID_PTR xLID=NULL;
  ZOLTAN_ID_PTR xGID=NULL;
  ZOLTAN_ID_PTR yGID=NULL;
  ZOLTAN_ID_PTR pinID=NULL;
  float *xwgt = NULL, *ywgt = NULL;
  int * Input_Parts=NULL;
  struct Zoltan_DD_Struct *dd = NULL;

/*   int final_output = 0; */

  ZOLTAN_TRACE_ENTER(zz, yo);

  memset (matrix, 0, sizeof(Zoltan_matrix)); /* Set all fields to 0 */

  /**************************************************/
  /* Obtain vertex information from the application */
  /**************************************************/

  ierr = Zoltan_Get_Obj_List(zz, &nX, &xGID, &xLID,
			     zz->Obj_Weight_Dim, &xwgt,
			     &Input_Parts);
  ZOLTAN_FREE(&xLID);

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

  ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, xGNO, nX,
					   randomizeInitDist, &matrix->globalX);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error assigning global numbers to vertices");
    goto End;
  }

  ierr = Zoltan_DD_Create (&dd, zz->Communicator, zz->Num_GID, 0, 1, nX, 0);

  /* Make our new numbering public */
  Zoltan_DD_Update (dd, xGID, NULL, (ZOLTAN_ID_PTR) xGNO, NULL, nX);

  /* I store : xGNO, xGID, xwgt, Input_Part */
  ierr = Zoltan_DD_Create (&matrix->ddX, zz->Communicator, 1, zz->Num_GID,
			   zz->Obj_Weight_Dim*sizeof(float)/sizeof(int), matrix->globalX/zz->Num_Proc, 0);
  /* Hope a linear assignment will help a little */
  Zoltan_DD_Set_Neighbor_Hash_Fn1(matrix->ddX, matrix->globalX/zz->Num_Proc);
  /* Associate all the data with our xGNO */
  Zoltan_DD_Update (matrix->ddX, (ZOLTAN_ID_PTR)xGNO, xGID, (ZOLTAN_ID_PTR) xwgt, Input_Parts, nX);

  ZOLTAN_FREE(&xGNO);
  ZOLTAN_FREE(&xGID);
  ZOLTAN_FREE(&xwgt);
  ZOLTAN_FREE(&Input_Parts);

  /*
   * Each processor:
   *   owns a set of pins (nonzeros)
   *   may provide some edge weights
   *
   * We assume that no two processes will supply the same pin.
   * But more than one process may supply pins for the same edge.
   */

  ierr = Zoltan_Hypergraph_Queries(zz, &matrix->nY,
				   &matrix->nPins, &yGID, &matrix->ystart,
				   &pinID);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }
  matrix->yend = matrix->ystart + 1;

  matrix->pinGNO = (int*)ZOLTAN_CALLOC(matrix->nPins, sizeof(int));
  if ((matrix->nPins > 0) && (matrix->pinGNO == NULL)) MEMORY_ERROR;
  /* Convert pinID to pinGNO using the same translation as x */
  Zoltan_DD_Find (dd, pinID, NULL, (ZOLTAN_ID_PTR)(matrix->pinGNO), NULL,
		  matrix->nPins, NULL);
  ZOLTAN_FREE(&pinID);
  Zoltan_DD_Destroy(&dd);
  dd = NULL;

  if (enforceSquare) {
    /* Convert yGID to yGNO using the same translation as x */
    /* Needed for graph : rowID = colID */
    matrix->globalY = matrix->globalX;
    matrix->ddY = matrix->ddX;
  }
  else { /* Hyperedges name translation is different from the one of vertices */
    matrix->yGNO = (int*)ZOLTAN_CALLOC(matrix->nY, sizeof(int));
    if (matrix->nY && matrix->yGNO == NULL) MEMORY_ERROR;

    /*     int nGlobalEdges = 0; */
    ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, matrix->yGNO, matrix->nY,
					     randomizeInitDist, &matrix->globalY);

    ierr = Zoltan_DD_Create (&dd, zz->Communicator, zz->Num_GID, 0, 1, matrix->nY, 0);

    /* Make our new numbering public */
    Zoltan_DD_Update (dd, yGID, NULL, (ZOLTAN_ID_PTR) matrix->yGNO, NULL, matrix->nY);

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
    Zoltan_DD_Find (dd, yGID, NULL, (ZOLTAN_ID_PTR)(matrix->yGNO), NULL,
		    matrix->nY, NULL);
    Zoltan_DD_Destroy(&dd);
    dd = NULL;

    /* We have to define ddY : yGNO, yGID, ywgt */
    ierr = Zoltan_DD_Create (&matrix->ddY, zz->Communicator, 1, zz->Num_GID,
			     matrix->ywgtdim*sizeof(float)/sizeof(int), matrix->globalY/zz->Num_Proc, 0);
    /* Hope a linear assignment will help a little */
    Zoltan_DD_Set_Neighbor_Hash_Fn1(matrix->ddY, matrix->globalY/zz->Num_Proc);
    /* Associate all the data with our xGNO */
    Zoltan_DD_Update (matrix->ddY, (ZOLTAN_ID_PTR)matrix->yGNO,
		      yGID, (ZOLTAN_ID_PTR) ywgt, NULL, matrix->nY);
  }


 End:
  ZOLTAN_FREE(&xLID);
  ZOLTAN_FREE(&xGNO);
  ZOLTAN_FREE(&xGID);
  ZOLTAN_FREE(&xwgt);
  ZOLTAN_FREE(&Input_Parts);
  if (dd != NULL)
    Zoltan_DD_Destroy(&dd);
  ZOLTAN_FREE(&yGID);
  ZOLTAN_FREE(&ywgt);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}




#ifdef __cplusplus
}
#endif
