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
  float *xWeight;
  int *Input_Parts = NULL, *xGNO = NULL;
  ZOLTAN_ID_PTR xGID=NULL, yGID=NULL;
  ZOLTAN_ID_PTR xLID=NULL;
  ZOLTAN_ID_PTR pinID=NULL;

/*   int final_output = 0; */

  ZOLTAN_TRACE_ENTER(zz, yo);

  memset (matrix, 0, sizeof(Zoltan_matrix)); /* Set all fields to 0 */

  /**************************************************/
  /* Obtain vertex information from the application */
  /**************************************************/

  ierr = Zoltan_Get_Obj_List(zz, &nX, &xGID, &xLID,
			     zz->Obj_Weight_Dim, &xWeight,
			     &Input_Parts);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
    goto End;
  }

  ZOLTAN_FREE(&Input_Parts); /* TODO: change --- Not used at this time */

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

  ierr = Zoltan_DD_Create (&matrix->ddX, zz->Communicator, zz->Num_GID, zz->Num_LID, 1, nX, 0);

  /* Make our new numbering public */
  Zoltan_DD_Update (matrix->ddX, xGID, xLID, (ZOLTAN_ID_PTR) xGNO, NULL, nX);

  ZOLTAN_FREE(&xGID);
  ZOLTAN_FREE(&xLID);
  ZOLTAN_FREE(&xGNO);

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

  matrix->yGNO = (int*)ZOLTAN_CALLOC(matrix->nY, sizeof(int));
  matrix->pinGNO = (int*)ZOLTAN_CALLOC(matrix->nPins, sizeof(int));
  if (((matrix->nY > 0 ) && (matrix->yGNO == NULL))
      || ((matrix->nPins > 0) && (matrix->pinGNO == NULL))){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* Convert pinID to pinGNO using the same translation as x */
  Zoltan_DD_Find (matrix->ddX, pinID, NULL, (ZOLTAN_ID_PTR)(matrix->pinGNO), NULL,
		  matrix->nPins, NULL);
  ZOLTAN_FREE(&pinID);

  if (enforceSquare) {
    /* Convert yGID to yGNO using the same translation as x */
    /* Needed for graph : rowID = colID */
    matrix->globalY = matrix->globalX;
    matrix->ddY = matrix->ddX;
  }
  else { /* Hyperedges name translation is different from the one of vertices */
/*     int nGlobalEdges = 0; */

    ierr = Zoltan_PHG_GIDs_to_global_numbers(zz, matrix->yGNO, matrix->nY,
					     randomizeInitDist, &matrix->globalY);

    ierr = Zoltan_DD_Create (&matrix->ddY, zz->Communicator, zz->Num_GID, 0, 1, matrix->nY, 0);

    /* Make our new numbering public */
    Zoltan_DD_Update (matrix->ddY, yGID, NULL, (ZOLTAN_ID_PTR) matrix->yGNO, NULL, matrix->nY);

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
  }

  Zoltan_DD_Find (matrix->ddY, yGID, NULL, (ZOLTAN_ID_PTR)(matrix->yGNO), NULL,
		  matrix->nY, NULL);

 End:
  ZOLTAN_FREE(&xGID);
  ZOLTAN_FREE(&xLID);
  ZOLTAN_FREE(&xGNO);
  ZOLTAN_FREE(&yGID);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}




#ifdef __cplusplus
}
#endif
