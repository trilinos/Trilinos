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
 *****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "zz_sort.h"
#include "hsfc_hilbert_const.h"


int Zoltan_LocalHSFC_Order(
			   ZZ *zz,               /* Zoltan structure */
			   int num_obj,          /* Number of (local) objects to order. */
			   ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
			                         /* The application must allocate enough space */
			   ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
			                         /* The application must allocate enough space */
			   int *rank,            /* rank[i] is the rank of gids[i] */
			   int *iperm,
			   ZOOS *order_opt       /* Ordering options, parsed by Zoltan_Order */
                           )
{

  static char *yo = "Zoltan_LocalHSFC_Order";

  int n, ierr=ZOLTAN_OK;

  double (*fhsfc)(ZZ*, double*);  /* space filling curve function */

  int wgt_dim=0; 
  float *obj_wgts=0;
  int *parts=0;

  int numGeomDims=0;
  double *geomArray=0;

  /* Variables for bounding box */
  double *minValInDim;
  double *maxValInDim;
  double *widthDim;


  double *hsfcKey=0;
  int *coordIndx=0;

  /* Counters */
  int objNum;
  int dimNum;

  int offset=0;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  ZOLTAN_TRACE_ENTER(zz, yo);

  /******************************************************************/
  /* If for some reason order_opt is NULL, allocate a new ZOOS here. */
  /* This should really never happen. */
  /******************************************************************/
  if (!order_opt)
  {
    order_opt = (ZOOS *) ZOLTAN_MALLOC(sizeof(ZOOS));
    strcpy(order_opt->method,"LOCAL_HSFC");
  }
  /******************************************************************/

  /* local HSFC only computes the rank vector */
  order_opt->return_args = RETURN_RANK; 


  /******************************************************************/
  /* Check that num_obj equals the number of objects on this proc. */
  /* This constraint may be removed in the future. */
  /******************************************************************/
  n = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  if ((ierr!= ZOLTAN_OK) && (ierr!= ZOLTAN_WARN))
  {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Get_Num_Obj returned error.");
    return(ZOLTAN_FATAL);
  }
  if (n != num_obj)
  {
    /* Currently this is a fatal error. */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input num_obj does not equal the number of objects.");
    return(ZOLTAN_FATAL);
  }
  /******************************************************************/

  /******************************************************************/
  /* Get lists of objects                                           */
  /******************************************************************/
  ierr = Zoltan_Get_Obj_List(zz, &n, &gids, &lids, wgt_dim, &obj_wgts, &parts);
  if ((ierr!= ZOLTAN_OK) && (ierr!= ZOLTAN_WARN))
  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Get_Obj_List returned error.");
      return(ZOLTAN_FATAL);
  }
  /******************************************************************/

  /******************************************************************/
  /* Get geometry for objects*/
  /******************************************************************/
  ierr = Zoltan_Get_Coordinates(zz, n, gids, lids, &numGeomDims,
			       &geomArray);
  if (ierr != 0)
  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_Get_Coordinates returned error.");
      return(ZOLTAN_FATAL);
  }
  /******************************************************************/

  /******************************************************************/
  /* Place coords in bounding box                                   */
  /******************************************************************/
  minValInDim =  (double *) malloc(numGeomDims * sizeof (double));
  maxValInDim =  (double *) malloc(numGeomDims * sizeof (double));
  widthDim =  (double *) malloc(numGeomDims * sizeof (double));

  for(dimNum=0; dimNum<numGeomDims; dimNum++)
  {
    minValInDim[dimNum] = HUGE_VAL;
    maxValInDim[dimNum] = -HUGE_VAL;
  }

  /*************************************************************/
  /* Determine min, max, and width for each dimension          */
  /*************************************************************/
  for (objNum=0; objNum<n; objNum++)
  {
    for(dimNum=0; dimNum<numGeomDims; dimNum++)
    {
      if (geomArray[objNum * numGeomDims + dimNum] < minValInDim[dimNum])
      {
        minValInDim[dimNum] = geomArray[objNum * numGeomDims + dimNum];
      }
      if (geomArray[objNum * numGeomDims + dimNum] > maxValInDim[dimNum])
      {
        maxValInDim[dimNum] = geomArray[objNum * numGeomDims + dimNum];
      }
    }
  }

  for(dimNum=0; dimNum<numGeomDims; dimNum++)
  {
    widthDim[dimNum] = maxValInDim[dimNum] - minValInDim[dimNum]; 
  }
  /*************************************************************/

  /*************************************************************/
  /* Rescale values to fit in bounding box                     */
  /*************************************************************/
  for (objNum=0; objNum<n; objNum++)
  {
    for(dimNum=0; dimNum<numGeomDims; dimNum++)
    {
      geomArray[objNum * numGeomDims + dimNum] -= minValInDim[dimNum];
      geomArray[objNum * numGeomDims + dimNum] /= widthDim[dimNum];
    }
  }
  /*************************************************************/

  free(minValInDim); minValInDim=0;
  free(maxValInDim); maxValInDim=0;
  free(widthDim); widthDim=0;
  /******************************************************************/

  /******************************************************************/   
  /* Specify which HSFC function to use (based on dim) */
  /******************************************************************/
  if (numGeomDims==1)
  {
    fhsfc = Zoltan_HSFC_InvHilbert1d;
  }
  else if (numGeomDims==2)
  {
    fhsfc = Zoltan_HSFC_InvHilbert2d;
  }
  else if (numGeomDims==3)
  {
    fhsfc = Zoltan_HSFC_InvHilbert3d;
  }
  else /* this error should have been previously caught */
  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Geometry should be of dimension 1, 2, or 3.");
      return(ZOLTAN_FATAL);
  }
  /******************************************************************/

  /******************************************************************/
  /* Generate hsfc keys and indices to be sorted                    */
  /******************************************************************/
  hsfcKey = (double *) malloc(n * sizeof (double));
  coordIndx = (int *) malloc(n *sizeof(int));
  for (objNum=0; objNum<n; objNum++)
  {
    hsfcKey[objNum] = fhsfc(zz, &(geomArray[objNum * numGeomDims]) );
    coordIndx[objNum] = objNum;
  }
  /******************************************************************/

  /******************************************************************/
  /* Sort indices based on keys                                     */
  /******************************************************************/
  Zoltan_quicksort_pointer_dec_double (coordIndx, hsfcKey, 0, n-1);
  /******************************************************************/


  /******************************************************************/
  /* get ranks                                                      */
  /******************************************************************/

  /******************************************************/
  /* Determine offsets                                  */
  /******************************************************/
  MPI_Scan(&n, &offset, 1, MPI_INT, MPI_SUM, zz->Communicator);
  offset -= n; /* MPI_Scan is inclusive, so subtract off local size */
  /******************************************************/

  for(objNum=0; objNum<n; objNum++)
  {
    /*MMW temporary hack to make Cedric's interface give me want I need */
    /*rank[coordIndx[objNum]] = objNum + offset; */
    rank[objNum] = coordIndx[objNum] + offset; 
  }

  /******************************************************************/

  /* iperm is to be deprecated so not calculated*/

  free(hsfcKey);
  free(coordIndx);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ZOLTAN_OK);

}


#ifdef __cplusplus
}
#endif
