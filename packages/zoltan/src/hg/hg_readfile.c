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


#include "hypergraph.h"

#define BUF_LEN 1000000

int Zoltan_HG_Readfile   (int, FILE*, int*, int*, int*, int**, int**, int*, float**, int*, float**) ;


int HG_Readfile (ZZ *zz, HGraph *hg, char *hgraphfile)
   {
   FILE *f;
   char line1[81], errstr[200] ;
   char *yo = "Zoltan_HG_Readfile" ;

   Zoltan_HG_HGraph_Init(hg);

   f = fopen (hgraphfile, "r") ;
   if (!f)
      {
      sprintf(errstr, "ERROR...not able to open file %s!\n",hgraphfile);
      ZOLTAN_PRINT_WARN (zz->Proc, yo, errstr) ;
      return ZOLTAN_WARN;
      }

   Zoltan_HG_Readfile (0, f, &hg->nVtx, &hg->nEdge, &hg->nPin, &hg->hindex,
    &hg->hvertex, &hg->VertexWeightDim, &hg->vwgt, &hg->EdgeWeightDim, &hg->ewgt) ;

   if (Zoltan_HG_Create_Mirror (zz, hg))
      return ZOLTAN_WARN;
   if (fclose(f))
      return ZOLTAN_WARN;
   return ZOLTAN_OK;
   }


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
