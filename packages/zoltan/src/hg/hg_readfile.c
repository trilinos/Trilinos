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

/*  This file is included only in hg_test, not in Zoltan. */

#include "hypergraph.h"

#define BUF_LEN 1000000


int hg_readfile (ZZ *zz, HGraph *hg, char *hgraphfile)
   {
   int ierr;
   FILE *f;
   char line1[81], errstr[200] ;
   char *yo = "hg_readfile" ;

   Zoltan_HG_HGraph_Init(hg);

   f = fopen (hgraphfile, "r") ;
   if (!f)
      {
      sprintf(errstr, "ERROR...not able to open file %s!\n",hgraphfile);
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, errstr) ;
      return ZOLTAN_FATAL;
      }

   ierr = Zoltan_HG_Readfile (0, f, &hg->nVtx, &hg->nEdge, &hg->nPin, 
    &hg->hindex, &hg->hvertex, &hg->VertexWeightDim, &hg->vwgt, 
    &hg->EdgeWeightDim, &hg->ewgt) ;
   if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      return ierr;

   ierr = Zoltan_HG_Create_Mirror (zz, hg);
   if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      return ierr;
   if (fclose(f))
      return ZOLTAN_WARN;
   return ZOLTAN_OK;
   }


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
