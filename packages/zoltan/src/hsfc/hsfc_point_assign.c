/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/



#include "hsfc.h"

/* For a detailed explaination of this module, please see the Developers
   Guide.  For instructions on its use, please see the Users Guide.

   Currently, this routine returns an error for points outside of the
   bounding box created during the load balancing phase.  This limitation
   is expected to be removed in the future.
*/



/* Point drop for refinement after above partitioning */
int Zoltan_HSFC_Point_Assign (ZZ *zz, double *x, int *proc)
   {
   double     scaled[3] ;
   double     fsfc ;
   Partition *p ;
   int        i ;
   HSFC_Data *d ;
   int        err ;
   char *yo = "Zoltan_HSFC_Point_Assign" ;

   ZOLTAN_TRACE_ENTER (zz, yo) ;
   d = (HSFC_Data *) zz->LB.Data_Structure ;
   if (d == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, 
          "No Decomposition Data available; use KEEP_CUTS parameter.");


   /* Insure that point is in bounding box */
   for (i = 0 ; i < d->ndimension ; i++)
      if ((x[i] > d->bbox_hi[i]) || (x[i] < d->bbox_lo[i]))
         ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "Point outside bounding box") ;

   /* Calculate scaled coordinates, calculate HSFC coordinate */
   for (i = 0 ; i < d->ndimension ; i++)
      scaled[i] = (x[i] - d->bbox_lo[i]) / d->bbox_extent[i] ;
   fsfc = d->fhsfc (scaled) ;           /* Note, this is a function call */

   /* Find partition containing point and return its number */
   p = (Partition *) bsearch (&fsfc, d->final_partition, zz->Num_Proc,
    sizeof (Partition), Zoltan_HSFC_compare) ;
   if (p == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "programming error, shouldn't happen") ;
   *proc = p->index ;
   err = ZOLTAN_OK ;

free:
   ZOLTAN_TRACE_EXIT (zz, yo) ;
   return err ;
   }

