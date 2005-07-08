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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "hsfc.h"

/* For a detailed explaination of this module, please see the Developers
   Guide.  For instructions on its use, please see the Users Guide.   */


/* Point drop for refinement after above partitioning */
int Zoltan_HSFC_Point_Assign (
   ZZ *zz, 
   double *x, 
   int *proc,
   int *part)
   {
   double     scaled[3];
   double     pt[3];
   double     fsfc;
   Partition *p;
   int        i;
   int        dim;
   HSFC_Data *d;
   int        err;
   char *yo = "Zoltan_HSFC_Point_Assign";

   ZOLTAN_TRACE_ENTER (zz, yo);
   d = (HSFC_Data *) zz->LB.Data_Structure;
   if (d == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL,
       "No Decomposition Data available; use KEEP_CUTS parameter.");

   if (d->Skip_Dimensions > 0){
     pt[1] = pt[2] = 0.0;
     pt[0] = d->Transformation[0][0]*x[0] +
             d->Transformation[0][1]*x[1] +
             d->Transformation[0][2]*x[2];
     if (d->Skip_Dimensions == 1){
       pt[1] = d->Transformation[1][0]*x[0] +
               d->Transformation[1][1]*x[1] +
               d->Transformation[1][2]*x[2];
     }
     dim = d->ndimension - d->Skip_Dimensions;
   }
   else{
     pt[0] = x[0];
     pt[1] = x[1];
     pt[2] = x[2];
     dim = d->ndimension;
   }

   /* Calculate scaled coordinates, calculate HSFC coordinate */
   for (i = 0; i < dim; i++)
      {
      scaled[i] = (pt[i] - d->bbox_lo[i]) / d->bbox_extent[i];
      if (scaled[i] < HSFC_EPSILON)         scaled[i] = HSFC_EPSILON;
      if (scaled[i] > 1.0 - HSFC_EPSILON)   scaled[i] = 1.0 - HSFC_EPSILON;
      }
   fsfc = d->fhsfc (zz, scaled);           /* Note, this is a function call */

   /* Find partition containing point and return its number */
   p = (Partition *) bsearch (&fsfc, d->final_partition, zz->LB.Num_Global_Parts,
    sizeof (Partition), Zoltan_HSFC_compare);
   if (p == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "programming error, shouldn't happen");
   if (part != NULL) {
      if (zz->LB.Remap)
         *part = zz->LB.Remap[p->index];
      else
         *part = p->index;
      }
   if (proc != NULL) {
      if (zz->LB.Remap) 
         *proc = Zoltan_LB_Part_To_Proc(zz, zz->LB.Remap[p->index], NULL);
      else
         *proc = Zoltan_LB_Part_To_Proc(zz, p->index, NULL);
      }
   err = ZOLTAN_OK;

End:
   ZOLTAN_TRACE_EXIT (zz, yo);
   return err;
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
