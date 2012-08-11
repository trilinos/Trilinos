/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "hsfc.h"
#include "zz_const.h"
#include "zz_util_const.h"

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

   for (i=0; i<d->ndimension; i++){
     pt[i] = x[i];  /* we don't want to change caller's "x" */
   }

   if (d->tran.Target_Dim > 0){   /* degenerate geometry */
     dim = d->tran.Target_Dim;
     Zoltan_Transform_Point(pt, d->tran.Transformation, d->tran.Permutation,
       d->ndimension, dim, pt);
   }
   else{
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
