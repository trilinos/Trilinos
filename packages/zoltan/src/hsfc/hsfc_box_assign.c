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

/* For a detailed description of the following algorithm, please see the
   Developers Guide.  For instructions on use, please see the Users
   Guide.

   This is a temporary algorithm which suffers a limitation that will
   be removed in the future.  It is an approximation in the sense that
   every processor reported is actually in a partition that falls (in part)
   within the user's input box, but some processors may be missed!  The
   algorithm places a lattice of points withinn the user's specified box
   and determines the partition in which each box belongs. A partition may
   have volume that falls within the box but has no lattice point falling
   within.  Hence, it may be missed.  */



/****************************************************************************/

/* returns list of processors and partitions falling within user's box */
int Zoltan_HSFC_Box_Assign (
 ZZ *zz, double xlo, double ylo, double zlo,
         double xhi, double yhi, double zhi, 
         int *procs, int *proc_count, int *parts, int *part_count)
   {
   double     x[3], xintl[3], xinth[3];
   double     xdelta, ydelta, zdelta;
   int       *part_array = NULL;
   int       *proc_array = NULL;
   int        tmp_part;
   HSFC_Data *d;
   int        n, i, loop, offset;
   int        oldpartcount = 0;
   int        first_proc, last_proc;
   const int  MAXLOOP  = 15;
   int        err = ZOLTAN_OK;
   int prime[] = {3, 7, 17, 37, 79, 163, 331, 673, 1361, 2729, 5471, 10949,
    21911, 43853, 87719, 175447, 350899, 701819, 1403641, 2807303, 5614657,
    11229331, 22458671, 44917381};
   char *yo = "Zoltan_HSFC_Box_Assign";

   ZOLTAN_TRACE_ENTER (zz, yo);
   d = (HSFC_Data *) zz->LB.Data_Structure;
   if (d == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL,
          "No Decomposition Data available; use KEEP_CUTS parameter.");

   *part_count = *proc_count = 0;
   part_array = (int *) ZOLTAN_MALLOC ((zz->LB.Num_Global_Parts + zz->Num_Proc)
    * sizeof (int));
   if (part_array == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Memory error.");
   proc_array = part_array + zz->LB.Num_Global_Parts;
   memset (part_array, 0, (zz->LB.Num_Global_Parts + zz->Num_Proc)*sizeof (int));

   if (d->ndimension == 1) {
      Zoltan_HSFC_Point_Assign (zz, &xlo, NULL, &n);
      Zoltan_HSFC_Point_Assign (zz, &xhi, NULL, &loop);
      for (i = n; i <= loop; i++)
         part_array[i] = 1;
      goto fini;
      }
   
   /* set initial resolution at about equal work for each dimension */
   offset = (d->ndimension == 2) ? 2 : 1;

   /* determine intersection of bounding box and dropped box */
   if      (xhi < d->bbox_lo[0])  xintl[0] = xinth[0] = d->bbox_lo[0];
   else if (xlo > d->bbox_hi[0])  xintl[0] = xinth[0] = d->bbox_hi[0];
   else {
        xintl[0] = (d->bbox_lo[0] < xlo) ? xlo           : d->bbox_lo[0];
        xinth[0] = (d->bbox_hi[0] < xhi) ? d->bbox_hi[0] : xhi;
        }

   if      (yhi < d->bbox_lo[1])  xintl[1] = xinth[1] = d->bbox_lo[1];
   else if (ylo > d->bbox_hi[1])  xintl[1] = xinth[1] = d->bbox_hi[1];
   else {
        xintl[1] = (d->bbox_lo[1] < ylo) ? ylo           : d->bbox_lo[1];
        xinth[1] = (d->bbox_hi[1] < yhi) ? d->bbox_hi[1] : yhi;
        }

   if      (zhi < d->bbox_lo[2])  xintl[2] = xinth[2] = d->bbox_lo[2];
   else if (zlo > d->bbox_hi[2])  xintl[2] = xinth[2] = d->bbox_hi[2];
   else {
        xintl[2] = (d->bbox_lo[2] < zlo) ? zlo           : d->bbox_lo[2];
        xinth[2] = (d->bbox_hi[2] < zhi) ? d->bbox_hi[2] : zhi;
        }

   for (loop = 0; loop < MAXLOOP; loop++) {
      oldpartcount = *part_count;

      xdelta = (xinth[0] - xintl[0])/ (double) prime[loop + offset];
      ydelta = (xinth[1] - xintl[1])/ (double) prime[loop + offset];
      zdelta = (xinth[2] - xintl[2])/ (double) prime[loop + offset];

      if (xdelta < REFINEMENT_LIMIT)  xdelta = 1.0;
      if (ydelta < REFINEMENT_LIMIT)  ydelta = 1.0;
      if (zdelta < REFINEMENT_LIMIT)  zdelta = 1.0;

      /* create regular lattice in given box, then look up associated processors */
      if (d->ndimension == 3) {
        for     (x[0] = xintl[0]; x[0] <= xinth[0]; x[0] += xdelta)
          for   (x[1] = xintl[1]; x[1] <= xinth[1]; x[1] += ydelta)
            for (x[2] = xintl[2]; x[2] <= xinth[2]; x[2] += zdelta)
               if (Zoltan_HSFC_Point_Assign(zz, x, NULL, &tmp_part) == ZOLTAN_OK)
                   part_array[tmp_part] = 1;
        }

      else if (d->ndimension == 2) {
         for   (x[0] = xintl[0]; x[0] <= xinth[0]; x[0] += xdelta)
           for (x[1] = xintl[1]; x[1] <= xinth[1]; x[1] += ydelta)
              if (Zoltan_HSFC_Point_Assign (zz, x, NULL, &tmp_part) == ZOLTAN_OK)
                  part_array[tmp_part] = 1;
         }

      /* test for early exit */
      for (i = 0; i < zz->LB.Num_Global_Parts; i++)
         if (part_array[i] > 0)
            *part_count++;
      if (*part_count == oldpartcount)
         break;
      }

   /* move results to user supplied array(s) & hope it(they) is(are) big enough */
fini:
   *part_count = 0;
   if (parts) {
      *part_count = 0;
      for (i = 0; i < zz->LB.Num_Global_Parts; i++)
         if (part_array[i] > 0)
            parts[(*part_count)++] = i;
      }

   if (procs) {
      for (i = 0; i < zz->LB.Num_Global_Parts; i++)
         if (part_array[i] > 0) {
            first_proc = Zoltan_LB_Part_To_Proc(zz, i, NULL);
            proc_array[first_proc] = 1;
            if (!zz->LB.Single_Proc_Per_Part) {
               /* Part may be spread across multiple procs. Include them all. */
               int j;
               if (first_proc < zz->LB.Num_Global_Parts - 1)
                  last_proc = Zoltan_LB_Part_To_Proc(zz, first_proc+1, NULL);
               else
                  last_proc = zz->Num_Proc;

               for (j = first_proc + 1; j < last_proc; j++)
                  proc_array[j] = 1;
               }
            }

      *proc_count = 0;
      for (i = 0; i < zz->Num_Proc; i++)
         if (proc_array[i] > 0)
            procs[(*proc_count)++] = i;
      }

free:
   ZOLTAN_FREE (&part_array);
   ZOLTAN_TRACE_EXIT (zz, yo);
   return err;
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
