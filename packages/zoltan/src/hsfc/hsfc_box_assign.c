/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
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

/* For a detailed description of the following algorithm, please see the
   Developers Guide.  For instructions on use, please see the Users
   Guide.

   This is a temporary algorithm which suffers from two limitations that will
   be removed in the future.  First, it will not accept points outside the
   bounding box determined during the load balance phase. Second, it is
   approximate in the sense that every processor reported is actually in a
   partition that falls (in part) within the user's input box, but some
   processors may be missed!  The algorithm places a lattice of points within
   the user's specified box and determines the partition in which each box
   belongs. A partition may have volume that falls within the box but has no
   lattice point falling within.  Hence, it may be missed.  */


/* returns list of processors with partitions falling within user's box */
int Zoltan_HSFC_Box_Assign (
 ZZ *zz, double xlo, double ylo, double zlo,
         double xhi, double yhi, double zhi, int *procs,int *count)
   {
   double     x[3] ;
   double     xdelta, ydelta, zdelta ;
   int       *array = NULL ;
   int        proc ;
   HSFC_Data *d ;
   int        n, i, loop ;
   const int  NN       = 4 ;    /* determines lattice spacing */
   int        oldcount = 0 ;
   const int  MAXLOOP  = 5 ;
   int        err ;
   char *yo = "Zoltan_HSFC_Box_Assign" ;

   ZOLTAN_TRACE_ENTER (zz, yo) ;
   d = (HSFC_Data *) zz->LB.Data_Structure ;
   if (d == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, 
          "No Decomposition Data available; use KEEP_CUTS parameter.");

   if (xlo < d->bbox_lo[0])   xlo = d->bbox_lo[0] ;
   if (xhi > d->bbox_hi[0])   xhi = d->bbox_hi[0] ;
   xdelta = (xlo == xhi) ? 1.0L : xhi - ylo ;

   if (d->ndimension == 1)
     {
     x[0] = xlo ;
     if (ZOLTAN_OK != Zoltan_HSFC_Point_Assign (zz, x, &n))
        ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "invalid low coordinate") ;
     x[0] = xhi ;
     if (ZOLTAN_OK != Zoltan_HSFC_Point_Assign (zz, x, &loop))
        ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "invalid high coordinate") ;
     *count = 0 ;
     for (i = n ; i <= loop ; i++)
        procs[(*count)++] = i ;
     goto free ;
     }

   if (ylo < d->bbox_lo[1])   ylo = d->bbox_lo[1] ;
   if (yhi > d->bbox_hi[1])   yhi = d->bbox_hi[1] ;
   ydelta = (ylo == yhi) ? 1.0L : yhi - ylo ;

   array = (int *) ZOLTAN_MALLOC (zz->Num_Proc * sizeof (int)) ;
   if (array == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to malloc proc list") ;
   memset (array, 0, zz->Num_Proc * sizeof (int)) ;   /* Clear processor array */

   for (loop = 0, n = NN ; loop < MAXLOOP ; loop++, n = 2*n)
      {
      /* create regular lattice in given box, then look up associated processors */
      if (d->ndimension == 3)
        {
        if (zlo < d->bbox_lo[2])   zlo = d->bbox_lo[2] ;
        if (zhi > d->bbox_hi[2])   zhi = d->bbox_hi[2] ;
        zdelta = (zlo == zhi) ? 1.0L : zhi - zlo ;

         for       (x[0] = xlo ; x[0] <= xhi ; x[0] += xdelta/n)
            for    (x[1] = ylo ; x[1] <= yhi ; x[1] += ydelta/n)
               for (x[2] = zlo ; x[2] <= zhi ; x[2] += zdelta/n)
                  if (ZOLTAN_OK == Zoltan_HSFC_Point_Assign (zz, x, &proc))
                     ++array[proc] ;
        }

      if (d->ndimension == 2)
         for    (x[0] = xlo ; x[0] <= xhi ; x[0] += xdelta/n)
            for (x[1] = ylo ; x[1] <= yhi ; x[1] += ydelta/n)
               if (ZOLTAN_OK == Zoltan_HSFC_Point_Assign (zz, x, &proc))
                  ++array[proc] ;

      /* move results to user supplied array & hope it is big enough */
      *count = 0 ;
      for (i = 0 ; i < zz->Num_Proc ; i++)
         if (array[i] > 0)
            procs[(*count)++] = i ;

      if (*count == oldcount || *count == zz->Num_Proc)
         break ;
      oldcount = *count ;
      }
   err = ZOLTAN_OK ;

free:
   ZOLTAN_FREE (&array) ;
   ZOLTAN_TRACE_EXIT (zz,yo) ;
   return err;
   }
