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

static void add_to_list(ZZ *, int, int, int *, int *, int, int) ;

/****************************************************************************/

/* returns list of processors and partitions falling within user's box */
int Zoltan_HSFC_Box_Assign (
 ZZ *zz, double xlo, double ylo, double zlo,
         double xhi, double yhi, double zhi, 
         int *procs, int *proc_count, int *parts, int *part_count)
   {
   double     x[3] ;
   double     xdelta, ydelta, zdelta ;
   int       *part_array = NULL ;
   int       *proc_array = NULL ;
   int        tmp_part ;
   int        tmp_proc ;
   HSFC_Data *d ;
   int        n, i, loop ;
   const int  NN       = 4 ;    /* determines lattice spacing */
   int        oldpartcount = 0 ;
   int        oldproccount = 0 ;
   const int  MAXLOOP  = 5 ;
   int        err ;
   char *yo = "Zoltan_HSFC_Box_Assign" ;
   int        include_procs = (procs != NULL);
   int        include_parts = (parts != NULL);
   int        size;

   ZOLTAN_TRACE_ENTER (zz, yo) ;
   d = (HSFC_Data *) zz->LB.Data_Structure ;
   if (d == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL,
          "No Decomposition Data available; use KEEP_CUTS parameter.");

   *part_count = *proc_count = 0;

   size = (include_parts ? zz->LB.Num_Global_Parts : 0) + 
          (include_procs ? zz->Num_Proc : 0) ;
   if (size == 0)
      goto free ;

   part_array = (int *) ZOLTAN_MALLOC (size * sizeof (int)) ;
   if (part_array == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Memory error.") ;

   proc_array = part_array + (include_parts ? zz->LB.Num_Global_Parts : 0) ;

   /* Clear processor array */
   memset (part_array, 0, size * sizeof (int)) ;  

   xdelta = (xlo == xhi) ? 1.0 : xhi - xlo ;
   ydelta = (ylo == yhi) ? 1.0 : yhi - ylo ;
   for (loop = 0, n = NN ; loop < MAXLOOP ; loop++, n = 2*n)
      {
      /* create regular lattice in given box, then look up associated processors */
      if (d->ndimension == 3)
         {
         zdelta = (zlo == zhi) ? 1.0L : zhi - zlo ;
         for       (x[0] = xlo ; x[0] <= xhi ; x[0] += xdelta/n)
            for    (x[1] = ylo ; x[1] <= yhi ; x[1] += ydelta/n)
               for (x[2] = zlo ; x[2] <= zhi ; x[2] += zdelta/n)
                  if (Zoltan_HSFC_Point_Assign (zz, x, &tmp_proc, &tmp_part) == 
                      ZOLTAN_OK) 
                     add_to_list(zz, include_procs, include_parts,
                                 proc_array, part_array, tmp_proc, tmp_part);
         }

      else if (d->ndimension == 2) 
         {
         for    (x[0] = xlo ; x[0] <= xhi ; x[0] += xdelta/n)
            for (x[1] = ylo ; x[1] <= yhi ; x[1] += ydelta/n)
               if (Zoltan_HSFC_Point_Assign (zz, x, &tmp_proc, &tmp_part) ==
                   ZOLTAN_OK)
                  add_to_list(zz, include_procs, include_parts,
                              proc_array, part_array, tmp_proc, tmp_part);
         }

      else if (d->ndimension == 1)
         {
         Zoltan_HSFC_Point_Assign (zz, &xlo, NULL, &n) ;
         Zoltan_HSFC_Point_Assign (zz, &xhi, NULL, &loop) ;
         for (i = n ; i <= loop ; i++) 
            {
            tmp_proc = Zoltan_LB_Part_To_Proc(zz, i, NULL);
            add_to_list(zz, include_procs, include_parts,
                        proc_array, part_array, tmp_proc, i);
            }
         }

      /* move results to user supplied array & hope it is big enough */
      if (include_parts)
         {
         *part_count = 0 ;
         for (i = 0 ; i < zz->LB.Num_Global_Parts ; i++)
            if (part_array[i] > 0)
               parts[(*part_count)++] = i ;
         }
     
      if (include_procs)
         {
         *proc_count = 0 ;
         for (i = 0 ; i < zz->Num_Proc ; i++)
            if (proc_array[i] > 0)
               procs[(*proc_count)++] = i ;
         }

      if (d->ndimension == 1 ||
          ((*part_count==oldpartcount || *part_count==zz->LB.Num_Global_Parts)
        && (*proc_count==oldproccount || *proc_count==zz->Num_Proc)))
         break ;

      oldpartcount = *part_count ;
      oldproccount = *proc_count ;
      }
   err = ZOLTAN_OK ;

free:
   ZOLTAN_FREE (&part_array) ;
   ZOLTAN_TRACE_EXIT (zz,yo) ;
   return err;
   }

/****************************************************************************/
static void add_to_list(
   ZZ *zz,
   int include_procs,    /* Flag:  add to processor list */
   int include_parts,    /* Flag:  add to partition list */
   int *proc_array,      /* Values to be incremented for procs to be added */
   int *part_array,      /* Values to be incremented for parts to be added */
   int proc,             /* First processor owning partition part. */
   int part              /* Partition to be added. */
)
   {
   /* Add a partition and ALL associated processors to the part_array and
    * proc_array lists. 
    * Note:  Special logic handles case where a partition is spread across
    * many processors (!zz->LB.Single_Proc_Per_Part).  
    * Want to include all processors in this case. 
    */
   int i;
   int last_proc;    /* Starting processor of next partition */

   if (include_parts) ++(part_array[part]) ;
   if (include_procs) 
      {
      ++(proc_array[proc]) ;
      if (!zz->LB.Single_Proc_Per_Part)
         {
         /* Partition may be spread across multiple procs.
            Include them all. */
         if (part < zz->LB.Num_Global_Parts - 1) 
            last_proc = Zoltan_LB_Part_To_Proc(zz, part+1, NULL);
         else
            last_proc = zz->Num_Proc;
                  
         for (i = proc+1 ; i < last_proc ; i++)
            ++(proc_array[i]) ;
         }
      }
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
