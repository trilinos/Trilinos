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

/* Andy's Space Filling Curve (SFC) partitioning modified by Bob */

#include "hsfc.h"

int Zoltan_HSFC( /* Zoltan_HSFC - Load Balance: Hilbert Space Filling Curve */
 ZZ            *zz,
 int           *num_import,        /* set to -1 to indicate imports not used */
 ZOLTAN_ID_PTR *import_gids,       /* ignored */
 ZOLTAN_ID_PTR *import_lids,       /* ignored */
 int          **import_procs,      /* ignored */
 int           *num_export,        /* number of objects to migrate for load 
                                      balance */
 ZOLTAN_ID_PTR *export_gids,       /* list of Global IDs of migrating objects */
 ZOLTAN_ID_PTR *export_lids,       /* optional list of Local IDs of migrating 
                                      objects */
 int          **export_procs)      /* list of corresponding processor 
                                      destinations */
   {
   /* malloc'd arrays that need to be freed before completion */
   Dots      *dots = NULL ;
   ZOLTAN_ID_PTR  gids = NULL ;
   ZOLTAN_ID_PTR  lids = NULL ;
   Partition *grand_partition = NULL ;   /* fine partition of [0,1] */
   double    *partition       = NULL ;   /* course partition of [0,1] */
   float     *grand_weight    = NULL ;   /* "binned" weights on grand partition */
   float     *temp_weight     = NULL ;   /* local binning before MPI_Allreduce */
   float     *weights         = NULL ;   /* to get original local dots' w vec */
   float     *target          = NULL ;   /* desired weight in each partition */
   float     *work_fraction   = NULL ;   /* percentage of load in each partition */
   double    *delta           = NULL ;   /* refinement interval */
   HSFC_Data *d               = NULL ;   /* pointer to persistant data storage */

   /* other (non malloc'd) variables */
   int     ndots ;                 /* number of dots on this machine */
   int     pcount ;                /* number of partitions in grand partition */
   double  start_time, end_time ;  /* used to time execution */
   double  total_weight ;

   /* temporary variables, loop counters, etc. */
   int           i, j, k ;            /* loop counters */
   double        sum ;
   Partition    *p ;
   int           loop, max_loop = 1 ;        /* max_loop reset later */
   double        actual, desired;
   double        imbalance, corrected ;
   int           err ;
   unsigned int  key[3] ;    /* worst case (largest) dimension */
   char         *yo = "Zoltan_HSFC" ;

   /* begin program with trace, timing, and debug startup prints */
   ZOLTAN_TRACE_ENTER (zz, yo) ;
   start_time = Zoltan_Time(zz->Timer) ;
   if (zz->Debug_Level > 5)
      ZOLTAN_PRINT_INFO (zz->Proc, yo, "HSFC starting") ;
   *num_export = *num_import = -1 ;     /* in case of early error exit */

   /* allocate persistant storage required by box drop and point drop */
   Zoltan_HSFC_Free_Structure (zz) ;
   zz->LB.Data_Structure = (void *) ZOLTAN_MALLOC (sizeof (HSFC_Data)) ;
   if (zz->LB.Data_Structure == NULL)
      ZOLTAN_HSFC_ERROR(LB_MEMERR,"Error returned by malloc") ;
   d = (HSFC_Data *) zz->LB.Data_Structure ;

   /* Determine if dots have 2 or 3 dimensions */
   d->ndimension = zz->Get_Num_Geom(zz->Get_Num_Geom_Data, &err) ;
   if ((d->ndimension != 2 && d->ndimension != 3) || err)
      ZOLTAN_HSFC_ERROR(ZOLTAN_HSFC_ILLEGAL_DIMENSION,
       "Illegal dimension or error returned by Get_Num_Geom.") ;

   /* Determine how many local objects (dots) are on this processor */
   ndots = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &err) ;
   if (err || ndots < 0)
      ZOLTAN_HSFC_ERROR(ZOLTAN_HSFC_FATAL_ERROR, "Get_Num_Obj returned error.");

   /* allocate storage for dots and their corresponding gids, lids, weights */
   if (ndots > 0)
      {
      dots    = (Dots *) ZOLTAN_MALLOC (sizeof(Dots)  * ndots) ;
      gids    =          ZOLTAN_MALLOC_GID_ARRAY   (zz, ndots) ;
      lids    =          ZOLTAN_MALLOC_LID_ARRAY   (zz, ndots) ;
      weights = (float*) ZOLTAN_MALLOC (sizeof(float) * ndots
              * ((zz->Obj_Weight_Dim > 0) ? zz->Obj_Weight_Dim : 1)) ;

      if (weights == NULL || dots == NULL || (gids == NULL && zz->Num_GID)
       || (lids == NULL && zz->Num_LID))
          ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to ZOLTAN_MALLOC IDs.");

      /* obtain dot information: gids, lids, weights  */
      Zoltan_Get_Obj_List (zz, gids, lids, zz->Obj_Weight_Dim, weights, &err) ;
      if (err)
         ZOLTAN_HSFC_ERROR (ZOLTAN_HSFC_FATAL_ERROR,
          "Error returned by Zoltan_Get_Obj_List.") ;
      }

   /* set dot weights from object weights, if none set to one */
   if (zz->Obj_Weight_Dim > 1)
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Only one weight per Dot can be processed") ;
   for (i = 0 ; i < ndots ; i++)
      dots[i].weight = (zz->Obj_Weight_Dim < 1) ? 1.0
       : weights[i * zz->Obj_Weight_Dim] ;       /* use first dimension only */
   ZOLTAN_FREE (&weights) ;

   /* get dots' coordinates */
   for (i = 0 ; i < ndots ; i++)
      {
      zz->Get_Geom (zz->Get_Geom_Data, zz->Num_GID, zz->Num_LID,
       gids + i*zz->Num_GID, lids + i*zz->Num_LID, dots[i].x, &err) ;
      if (err != 0)
         ZOLTAN_HSFC_ERROR(ZOLTAN_HSFC_FATAL_ERROR, "Error in Get_Geom.") ;
      }

   /* allocate storage for partitions and "binned" weights */
   pcount = N * (zz->Num_Proc - 1) + 1 ;
   work_fraction   = (float *)   ZOLTAN_MALLOC (sizeof (float)     * zz->Num_Proc);
   target          = (float *)   ZOLTAN_MALLOC (sizeof (float)     * zz->Num_Proc);
   partition       = (double *)  ZOLTAN_MALLOC (sizeof (double)    * zz->Num_Proc);
   delta           = (double *)  ZOLTAN_MALLOC (sizeof (double)    * zz->Num_Proc);
   grand_partition = (Partition*)ZOLTAN_MALLOC (sizeof (Partition) * pcount) ;
   grand_weight    = (float *)   ZOLTAN_MALLOC (sizeof (float)     * (pcount+1)) ;
   temp_weight     = (float *)   ZOLTAN_MALLOC (sizeof (float)     * (pcount+1)) ;

   if (partition == NULL || grand_weight == NULL || grand_partition == NULL
    || temp_weight == NULL || target == NULL || work_fraction == NULL
    || delta == NULL)
      ZOLTAN_HSFC_ERROR (LB_MEMERR, "Malloc error for partitions, targets") ;

   {  /* Get bounding box, smallest coordinate aligned box containing all dots */
   double in[6], out[6], temp ;
   for (i = 0 ; i < d->ndimension ; i++)
      {
      in[i]                 = -HUGE_VAL ;
      in[i + d->ndimension] =  HUGE_VAL ;
      }
   for (i = 0 ; i < ndots ; i++)
      for (j = 0 ; j < d->ndimension ; j++)
        {
        if (dots[i].x[j]>in[j])               in[j]              =dots[i].x[j];
        if (dots[i].x[j]<in[j+d->ndimension]) in[j+d->ndimension]=dots[i].x[j];
        }

   for (i = 0 ; i < d->ndimension ; i++)  /* Andy's trick to save MPI call */
      in[i+d->ndimension] = -in[i+d->ndimension] ;
   MPI_Allreduce (in,out, 2*d->ndimension,MPI_DOUBLE,MPI_MAX,zz->Communicator);

   /* Enlarge bounding box to make points on faces become interior (Andy) */
   for (i = 0 ; i < d->ndimension ; i++)
      {
      temp = (out[i] - (-out[i+d->ndimension])) * ZOLTAN_HSFC_EPSILON * 0.5L ;
      d->bbox_hi[i]     =  out[i]               + temp ;
      d->bbox_lo[i]     = -out[i+d->ndimension] - temp ;
      d->bbox_extent[i] =  d->bbox_hi[i]        - d->bbox_lo[i] ;
      }
   }

   /* Compute & save HSFC value for all objects (dots) on this machine */
   d->fhsfc = (d->ndimension == 2) ? Zoltan_HSFC_fhsfc2d : Zoltan_HSFC_fhsfc3d ;
   for (i = 0 ; i < ndots ; i++)
      {
      double scaled[3] ;
      for (j = 0 ; j < d->ndimension ; j++)
         scaled[j] = (dots[i].x[j] - d->bbox_lo[j]) / d->bbox_extent[j] ;
      d->fhsfc (scaled, (int*)&TWO, key) ;
      dots[i].fsfc = convert_key (key) ;
      }

   /* Initalize work fraction (should be user specified vector in future) */
   for (i = 0 ; i < zz->Num_Proc ; i++)
      work_fraction[i] = 1.0L / (double) zz->Num_Proc ;

   /* Initialize grand partition to equally spaced intervals on [0,1] */
   for (i = 0 ; i < pcount ; i++)
      {
      grand_partition[i].index =  i ;
      grand_partition[i].r     = (i+1) / (double)pcount ;
      grand_partition[i].l     =  i    / (double)pcount ;
      }



   /* This loop is the real guts of the partitioning algorithm */
   for (loop = 0 ; loop < max_loop ; loop++)
      {
      /* globally bin weights for all dots using current grand partition */
      for (i = 0 ; i < pcount ; i++)
         temp_weight[i] = 0.0 ;
      temp_weight[pcount] = (float) ndots ;   /* holds total dots is system */

      for (i = 0 ; i < ndots ; i++)
         {
         p = bsearch (&dots[i].fsfc, grand_partition, pcount, sizeof(Partition),
          compare) ;
         if (p)
            temp_weight [p->index] += dots[i].weight ;
         else    /* programming error! */
            ZOLTAN_PRINT_INFO(zz->Proc,yo,"BSEARCH in LOOP ERROR") ;
         }
      MPI_Allreduce(temp_weight, grand_weight, pcount+1, MPI_FLOAT, MPI_SUM,
       zz->Communicator) ;

      /* Allocate ideal weights per partition, calculate max_loop -- once only */
      if (loop == 0)
         {
         max_loop = (int) ceil (log(grand_weight[pcount])/log(N-1)) ;   /* typically around 5 */
         total_weight = 0.0L ;
         for (i = 0 ; i < pcount ; i++)
            total_weight += grand_weight[i] ;
         for (i = 0 ; i < zz->Num_Proc ; i++)
            target[i] = work_fraction[i] * total_weight ;
         }

      /* debug: print target and grand partition values */
      if (zz->Debug_Level > 6 && zz->Proc == 0)
         {
         printf ("\n\nTarget %.2f, loop %d\n", target[0], loop) ;
         for (i = 0 ; i < pcount ; i++)
         printf ("Grand Partition %3d,   count %3.0f,    right= %0.6f\n",
          i, grand_weight[i], grand_partition[i].r) ;
         }

      /* test if convergence is satisfactory at nominal last loop */
      if (loop == max_loop-1)
         {
         j = 0 ;
         for (i = 0 ; i < pcount ; i++)
            if (grand_weight[i] > 0)
               j++ ;
         if (j > 2.1 * zz->Num_Proc)
            max_loop++ ;
         }

      /* don't need new partition from grand partition for last loop */
      if (loop == max_loop-1
       || loop > 2 * (int) ceil(log(grand_weight[pcount])/log(N-1)))
          break ;

      /* set new partition on [0,1] by summing grand_weights until targets met */
      /* k counts partitions and i counts grand partitions */
      {
      int adjust = 0 ;
      sum = 0.0L ;
      for (k = i = 0 ; i < pcount  &&  k < zz->Num_Proc ; )
         if (sum +  grand_weight[i] <= target[k])
             sum += grand_weight[i++] ;
         else
            {
            /* stop current summing to set new partition coordinate */
            partition[k] = grand_partition[i].l ;
            delta[k] = (grand_partition[i].r - grand_partition[i].l)/(N-1) ;

            /* If a "bin" holds multiple processors worth of dots: */
            if (k > 0 && partition[k-1] >= partition[k])
               {
               partition[k] = partition[k-1] + 0.5L
                * (grand_partition[i].r - partition[k-1]) ;
               delta[k] = delta[k-1] = delta[k-1] * 0.5L ;
               adjust++ ;                  /* need extra time for convergence */
               }
            grand_weight[i] -= (target[k] - sum) ; /* prevents systematic bias */
            sum = 0.0L ;
            k++ ;
            }
      if (adjust > 0)
         max_loop++ ;       /* need extra loop to refine overstuffed bins */
      }

      /* determine new grand partition by refining new partition boundary */
      for (k = i = 0 ; i < zz->Num_Proc - 1 ; i++)
         for (j = 0 ; j < N ; j++)
            {
            grand_partition[k].r = partition[i] + j * delta[i] ;
            grand_partition[k].l = (k == 0) ? 0.0 : grand_partition[k-1].r ;
            k++ ;
            }
      grand_partition[pcount-1].r = 1.0L ;    /* prevent roundoff error */
      grand_partition[pcount-1].l = grand_partition[pcount-2].r ;
      } /* end of loop */


   ZOLTAN_FREE (&temp_weight) ;
   ZOLTAN_FREE (&partition) ;
   ZOLTAN_FREE (&delta) ;
   ZOLTAN_FREE (&work_fraction) ;

   /* set new partition on [0,1] by summing grand_weights until targets met */
   /* k counts partitions and i counts grand partitions */
   d->final_partition= (Partition *)ZOLTAN_MALLOC(sizeof(Partition) * zz->Num_Proc);
   if (d->final_partition == NULL)
      ZOLTAN_HSFC_ERROR (LB_MEMERR, "Unable to malloc final_partition") ;

   sum       = 0.0L ;
   imbalance = 0.0L ;
   corrected = target[0] ;
   actual = desired = total_weight ;
   for (k = i = 0 ; i < pcount /* &&  k < zz->Num_Proc */ ; i++)
      if (sum +  grand_weight[i] <= corrected)
          sum += grand_weight[i] ;
      else
         {
         /* stop current summing to set new partition, greedy algorithm */
         if ((corrected - sum) >  (sum + grand_weight[i] - corrected))
            sum += grand_weight[i] ;
         else
            i-- ;

         d->final_partition[k].r = grand_partition[i].r ;
         d->final_partition[k].l = (k == 0) ? 0.0 : d->final_partition[k-1].r ;
         d->final_partition[k].index = k ;
         if (sum > (imbalance * target[k])  && target[k] > ZOLTAN_HSFC_EPSILON)
            imbalance = sum/target[k] ;

         /* correct target[]s for cumulative partitioning errors (Bruce H.) */
         actual   -= sum ;
         desired  -= target[k] ;
         corrected = target[++k] * (actual/desired) ;
         sum       = 0.0L ;                         /* preparing for next loop */
         }
   d->final_partition[zz->Num_Proc-1].r = 1.0 ;
   d->final_partition[zz->Num_Proc-1].l = d->final_partition[zz->Num_Proc-2].r ;
   d->final_partition[zz->Num_Proc-1].index = zz->Num_Proc-1 ;
   if (sum > (imbalance * target[k])  &&  target[k] > ZOLTAN_HSFC_EPSILON)
      imbalance = sum/target[k] ;

   /* Count the number of objects to export from this processor */
   *num_export = 0 ;
   for (i = 0 ; i < ndots ; i++)
      {
      dots[i].proc = -1 ;  /* in case of (programming) error returned by bsearch */
      p = bsearch (&dots[i].fsfc, d->final_partition, zz->Num_Proc,
       sizeof(Partition), compare) ;
      if (p)
         {
         dots[i].proc = p->index ;
         if (p->index != zz->Proc)
            ++(*num_export) ;
         }
      else
         ZOLTAN_PRINT_INFO(zz->Proc, yo, "BSEARCH ERROR"); /* programming error */
      }

   /* free stuff before next required allocations */
   ZOLTAN_FREE (&grand_weight) ;
   ZOLTAN_FREE (&grand_partition) ;
   ZOLTAN_FREE (&target) ;

   /* allocate storage for export information and fill in data */
   err = Zoltan_Special_Malloc (zz, (void**) export_gids, *num_export,
    ZOLTAN_SPECIAL_MALLOC_GID) ;
   if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to malloc global ids") ;

   if (zz->Num_LID > 0)
      {
      err = Zoltan_Special_Malloc (zz, (void**) export_lids, *num_export,
       ZOLTAN_SPECIAL_MALLOC_LID) ;
      if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
         ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to malloc local ids") ;
      }

   err = Zoltan_Special_Malloc (zz, (void**) export_procs, *num_export,
    ZOLTAN_SPECIAL_MALLOC_INT) ;
   if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to malloc proc list") ;

   /* Fill in export arrays */
   for (j = i = 0 ; i < ndots ; i++)
      if (dots[i].proc != zz->Proc  &&  dots[i].proc != -1)
         {
         *((*export_procs)+j) = dots[i].proc ;
         ZOLTAN_SET_GID(zz, *export_gids + j*zz->Num_GID, gids + i*zz->Num_GID);
         if (zz->Num_LID > 0)
            ZOLTAN_SET_LID(zz, *export_lids+j*zz->Num_LID, lids+i*zz->Num_LID);
         j++ ;
         }

   /* DEBUG: print useful information */
   end_time = Zoltan_Time(zz->Timer) ;
   if (zz->Debug_Level > 2  && zz->Proc == 0)
      printf ("TOTAL Processing Time is %.4f seconds\n", end_time-start_time) ;
   if (zz->Debug_Level > 5)
      printf ("<%d> Number of loops = %d\n", zz->Proc, loop+1) ;
   if (zz->Debug_Level > 5)
      printf ("<%d> export count %d, ndots %d\n", zz->Proc, *num_export, ndots) ;
   if (zz->Debug_Level > 7)
      for (i = 0 ; i < ndots ; i++)
         {
         printf ("<%d> GID: %u  LID: %u  Proc: %d  Weight %.1f  HSFC  %.4f\n", zz->Proc,
          gids[i*zz->Num_GID], lids[i*zz->Num_LID], dots[i].proc, dots[i].weight, dots[i].fsfc) ;
         printf ("PROC %d DOT %03u\n", p->index, gids[i]) ;
         }

#ifdef ZOLTAN_INTERNAL_TEST
   /* (optional) test point drop functionality */
   j = 0 ;
   for (i = 0 ; i < ndots ; i++)
      if (dots[i].proc != Zoltan_HSFC_Point_Drop (zz, dots[i].x))
         j++ ;
   printf ("<%d> Point Drop Test %s\n", zz->Proc, j ? "FAILED" : "PASSED");

   /* (optional) test box drop functionality */
   if (zz->Proc == 0)
      {
      int *proclist ;
      double hi[3], lo[3] ;

      hi[0] = d->bbox_hi[0] - ZOLTAN_HSFC_EPSILON ;
      hi[1] = d->bbox_hi[1] - ZOLTAN_HSFC_EPSILON ;
      hi[2] = d->bbox_hi[2] - ZOLTAN_HSFC_EPSILON ;
      lo[0] = d->bbox_lo[0] + ZOLTAN_HSFC_EPSILON ;
      lo[1] = d->bbox_lo[1] + ZOLTAN_HSFC_EPSILON ;
      lo[2] = d->bbox_lo[2] + ZOLTAN_HSFC_EPSILON ;

      proclist = (int *) ZOLTAN_MALLOC (sizeof(int) * zz->Num_Proc);
      Zoltan_HSFC_Box_Drop (zz, proclist, hi, lo, &i) ;
      printf ("Box Drop Test %s\n", (i==zz->Num_Proc) ? "PASSED" : "FAILED") ;
      ZOLTAN_FREE (&proclist) ;
      }
#endif

   /* done, do we keep data structure for box drop and point drop? */
   Zoltan_Bind_Param (HSFC_params, "KEEP_CUTS", (void *) &i) ;
   i = 0 ;
   Zoltan_Assign_Param_Vals (zz->Params, HSFC_params, zz->Debug_Level, zz->Proc,
    zz->Debug_Proc) ;
   if (i == 0)
      Zoltan_HSFC_Free_Structure (zz) ;

   /* really done now, now free dynamic storage and exit with return status */
   err = (imbalance > zz->LB.Imbalance_Tol) ? ZOLTAN_WARN : ZOLTAN_HSFC_OK ;
   if (zz->Proc == 0 && err == ZOLTAN_WARN && zz->Debug_Level > 0)
      ZOLTAN_PRINT_WARN (zz->Proc, yo, "HSFC: Imbalance exceeds user specificatation");

free:
   ZOLTAN_FREE (&dots) ;
   ZOLTAN_FREE (&gids) ;
   ZOLTAN_FREE (&lids) ;
   ZOLTAN_FREE (&grand_partition) ;
   ZOLTAN_FREE (&partition) ;
   ZOLTAN_FREE (&grand_weight) ;
   ZOLTAN_FREE (&temp_weight) ;
   ZOLTAN_FREE (&weights) ;
   ZOLTAN_FREE (&target) ;
   ZOLTAN_FREE (&work_fraction) ;
   ZOLTAN_FREE (&delta) ;

   if (zz->Debug_Level > 5)
      ZOLTAN_PRINT_INFO (zz->Proc, yo, "Exiting HSFC");
   ZOLTAN_TRACE_EXIT (zz, yo) ;
   return err ;
   }



/* Point drop for refinement after above partitioning */
int Zoltan_HSFC_Point_Drop (ZZ *zz, double *x)
   {
   double     scaled[3] ;
   double     fsfc ;
   Partition *p ;
   int        i ;
   HSFC_Data *d ;
   unsigned int key[3] ;
   char *yo = "Zoltan_HSFC_Point_Drop" ;

   d = (HSFC_Data *) zz->LB.Data_Structure ;
   if (d == NULL)
      return ZOLTAN_HSFC_FATAL_ERROR ;

   /* Insure that point is in bounding box */
   for (i = 0 ; i < d->ndimension ; i++)
      if ((x[i] > d->bbox_hi[i]) || (x[i] < d->bbox_lo[i]))
         return ZOLTAN_HSFC_FATAL_ERROR ;

   /* Calculate scaled coordinates, calculate HSFC coordinate */
   for (i = 0 ; i < d->ndimension ; i++)
      scaled[i] = (x[i] - d->bbox_lo[i]) / d->bbox_extent[i] ;
   d->fhsfc (scaled, (int *)&TWO, key) ;
   fsfc = convert_key (key) ;

   /* Find partition containing point and return its number */
   p = (Partition *) bsearch (&fsfc, d->final_partition, zz->Num_Proc,
    sizeof (Partition), compare) ;
   if (p)
      return p->index ;
   return ZOLTAN_HSFC_POINT_NOT_FOUND ;    /* programming error!! */
   }



/* temporary version of Box Drop, provides approximate answer */
void Zoltan_HSFC_Box_Drop (ZZ *zz, int *procs,double *hi,double *lo,int *count)
   {
   double     x[3] ;
   int       *array ;
   int        proc ;
   HSFC_Data *d ;
   int        n, i, loop ;
   const int  NN       = 4 ;
   int        oldcount = 0 ;
   const int  MAXLOOP  = 5 ;
   char *yo = "Zoltan_HSFC_Box_Drop" ;

   d = (HSFC_Data *) zz->LB.Data_Structure ;
   if (d == NULL)
      return ;

   array = (int *) ZOLTAN_MALLOC (zz->Num_Proc * sizeof (int)) ;
   if (array == NULL)
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Failed to malloc proc list") ;

   for (loop = 0, n = NN ; loop < MAXLOOP ; loop++, n = 2*n)
      {
      memset(array, 0, zz->Num_Proc * sizeof(int)) ; /* Clear processor array */

      /* create regular lattice in given box, then look up associated processors */
      if (d->ndimension == 3)
         for       (x[0] = lo[0] ; x[0] <= hi[0] ; x[0] += (hi[0]-lo[0])/n)
            for    (x[1] = lo[1] ; x[1] <= hi[1] ; x[1] += (hi[1]-lo[1])/n)
               for (x[2] = lo[2] ; x[2] <= hi[2] ; x[2] += (hi[2]-lo[2])/n)
                  {
                  proc = Zoltan_HSFC_Point_Drop (zz, x) ;
                  if (proc >= 0)
                     ++array[proc] ;
                  }

      if (d->ndimension == 2)
         for    (x[0] = lo[0] ; x[0] <= hi[0] ; x[0] += (hi[0]-lo[0])/n)
            for (x[1] = lo[1] ; x[1] <= hi[1] ; x[1] += (hi[1]-lo[1])/n)
               {
               proc = Zoltan_HSFC_Point_Drop (zz, x) ;
               if (proc >= 0)
                  ++array[proc] ;
               }

      /* move results to user supplied array & hope it is big enough */
      *count = 0 ;
      for (i = 0 ; i < zz->Num_Proc ; i++)
         if (array[i] > 0)
            procs[(*count)++] = i ;

      if (*count == oldcount || *count == zz->Num_Proc)
         break ;
      oldcount = *count ;
      }
   ZOLTAN_FREE (&array) ;
   }



/* routine for binary search, bsearch, to locate partition containing dot */
static int compare (const void *key, const void *arg)
   {
   if ( *(double *) key >  ((Partition *)arg)->r)
      return 1 ;
   if ( *(double *) key <= ((Partition *)arg)->l)
      return -1 ;
   return 0 ;
   }



/* convert 2 or 3 integer word HSFC value (key) to double precision */
static double convert_key (unsigned int *key)
   {
   return ldexp( (double)key[0], -32) + ldexp( (double)key[1], -64) ;
   }



/* free state memory needed by point & box drop routines */
void Zoltan_HSFC_Free_Structure (ZZ *zz)
   {
   HSFC_Data *data = (HSFC_Data *)zz->LB.Data_Structure ;
   if (data != NULL)
      ZOLTAN_FREE (&(data->final_partition)) ;
   ZOLTAN_FREE (&(zz->LB.Data_Structure)) ;
   }



/* function to read  "KEEP_CUTS" parameter: */
int Zoltan_HSFC_Set_Param (char *name, char *val)
   {
   int status, index ;
   PARAM_UTYPE result ;

   status = Zoltan_Check_Param (name, val, HSFC_params, &result, &index) ;
   return status ;
   }
