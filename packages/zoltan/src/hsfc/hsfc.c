// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Andy Bauer's Space Filling Curve (SFC) partitioning modified by Bob Heaphy */

/* This code uses the Inverse Hilbert Space-Filling Curve to map 1-3 dimensional
   problems to the unit interval, [0,1].  It then partitions the unit interval
   into non-overlapping segments each containing equal weights of computational
   objects.  For details concerning the algorithm, please read the Developers
   Guide. For instructions on using this partitioning algorithm, please read the
   Users Guide.
*/


#include "hsfc.h"
#include "hsfc_params.h"
#include "zz_const.h"
#include <float.h>

/****************************************************************************/

static int partition_stats(ZZ *zz, int ndots,
                   Dots *dots, int *obj_sizes,
                   float *work_fraction, int *parts, int *new_parts);

/****************************************************************************/

/* Zoltan_HSFC - Main routine, Load Balance: Hilbert Space Filling Curve */
int Zoltan_HSFC(
 ZZ             *zz,
 float          *part_sizes,        /* Input: Array of size zz->Num_Global_Parts
                                      containing the percentage of work to be
                                      assigned to each partition. */
 int            *num_import,        /* -1 indicates imports not used */
 ZOLTAN_ID_PTR  *import_gids,       /* ignored */
 ZOLTAN_ID_PTR  *import_lids,       /* ignored */
 int           **import_procs,      /* ignored */
 int           **import_to_part,    /* ignored */
 int            *num_export,        /* number of objects to export to lb */
 ZOLTAN_ID_PTR  *export_gids,       /* list of Global IDs of exported objects */
 ZOLTAN_ID_PTR  *export_lids,       /* optional: Local IDs of exported objects */
 int           **export_procs,      /* corresponding processor destinations */
 int           **export_to_parts    /* partition assignments for export */
)
   {
   /* malloc'd arrays that need to be freed before completion */
   Dots      *dots            = NULL;
   ZOLTAN_ID_PTR  gids        = NULL;
   ZOLTAN_ID_PTR  lids        = NULL;
   int       *parts           = NULL;   /* initial partitions for objects */
   int       *new_proc        = NULL;   /* new processor assignments for dots */
   int       *new_part        = NULL;   /* new partition assignments for dots */
   Partition *grand_partition = NULL;   /* fine partition of [0,1] */
   double    *partition       = NULL;   /* course partition of [0,1] */
   double    *grand_weight    = NULL;   /* "binned" weights on grand partition */
   double    *temp_weight     = NULL;   /* local binning before MPI_Allreduce */
   float     *weights         = NULL;   /* to get original local dots' wt vec */
   float     *target          = NULL;   /* desired weight in each partition */
   float     *work_fraction   = NULL;   /* % of load in each partition */
   double    *delta           = NULL;   /* refinement interval */
   double    *tsum            = NULL;
   double    *geom_vec        = NULL;   /* Temporary coordinates array. */
   HSFC_Data *d               = NULL;   /* pointer to persistant data storage */

   /* other (non malloc'd) variables */
   int     ndots;                 /* number of dots on this machine */
   int     pcount;                /* number of partitions in grand partition */
   int     new_map;               /* flag indicating whether parts were
                                     remapped */
   double  start_time=0.0L, end_time=0.0L;  /* used to time execution */
   double  start_stat_time=0.0L, total_stat_time=0.0L;
   double  total_weight = 0.0;

   /* temporary variables, loop counters, etc. */
   int        tmp;
   int        i, j, k;            /* loop counters */
   double     sum;
   int        done, out_of_tolerance=0;    /* binary flags */
   Partition *p = NULL;
   int        loop = 0;
   double     actual, desired, correction;
   double     temp, in[6], out[6];
   int        err;
   int        final_output;
   int        param;
   int        idummy;
   double     ddummy;
   int        dim;
   char      *yo = "Zoltan_HSFC";

   /* begin program with trace, timing, and initializations */
   ZOLTAN_TRACE_ENTER (zz, yo);
   if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
      MPI_Barrier(zz->Communicator);
      start_time = Zoltan_Time(zz->Timer);
      }
   *num_export = *num_import = -1;             /* in case of early error exit */
   Zoltan_Bind_Param (HSFC_params, "FINAL_OUTPUT", (void*) &final_output);
   Zoltan_Bind_Param (HSFC_params, "KEEP_CUTS", (void*) &param);
   Zoltan_Bind_Param (HSFC_params, "REDUCE_DIMENSIONS", (void*) &idummy);
   Zoltan_Bind_Param (HSFC_params, "DEGENERATE_RATIO", (void*) &ddummy);
   param = idummy = final_output = 0;
   ddummy = 0.0;
   Zoltan_Assign_Param_Vals (zz->Params, HSFC_params, zz->Debug_Level, zz->Proc,
    zz->Debug_Proc);

   if (sizeof (int) != 4) {
     ZOLTAN_HSFC_ERROR(ZOLTAN_FATAL,
                       "HSFC implemented only for 32-bit integers");
   }

   /* allocate persistent storage required by box assign and point assign */

   if (zz->LB.Data_Structure == NULL){
     zz->LB.Data_Structure = (void*) ZOLTAN_MALLOC (sizeof (HSFC_Data));
     if (zz->LB.Data_Structure == NULL){
       ZOLTAN_HSFC_ERROR(ZOLTAN_MEMERR, "Error returned by malloc");
     }

     d = (HSFC_Data*) zz->LB.Data_Structure;
     memset ((void*)d, 0, sizeof (HSFC_Data));
     Zoltan_Initialize_Transformation(&(d->tran));
   }
   else{
     d = (HSFC_Data*) zz->LB.Data_Structure;
     ZOLTAN_FREE (&d->final_partition);
   }

   /* obtain dot information: gids, lids, weights  */
   err = Zoltan_Get_Obj_List (zz, &ndots, &gids, &lids, zz->Obj_Weight_Dim,
    &weights, &parts);
   if (err)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "Error in Zoltan_Get_Obj_List.");

   MPI_Allreduce(&ndots, &i, 1, MPI_INT, MPI_MAX, zz->Communicator);
   if (i < 1){
     if (zz->Proc == 0)
       ZOLTAN_PRINT_WARN(zz->Proc, yo, "No objects to partition")
     *num_export = 0;
     goto EndReporting;
   }

   /* allocate storage for dots and their corresponding gids, lids, weights */
   if (ndots > 0) {
      dots = (Dots*) ZOLTAN_MALLOC (sizeof(Dots) * ndots);
      if (dots == NULL)
          ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to ZOLTAN_MALLOC dots.");
      }

   /* set dot weights from object weights, if none set to one */
   if (zz->Obj_Weight_Dim > 1)
      ZOLTAN_PRINT_WARN(zz->Proc, yo,"Only one weight per GID can be processed");
   if (zz->Obj_Weight_Dim < 1)
      for (i = 0; i < ndots; i++)
         dots[i].weight = DEFAULT_WEIGHT;
   else
      for (i = 0; i < ndots; i++)
         dots[i].weight = weights [i * zz->Obj_Weight_Dim]; /* 1 dimension only */
   ZOLTAN_FREE (&weights);

   /* get dots' coordinates */
   err = Zoltan_Get_Coordinates(zz, ndots, gids, lids, &(d->ndimension),
                                &geom_vec);

   if (err != 0)
      ZOLTAN_HSFC_ERROR(ZOLTAN_FATAL, "Error in Zoltan_Get_Coordinates.");

   if (d->tran.Target_Dim > 0){  /* degenerate geometry */
     dim = d->tran.Target_Dim;
   }
   else{
     dim = d->ndimension;
   }

   for (i = 0, tmp=0; i < ndots; i++, tmp += d->ndimension) {
      for (j = 0; j < d->ndimension; j++){
         dots[i].x[j] = geom_vec[tmp + j];
      }
   }

   ZOLTAN_FREE(&geom_vec);

   ZOLTAN_TRACE_DETAIL (zz, yo, "Obtained dot information");

   /* allocate storage for partitions and "binned" weights */
   pcount = N * (zz->LB.Num_Global_Parts - 1) + 1;
   target       =(float*) ZOLTAN_MALLOC(sizeof(float) *zz->LB.Num_Global_Parts);
   partition    =(double*)ZOLTAN_MALLOC(sizeof(double)*zz->LB.Num_Global_Parts);
   delta        =(double*)ZOLTAN_MALLOC(sizeof(double)*zz->LB.Num_Global_Parts);
   tsum         =(double*)ZOLTAN_MALLOC(sizeof(double)*zz->LB.Num_Global_Parts);
   grand_weight =(double*)ZOLTAN_MALLOC(sizeof(double) * pcount*3);/*max,min,sum*/
   temp_weight  =(double*)ZOLTAN_MALLOC(sizeof(double) * pcount*3);/*max,min,sum*/
   grand_partition = (Partition*)ZOLTAN_MALLOC (sizeof (Partition) * pcount);

   if (target == NULL  ||  partition == NULL  ||  delta == NULL || tsum == NULL
    || grand_weight == NULL  ||  temp_weight == NULL  ||  grand_partition == NULL)
       ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Malloc error for partitions, targets");

   if (zz->Obj_Weight_Dim <= 1)
      work_fraction = part_sizes;
   else {
      work_fraction =(float*)ZOLTAN_MALLOC(sizeof(float)*zz->LB.Num_Global_Parts);
      if (work_fraction == NULL)
         ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Malloc error for work_fraction");

      /* HSFC supports only the first weight; select the appropriate part_sizes */
      /* entries for the first weight. */
      for (i = 0; i < zz->LB.Num_Global_Parts; i++)
         work_fraction[i] = part_sizes[i*zz->Obj_Weight_Dim];
      }

   /* Get bounding box, smallest coordinate aligned box containing all dots */
   for (i = 0; i < dim; i++)
      out[i] = in[i] = -HUGE_VAL;
   for (i = dim; i < 2*dim; i++)
      out[i] = in[i] =  HUGE_VAL;
   for (i = 0; i < ndots; i++) {
     for (j = 0; j < dim; j++) {
       /* get maximum bound box coordinates: */
       if(dots[i].x[j]>in[j]) in[j]=dots[i].x[j];
       }
     for (j = dim; j < 2*dim; j++) {
       /* get minimum bound box coordinates: */
       if(dots[i].x[j-dim]<in[j]) in[j]=dots[i].x[j-dim];
       }
     }
   err = MPI_Allreduce(in,out,dim,MPI_DOUBLE,MPI_MAX,zz->Communicator);
   err = MPI_Allreduce(in+dim,out+dim,dim,MPI_DOUBLE,MPI_MIN,zz->Communicator);
   if (err != MPI_SUCCESS)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "Bounding box MPI_Allreduce error");

   /* Enlarge bounding box to make points on faces become interior (Andy) */
   for (i = 0; i < dim; i++) {
      temp = (out[i] - out[i+dim]) * HSFC_EPSILON;
      d->bbox_hi[i] = out[i] + temp;
      d->bbox_lo[i] = out[i+dim] - temp;
      d->bbox_extent[i] = d->bbox_hi[i] - d->bbox_lo[i];
      if (d->bbox_extent[i] == 0.0)
          d->bbox_extent[i]  = 1.0; /* degenerate axis, avoid divide by zero */
      }
   ZOLTAN_TRACE_DETAIL (zz, yo, "Determined problem bounding box");

   /* Save appropriate HSFC function for later use below and in point_assign() */
   if      (dim== 1)  d->fhsfc = Zoltan_HSFC_InvHilbert1d;
   else if (dim== 2)  d->fhsfc = Zoltan_HSFC_InvHilbert2d;
   else if (dim== 3)  d->fhsfc = Zoltan_HSFC_InvHilbert3d;

   /* Scale coordinates to bounding box, compute HSFC */
   for (i = 0; i < ndots; i++) {
      for (j = 0; j < dim; j++)
         out[j] = (dots[i].x[j] - d->bbox_lo[j]) / d->bbox_extent[j];
      dots[i].fsfc = d->fhsfc (zz, out);      /* Note, this is a function call */
      }

   /* Initialize grand partition to equally spaced intervals on [0,1] */
   for (i = 0; i < pcount; i++) {
      grand_partition[i].index =  i;
      grand_partition[i].r     = (i+1) / (double)pcount;
      grand_partition[i].l     =  i    / (double)pcount;
      }
   grand_partition[pcount-1].r += (2.0 * FLT_EPSILON) ;

   ZOLTAN_TRACE_DETAIL (zz, yo, "About to enter main loop\n");


   /* This loop is the real guts of the partitioning algorithm */
   for (loop = 0; loop < MAX_LOOPS; loop++) {
      /* initialize bins, DEFAULT_BIN_MAX is less than any possible max,... */
      for (i = 0;        i <   pcount; i++)
         grand_weight[i] = temp_weight[i] = 0.0; /* SUM */
      for (i =   pcount; i < 2*pcount; i++)
         grand_weight[i] = temp_weight[i] = DEFAULT_BIN_MAX;
      for (i = 2*pcount; i < 3*pcount; i++)
         grand_weight[i] = temp_weight[i] = DEFAULT_BIN_MIN;

      /* bin weights, max, min for all dots using current grand partition */
      for (i = 0; i < ndots; i++) {
         if (loop > 0
          && dots[i].fsfc <   grand_partition[dots[i].part].r
          && dots[i].fsfc >=  grand_partition[dots[i].part].l)
             p = &grand_partition[dots[i].part];   /* short cut if unchanged */
         else
             p = (Partition*) bsearch (&dots[i].fsfc, grand_partition, pcount,
              sizeof(Partition), Zoltan_HSFC_compare);
         if (p == NULL)            /* programming error if NULL */
            ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "BSEARCH RETURNED ERROR");

         dots[i].part = p->index;
         temp_weight [p->index] += dots[i].weight;             /* local wt sum */
         if (dots[i].fsfc > temp_weight[p->index + pcount])
            temp_weight [p->index + pcount] = dots[i].fsfc;    /*local indx max*/
         if (dots[i].fsfc < temp_weight[p->index + 2 * pcount])
            temp_weight [p->index + pcount * 2] = dots[i].fsfc;/*local indx min*/
         }
      err = MPI_Allreduce(temp_weight, grand_weight, pcount,
                          MPI_DOUBLE, MPI_SUM, zz->Communicator);
      err = MPI_Allreduce(temp_weight+pcount, grand_weight+pcount, pcount,
                          MPI_DOUBLE, MPI_MAX, zz->Communicator);
      err = MPI_Allreduce(temp_weight+2*pcount, grand_weight+2*pcount, pcount,
                          MPI_DOUBLE, MPI_MIN, zz->Communicator);
      if (err != MPI_SUCCESS)
         ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "MPI_Allreduce returned error");
      ZOLTAN_TRACE_DETAIL (zz, yo, "Complete main loop MPI_Allreduce");

      /* Allocate ideal weights per partition -- once is sufficient */
      if (loop == 0) {
         total_weight = 0.0;
         for (i = 0; i < pcount; i++)
            total_weight += grand_weight[i];
         for (i = 0; i < zz->LB.Num_Global_Parts; i++)
            target[i] = work_fraction[i] * total_weight;
         }

      /* debug: print target and grand partition values */
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL && zz->Proc == 0) {
         printf ("\n\nTarget %.2f, loop %d\n", target[0], loop);
         for (i = 0; i < pcount; i++)
            printf ("Grand Partition %3d,  weight %3.0f, left = %0.6f,   right"
             "= %0.6f,  delta = %19.14f\n", i, grand_weight[i],
             grand_partition[i].l, grand_partition[i].r, grand_weight[i+pcount]
             - grand_weight[i+2*pcount]);
         }

      /* create new partition by summing contiguous bins l to r until target */
      for (i = 0; i < pcount; i++)    /* since grand weight modified below */
         temp_weight[i] = grand_weight[i];
      for (k = 0; k < zz->LB.Num_Global_Parts; k++) {
         partition[k] = 1.0;
         delta[k]     = 0.0;
         }
      done = 1;                       /* flag, 1 means done is true */
      sum  = 0.0;
      for (k = i = 0; i < pcount && k < zz->LB.Num_Global_Parts - 1; )
         if (sum +  grand_weight[i] <= target[k])
             sum += grand_weight[i++];
         else {
            /* stop current summing to set new partition coordinate */
            partition[k] = grand_partition[i].l;
            delta[k] = (grand_partition[i].r - grand_partition[i].l)/(N-1);

            /* If a "bin" holds multiple processors worth of dots: */
            if (k > 0  &&  partition[k-1] >= partition[k]) {
               partition[k] = 0.5 * (grand_partition[i].r - partition[k-1])
                + partition[k-1];
               delta[k] = delta[k-1] = delta[k-1] * 0.5;
               }

            /* max - min test if current bin is capable of further refinement */
            temp = grand_weight[i+pcount] - grand_weight[i+2*pcount];
            if (i < pcount  &&  temp > REFINEMENT_LIMIT)
               done = 0;         /* not done, at least this bin is refinable */

            /* get ready for next loop */
            grand_weight[i] -= (target[k] - sum); /* prevents systematic bias */
            sum = 0.0;
            k++;
            }
      if (done)    /* this is the actual loop stopping criterion */
         break;    /* no bin with "cut" is refinable, time to quit */

      /* determine new grand partition by refining new partition boundary */
      for (k = i = 0; i < zz->LB.Num_Global_Parts - 1; i++)
         for (j = 0; j < N; j++, k++)
            grand_partition[k+1].l = grand_partition[k].r
             = partition[i] + j * delta[i];
      grand_partition[0].l        = 0.0;
      grand_partition[pcount-1].r = 1.0 + (2.0 * FLT_EPSILON) ;
      } /* end of loop */

   ZOLTAN_TRACE_DETAIL (zz, yo, "Exited main loop");
   Zoltan_Multifree (__FILE__, __LINE__, 3, &grand_weight, &partition, &delta);
   if (zz->Obj_Weight_Dim > 1)
      ZOLTAN_FREE (&work_fraction);

    d->final_partition = (Partition*) ZOLTAN_MALLOC(sizeof(Partition)
    * zz->LB.Num_Global_Parts);
   if (d->final_partition == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Unable to malloc final_partition");

   /* initializations required to start loop below */
   d->final_partition[0].l = 0.0;
   for (k = 0; k < zz->LB.Num_Global_Parts; k++) {
      tsum[k] = 0.0;
      d->final_partition[k].index = k;
      }
   actual = desired = total_weight;
   i = k = 0;
   correction = 1.0;

   /* set new partition on [0,1] by summing grand_weights until targets met */
   /* k counts partitions and i counts grand partitions */
   while (1)  {
      if (k >= zz->LB.Num_Global_Parts  ||  i >= pcount)
         break;

      /* case:  current partition should remain empty */
      if (target[k] == 0.0)  {
         d->final_partition[k].r = grand_partition[i].l;
         k++;
         if (k < zz->LB.Num_Global_Parts)
            d->final_partition[k].l = grand_partition[i].l;
         continue;
         }

      /* case: current bin weights fit into current partition */
      temp = correction * target[k];
      if (tsum[k] + temp_weight[i] <= temp)  {
         tsum[k] += temp_weight[i];
         i++;
         continue;
         }

      /* case: current bin weights overfill current partition */
      if (temp - tsum[k] > tsum[k] + temp_weight[i] - temp)  {
         tsum[k] += temp_weight[i];
         actual  -= tsum[k];
         desired -= target[k];
         d->final_partition[k].r = grand_partition[i].r;
         k++;
         if (k < zz->LB.Num_Global_Parts)
            d->final_partition[k].l = grand_partition[i].r;
         i++;
         }
      else    { /* don't include current bin weight in current partition */
         actual  -= tsum[k];
         desired -= target[k];
         d->final_partition[k].r = grand_partition[i].l;
         k++;
         if (k < zz->LB.Num_Global_Parts)
            d->final_partition[k].l = grand_partition[i].l;
         }

      /* correct target[]s for cumulative partitioning errors (Bruce H.) */
      correction = ((desired == 0) ? 1.0 : actual/desired);
      }
   /* check if we didn't close out last sum, fix right boundary if needed */
   if (i == pcount && k < zz->LB.Num_Global_Parts)
     {
     d->final_partition[k].r = 1.0 + (2.0 * FLT_EPSILON) ;
     k++;
     }

   /* if last partition(s) is empty, loop stops w/o setting final_partition(s) */
   for (i = k; i < zz->LB.Num_Global_Parts; i++)  {
      d->final_partition[i].r = 1.0 + (2.0 * FLT_EPSILON);
      d->final_partition[i].l = 1.0 + (2.0 * FLT_EPSILON);
      }

   out_of_tolerance = 0;
   for (k = 0; k < zz->LB.Num_Global_Parts; k++)
      if (tsum[k] > target[k] * zz->LB.Imbalance_Tol[0])
         out_of_tolerance = 1;
   ZOLTAN_TRACE_DETAIL (zz, yo, "Determined final partition");

   new_part = (int *) ZOLTAN_MALLOC(2 * ndots * sizeof(int));
   if (ndots && !new_part)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Memory error.");
   new_proc = new_part + ndots;

   /* Set the final part number for all dots */
   for (i = 0; i < ndots; i++) {
      j = dots[i].part / N; /* grand_partition is N times final_partition parts*/
      if (dots[i].fsfc <  d->final_partition[j].r
       && dots[i].fsfc >= d->final_partition[j].l)
          p = d->final_partition + j;
      else
          p = (Partition*) bsearch (&dots[i].fsfc, d->final_partition,
           zz->LB.Num_Global_Parts, sizeof(Partition), Zoltan_HSFC_compare);
      if (p == NULL)
          ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "BSEARCH RETURNED ERROR");

      dots[i].part = p->index;
      new_part[i] = dots[i].part;
      new_proc[i] = Zoltan_LB_Part_To_Proc(zz, p->index,
        gids + i * zz->Num_GID);
      }

   /* Remap partitions to reduce data movement. */
   if (zz->LB.Remap_Flag) {
      err = Zoltan_LB_Remap(zz, &new_map, ndots, new_proc, parts, new_part, 1);
      if (err < 0)
         ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL,"Error returned from Zoltan_LB_Remap");
      }

   total_stat_time = 0.0;
   if (final_output){
     int *objSizes = NULL;
     if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
        start_stat_time = Zoltan_Time(zz->Timer);
        }
     if ((ndots > 0) && ((zz->Get_Obj_Size_Multi) || (zz->Get_Obj_Size))){

       objSizes = (int *) ZOLTAN_MALLOC(ndots * sizeof(int));
       if (!objSizes){
         ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to ZOLTAN_MALLOC sizes.");
       }

       if (zz->Get_Obj_Size_Multi) {
         zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data,
                                zz->Num_GID, zz->Num_LID, ndots,
                                gids, lids, objSizes, &err);
         if (err < 0) {
           ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "Error returned from "
                           "ZOLTAN_OBJ_SIZE_MULTI function.");
         }
       }
       else if (zz->Get_Obj_Size) {
         for (i = 0; i < ndots; i++) {
           ZOLTAN_ID_PTR lid = (zz->Num_LID ? &(lids[i*zz->Num_LID]):NULL);
           objSizes[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data,
                                          zz->Num_GID, zz->Num_LID,
                                            &(gids[i*zz->Num_GID]),
                                            lid, &err);
           if (err < 0) {
             ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "Error returned from "
                           "ZOLTAN_OBJ_SIZE function.");
           }
         }
       }
     }
     err = partition_stats(zz, ndots, dots, objSizes, work_fraction,
                           parts, new_part);

     if (err != ZOLTAN_OK){
       ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "statistics");
     }
     if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
        total_stat_time = Zoltan_Time(zz->Timer) - start_stat_time;
     }

     ZOLTAN_FREE(&objSizes);
   }

   /* free stuff before next required allocations */
   Zoltan_Multifree (__FILE__, __LINE__, 3, &temp_weight, &grand_partition,
    &target);

   /* Count the number of objects to export from this processor */
   *num_export = 0;
   for (i = 0; i < ndots; i++) {
      if (new_part[i] != parts[i]  ||  zz->Proc != new_proc[i])
         ++(*num_export);
      }

EndReporting:

   if (!zz->LB.Return_Lists)
      *num_export = -1;
   else if (zz->LB.Return_Lists == ZOLTAN_LB_CANDIDATE_LISTS) {
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "Candidate Lists not supported in HSFC;"
                                       "change RETURN_LISTS parameter");
   }
   else if (*num_export > 0) {
      /* allocate storage for export information and fill in data */
      if (!Zoltan_Special_Malloc (zz, (void**) export_gids, *num_export,
                                  ZOLTAN_SPECIAL_MALLOC_GID))
         ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to malloc global ids");

      if (zz->Num_LID > 0) {
         if (!Zoltan_Special_Malloc (zz, (void**) export_lids, *num_export,
                                     ZOLTAN_SPECIAL_MALLOC_LID)){
            ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to malloc local ids");
            }
         }

      if (!Zoltan_Special_Malloc (zz, (void**) export_procs, *num_export,
                                  ZOLTAN_SPECIAL_MALLOC_INT))
         ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to malloc export proc list");

      if (!Zoltan_Special_Malloc (zz, (void**) export_to_parts, *num_export,
                                  ZOLTAN_SPECIAL_MALLOC_INT))
         ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Failed to malloc export part list");

      /* Fill in export arrays */
      for (j = i = 0; i < ndots; i++) {
         if (new_part[i] != parts[i] || zz->Proc != new_proc[i]) {
            *((*export_procs)   +j) = new_proc[i];
            *((*export_to_parts)+j) = new_part[i];
            if (zz->Num_GID > 0)
               ZOLTAN_SET_GID(zz, *export_gids+j*zz->Num_GID, gids+i*zz->Num_GID);
            if (zz->Num_LID > 0)
               ZOLTAN_SET_LID(zz, *export_lids+j*zz->Num_LID, lids+i*zz->Num_LID);
            j++;
            }
         }
      }
    ZOLTAN_TRACE_DETAIL (zz, yo, "Filled in export information");

   /* DEBUG: print useful information */
   if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL  &&  zz->Proc == 0)
      printf ("<%d> Number of loops = %d\n", zz->Proc, loop + 1);
   if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST  && zz->Proc == 0)
      printf ("<%d> export count %d, ndots %d\n", zz->Proc, *num_export, ndots);
   if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST)
      for (i = 0; i < ndots; i++) {
         printf ("<%d> GID: ", zz->Proc);
         ZOLTAN_PRINT_GID (zz, &gids[i*zz->Num_GID]);
         if (lids || zz->Num_LID > 0) {
            printf ("  LID: ");
            ZOLTAN_PRINT_LID (zz, &lids[i*zz->Num_LID]);
            }
         printf ("  Part: %d  Weight %.1f  HSFC  %.6f\n", new_part[i],
          dots[i].weight, dots[i].fsfc);
         printf ("PROC %d DOT " ZOLTAN_ID_SPEC "\n", p->index, gids[i]);
         }

   ZOLTAN_FREE(&new_part);

   /* done, do we keep data structure for box drop and point drop? */

   if (!param &&                    /* we don't need partitions */
       (d->tran.Target_Dim < 0)){   /* we don't need transformation */
     Zoltan_HSFC_Free_Structure (zz);
   }

   /* really done now, now free dynamic storage and exit with return status */
   err = ((out_of_tolerance) ? ZOLTAN_WARN : ZOLTAN_OK);
   if (zz->Proc == 0 && err == ZOLTAN_WARN && zz->Debug_Level >ZOLTAN_DEBUG_NONE)
      ZOLTAN_PRINT_WARN (zz->Proc, yo, "HSFC: Imbalance exceeds user limit");

End:

   if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      Zoltan_Special_Free(zz, (void **)export_gids, ZOLTAN_SPECIAL_MALLOC_GID);
      Zoltan_Special_Free(zz, (void **)export_lids, ZOLTAN_SPECIAL_MALLOC_LID);
      Zoltan_Special_Free(zz, (void **)export_procs, ZOLTAN_SPECIAL_MALLOC_INT);
      Zoltan_Special_Free(zz, (void **)export_to_parts,
                          ZOLTAN_SPECIAL_MALLOC_INT);
      }
   Zoltan_Multifree (__FILE__, __LINE__, 12, &dots, &gids, &lids, &partition,
    &grand_partition, &grand_weight, &temp_weight, &weights, &target, &delta,
    &parts, &tsum);
   if (zz->Obj_Weight_Dim > 1)
      ZOLTAN_FREE (&work_fraction);

   if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
      MPI_Barrier(zz->Communicator);
      end_time = Zoltan_Time(zz->Timer);
      if (zz->Debug_Proc == zz->Proc)
         printf ("HSFC Processing Time is %.6f seconds\n",
                  end_time-start_time-total_stat_time);
      }

   ZOLTAN_TRACE_EXIT (zz, yo);
   return err;
   }



/* routine for bsearch locating the partition segment holding key */
int Zoltan_HSFC_compare (const void *key, const void *arg)
   {
   if ( *(double*) key >=  ((Partition*) arg)->r)  return  1;
   if ( *(double*) key <   ((Partition*) arg)->l)  return -1;

   return 0;     /* key in arg interval [l,r) */
   }



/* free state memory needed by point & box drop routines */
void Zoltan_HSFC_Free_Structure (ZZ *zz)
   {
   HSFC_Data *data = (HSFC_Data*) zz->LB.Data_Structure;
   if (data != NULL)
      ZOLTAN_FREE (&data->final_partition);
   ZOLTAN_FREE (&zz->LB.Data_Structure);
   }

/* Copy an hsfc data structure */

int Zoltan_HSFC_Copy_Structure(ZZ *toZZ, ZZ const *fromZZ)
{
  char *yo = "Zoltan_HSFC_Copy_Structure";
  int len;
  HSFC_Data *to;
  HSFC_Data const *from;

  Zoltan_HSFC_Free_Structure(toZZ);
  from = (HSFC_Data const *)fromZZ->LB.Data_Structure;

  if (!from){
    return ZOLTAN_OK;
  }

  to = (HSFC_Data *)ZOLTAN_MALLOC(sizeof(HSFC_Data));
  if (!to){
    ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  toZZ->LB.Data_Structure = (void *)to;

  *to = *from;

  if (from->final_partition){
    len = sizeof(Partition) * fromZZ->LB.Num_Global_Parts;

    to->final_partition = (Partition*) ZOLTAN_MALLOC(len);

    if (!to->final_partition){
      Zoltan_HSFC_Free_Structure(toZZ);
      ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }

    memcpy(to->final_partition, from->final_partition, len);
  }

  return ZOLTAN_OK;
}

/* function to read HSFC parameters: */
int Zoltan_HSFC_Set_Param (char *name, char *val)
   {
   int index;
   PARAM_UTYPE result;

   return Zoltan_Check_Param (name, val, HSFC_params, &result, &index);
   }


/****************************************************************************/
static int partition_stats(ZZ *zz, int ndots,
                   Dots *dots, int *obj_sizes,
                   float *work_fraction, int *parts, int *new_parts)
{
int wgtDim = zz->Obj_Weight_Dim;
int proc = zz->Proc;
int i, j, max, numParts;
double move, gmove, bal, max_imbal, ib;
double lwtot, gwtot;
double *lpartWgt = NULL;
double *gpartWgt = NULL;

  lwtot = 0.0;
  for (i = 0, max=0; i < ndots; i++){
    if (new_parts[i] > max) max = new_parts[i];
    if (wgtDim){
      lwtot += dots[i * wgtDim].weight;
    }
    else{
      lwtot += 1.0;
    }
  }

  MPI_Reduce(&lwtot,&gwtot,1,MPI_DOUBLE, MPI_SUM, 0, zz->Communicator);

  MPI_Allreduce(&max,&numParts,1,MPI_INT, MPI_MAX, zz->Communicator);
  numParts++;

  lpartWgt = (double *)ZOLTAN_CALLOC(numParts, sizeof(double));
  gpartWgt = (double *)ZOLTAN_MALLOC(numParts * sizeof(double));

  if (numParts && (!lpartWgt || !gpartWgt)){
    ZOLTAN_FREE(&lpartWgt);
    ZOLTAN_FREE(&gpartWgt);
    ZOLTAN_PRINT_ERROR(zz->Proc, "partition_stats", "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  for (i = 0, j=0, move=0.0; i < ndots; i++, j += wgtDim){
    if (wgtDim)
      lpartWgt[new_parts[i]] +=  (double)dots[j].weight;
    else
      lpartWgt[new_parts[i]] += 1.0;

    if (parts[i] != new_parts[i]){
      if (obj_sizes)
        move += (double)obj_sizes[i];
      else
        move += 1;
    }
  }

  MPI_Reduce(lpartWgt, gpartWgt, numParts, MPI_DOUBLE, MPI_SUM,
             0, zz->Communicator);
  MPI_Reduce(&move, &gmove, 1, MPI_DOUBLE, MPI_SUM, 0, zz->Communicator);

  ZOLTAN_FREE(&lpartWgt);

  if (proc == 0) {
    static int nRuns=0;
    static double balsum, balmax, balmin;
    static double movesum, movemax, movemin;
    char *countType;

    max_imbal = 0.0;

    if (gwtot) {
      for (i = 0; i < numParts; i++){
        if (work_fraction[i]) {
          ib= (gpartWgt[i]-work_fraction[i]*gwtot)/(work_fraction[i]*gwtot);
          if (ib>max_imbal)
            max_imbal = ib;
        }
      }
    }

    bal = 1.0 + max_imbal;

    if (nRuns){
      if (gmove > movemax) movemax = gmove;
      if (gmove < movemin) movemin = gmove;
      if (bal > balmax) balmax = bal;
      if (bal < balmin) balmin = bal;
      movesum += gmove;
      balsum += bal;
    }
    else{
      movemax = movemin = movesum = gmove;
      balmax = balmin = balsum = bal;
    }

    countType = "moveCnt";
    if (obj_sizes){
      countType = "moveVol";
    }

    nRuns++;
    printf(" STATS Runs %d  bal  CURRENT %f  MAX %f  MIN %f  AVG %f\n",
            nRuns, bal, balmax, balmin, balsum/nRuns);
    printf(" STATS Runs %d  %s CURRENT %f  MAX %f  MIN %f  AVG %f\n",
            nRuns, countType, gmove, movemax, movemin, movesum/nRuns);
  }

  ZOLTAN_FREE(&gpartWgt);
  MPI_Barrier(zz->Communicator);

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
