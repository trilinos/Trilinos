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

#include "hsfc.h"
#include "zz_const.h"
#include "zz_util_const.h"
#include <math.h>


/* For a detailed description of the following algorithm, please see the
   Developers Guide.  For instructions on use, please see the Users
   Guide.  This code assumes a 32 bit integer length!!! For a note on
   increasing the precision see the end of this file. */

static double next_query_2d (ZZ*, double *lquerybox, double *hquerybox, double);
static double next_query_3d (ZZ*, double *lquerybox, double *hquerybox, double);


/* returns list of processors and partitions intersecting the user's query box */
int Zoltan_HSFC_Box_Assign (
 ZZ *zz, double xlo, double ylo, double zlo,     /* low query box vertex */
         double xhi, double yhi, double zhi,     /* high query box vertex */
         int *procs, int *proc_count, int *parts, int *part_count)
   {
   double     xintl[3], xinth[3];         /* low and high bounded query box */
   int       *part_array = NULL;
   int       *proc_array = NULL;
   HSFC_Data *d;                    /* HFSC's portion of Zoltan data structure */
   int        n, i, loop;                 /* loop counters */
   int        tmp, dim = 0;
   int        first_proc, last_proc;
   double     fsfc, starting;
   double lo[3], hi[3];
   Partition *p;
   const double FUZZY = 1.0e-9; /* to build region about a query point or line */
                                /* Erik suggested that this value be tied to */
                                /* the actual precision available (dimension */
                                /* specific - 2^18 in 3d, 2^27 in 2d. */
   int        err = ZOLTAN_OK;
   int        *remap;
   char      *yo = "Zoltan_HSFC_Box_Assign";

   ZOLTAN_TRACE_ENTER (zz, yo);
   d = (HSFC_Data *) zz->LB.Data_Structure;           /* persistant HSFC data */
   if (d == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL,
          "No Decomposition Data available; use KEEP_CUTS parameter.");

   /* allocate memory to store results */
   *part_count = *proc_count = 0;
   part_array = (int *) ZOLTAN_MALLOC ((zz->LB.Num_Global_Parts + zz->Num_Proc)
    * sizeof (int));
   if (part_array == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Memory allocation error.");
   proc_array = part_array + zz->LB.Num_Global_Parts;
   memset (part_array, 0, (zz->LB.Num_Global_Parts + zz->Num_Proc)*sizeof (int));

   /* One dimensional case is trival, do it and exit */
   if (d->ndimension == 1) {
      remap = zz->LB.Remap;            /* Don't let Point_Assign remap IDs.*/
      zz->LB.Remap = NULL;             /* We will do this at "fini".       */
      Zoltan_HSFC_Point_Assign (zz, &xlo, NULL, &n);
      Zoltan_HSFC_Point_Assign (zz, &xhi, NULL, &loop);
      for (i = n; i <= loop; i++)
         part_array[i] = 1;
      zz->LB.Remap = remap;
      goto fini;
      }

   if (d->tran.Target_Dim > 0){   /* It must be 1 or 2 */
     /* 
      * Degenerate geometry:
      * Transform query box into coordinates that were used for partitioning,
      * and place an axis aligned bounding box around it.  This box in the new 
      * coordinates may encompass more "dots" than did the original box, but 
      * it won't miss any dots.
      */
     lo[0] = xlo; lo[1] = ylo; lo[2] = zlo;
     hi[0] = xhi; hi[1] = yhi; hi[2] = zhi;

     Zoltan_Transform_Box(lo, hi, d->tran.Transformation, d->tran.Permutation,
       d->ndimension, d->tran.Target_Dim);

     xlo = lo[0]; xhi = hi[0];
     ylo = 0.0; yhi = 0.0;
     zlo = 0.0; zhi = 0.0;

     if (d->tran.Target_Dim == 2){
       ylo = lo[1]; yhi = hi[1];
       dim = d->tran.Target_Dim;
     }
     else if (d->tran.Target_Dim == 1){
       /* 
        * Don't let Point_Assign transform coordinates (we already
        * did that) or remap partition numbers (we'll do that at "fini").
        */
       dim = d->ndimension;
       remap = zz->LB.Remap;
       d->tran.Target_Dim = 0;   /* don't transform coordinates */
       d->ndimension = 1;   /* partitions were calculated as 1D */
       zz->LB.Remap = NULL; /* don't remap partition numbers */
       Zoltan_HSFC_Point_Assign (zz, &xlo, NULL, &n);
       Zoltan_HSFC_Point_Assign (zz, &xhi, NULL, &loop);
       for (i = n; i <= loop; i++)  /* loop < n */
          part_array[i] = 1;
       d->tran.Target_Dim = 1;
       d->ndimension = dim;       
       zz->LB.Remap = remap;
       goto fini;
     }
   }
   else{
     dim = d->ndimension;
   }

   /* determine intersection of bounding box and query box */
   xintl[0] = (xlo - d->bbox_lo[0]) / d->bbox_extent[0];
   if      (xintl[0] < FUZZY)       xintl[0] = FUZZY;
   else if (xintl[0] > 1.0-FUZZY)   xintl[0] = 1.0-FUZZY;

   xinth[0] = (xhi - d->bbox_lo[0]) / d->bbox_extent[0];
   if      (xinth[0] < FUZZY)       xinth[0] = FUZZY;
   else if (xinth[0] > 1.0-FUZZY)   xinth[0] = 1.0-FUZZY;

   xintl[1] = (ylo - d->bbox_lo[1]) / d->bbox_extent[1];
   if      (xintl[1] < FUZZY)       xintl[1] = FUZZY;
   else if (xintl[1] > 1.0-FUZZY)   xintl[1] = 1.0-FUZZY;

   xinth[1] = (yhi - d->bbox_lo[1]) / d->bbox_extent[1];
   if      (xinth[1] < FUZZY)       xinth[1] = FUZZY;
   else if (xinth[1] > 1.0-FUZZY)   xinth[1] = 1.0-FUZZY;

   /* application programs need to add "dots" even if query box doesn't actually
      intersect unit cube.  FUZZY forces closest virtual overlap. */
   if (xinth[0] - xintl[0] < FUZZY)  {
       xintl[0] -= (FUZZY/2.0);
       xinth[0] += (FUZZY/2.0);
       }

   if (xinth[1] - xintl[1] < FUZZY)  {
       xintl[1] -= (FUZZY/2.0);
       xinth[1] += (FUZZY/2.0);
       }

   if (dim == 2)
      for (i = 0; i < zz->LB.Num_Global_Parts; i++)  {
         if (d->final_partition[i].l == d->final_partition[i].r)
            continue;            /* ignore empty partitions */

         /* find next query result in range [0.0, 1.0] */
         starting = (d->final_partition[i].l == 0.0)
          ? 0.0 : (d->final_partition[i].l + FUZZY);
         fsfc = next_query_2d (zz, xintl, xinth, starting);
         if (fsfc < 0.0)
            break;        /* done indication, no more partitions to be found */

         /* Find next nonempty partition entering or already in query region */
         p = (Partition *) bsearch (&fsfc, d->final_partition,
          zz->LB.Num_Global_Parts, sizeof (Partition), Zoltan_HSFC_compare);
         if (p == NULL)
            ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "programming error, can't happen");
         if (d->final_partition[p->index].l != d->final_partition[p->index].r)
            part_array [p->index] = 1;         /* mark non empty partitions */
         i = p->index;
         }

   if (dim == 3)  {
      /* complete the z axis information, as above */
      xintl[2] = (zlo - d->bbox_lo[2]) / d->bbox_extent[2];
      if      (xintl[2] < 0.0)   xintl[2] = FUZZY;
      else if (xintl[2] > 1.0)   xintl[2] = 1.0-FUZZY;

      xinth[2] = (zhi - d->bbox_lo[2]) / d->bbox_extent[2];
      if      (xinth[2] < 0.0)   xinth[2] = FUZZY;
      else if (xinth[2] > 1.0)   xinth[2] = 1.0-FUZZY;

      if (xinth[2] - xintl[2] < FUZZY)  {
          xintl[2] -= (FUZZY/2.0);
          xinth[2] += (FUZZY/2.0);
          }

      for (i = 0; i < zz->LB.Num_Global_Parts; i++)  {
         if (d->final_partition[i].l == d->final_partition[i].r)
            continue;                    /* ignore empty partition */

         /* find next query result in range [0.0. 1.0] */
         starting = (d->final_partition[i].l == 0.0)
          ? 0.0 : (d->final_partition[i].l + FUZZY);
         fsfc = next_query_3d (zz, xintl, xinth, starting);
         if (fsfc < 0.0)
            break;                       /* done, no more to be found */

         /* Find partition containing point and return its number */
         p = (Partition *) bsearch (&fsfc, d->final_partition,
          zz->LB.Num_Global_Parts, sizeof (Partition), Zoltan_HSFC_compare);
         if (p == NULL)
            ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "programming error, can't happen");
         if (d->final_partition[p->index].l != d->final_partition[p->index].r)
            part_array [p->index] = 1;    /* ignore empty partitions */
         i = p->index;
         }
      }


   /* move results to user supplied arrays & hope they are big enough */
fini:
   *part_count = 0;
   if (parts)  {
      *part_count = 0;
      for (i = 0; i < zz->LB.Num_Global_Parts; i++)
         if (part_array[i] > 0) {
            if (zz->LB.Remap)
               parts[(*part_count)++] = zz->LB.Remap[i];
            else
               parts[(*part_count)++] = i;
            }
      }

   if (procs)  {
      for (i = 0; i < zz->LB.Num_Global_Parts; i++)
         if (part_array[i] > 0) {
            tmp = (zz->LB.Remap ? zz->LB.Remap[i] : i);
            first_proc = Zoltan_LB_Part_To_Proc(zz, tmp, NULL);
            proc_array[first_proc] = 1;
            if (!zz->LB.Single_Proc_Per_Part) {
               /* Part may be spread across multiple procs. Include them all. */
               int j;
               if (tmp < zz->LB.Num_Global_Parts - 1)
                  last_proc = Zoltan_LB_Part_To_Proc (zz, tmp + 1, NULL);
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

End:
   ZOLTAN_FREE (&part_array);
   ZOLTAN_TRACE_EXIT (zz, yo);
   return err;
   }



/* returns the first Hilbert coordinate in [s,1] to enter the user's query box */
static double next_query_2d (ZZ *zz, double *lquerybox, double *hquerybox,
 double s)
   {
   int state, newstate=0, savelevel, prune, backtrack;
   unsigned int temp, x, y;
   unsigned int savestate=0, savequad=0, savenpt=0;  /* saved info for backtracking */
   unsigned int qlox, qloy, qhix, qhiy;        /* query box bounds */
   unsigned int nptx, npty, keyx, keyy, newnpt=0;
   unsigned int start[2], startbits=0;
   unsigned int intersect_hi, intersect_lo;
   double t;
   unsigned i, level, quadrant;                    /* loop counters */
   static const unsigned int *dk[]
    = {data2d,  data2d  +4, data2d  +8, data2d  +12};
   static const unsigned int *st[]
    = {state2d, state2d +4, state2d +8, state2d +12};
   static const unsigned int MAXLEVEL = 28;  /* only 56 significant bits, 28 per dimension */

   /* convert floating normalized, intersected query box corners to integers */
   qlox = (unsigned int) (lquerybox[0] * (double) IMAX);
   qloy = (unsigned int) (lquerybox[1] * (double) IMAX);
   qhix = (unsigned int) (hquerybox[0] * (double) IMAX);
   qhiy = (unsigned int) (hquerybox[1] * (double) IMAX);

   /* convert floating minimum query hilbert coordinate to integer */
   start[1] = (unsigned int) (modf (s * (double) IMAX, &t) * (double) IMAX);
   start[0] = (unsigned int) t;
   
   /* initializations before starting main loop */
   state = 0;
   prune = 1;
   backtrack = 1;
   savelevel = -1;
   nptx = npty = 0;
   keyx = keyy = 0;

   /* now loop over each level of hilbert curve refinement */
   for (level = 0; level < MAXLEVEL; level++)  {
      if (prune)  {
         /* get next 2 bits of start for pruning from its shift register */
         startbits =  start[0] >> 30;                 /* extract top two bits */
         start[0]  = (start[0] << 2) | (start[1] >> 30);  /* shift hi word */
         start[1]  =  start[1] << 2;                      /* shift lo word */
         }

      /* compute current search space intersected with the query region */
      temp = IMAX >> level;        /* top of range ends in all 1 bits */
      x = ((nptx | temp) > qhix) ? qhix : (nptx | temp);
      y = ((npty | temp) > qhiy) ? qhiy : (npty | temp);
      intersect_hi = (((x >> (30-level)) & 2) | ((y >> (31-level)) & 1));

      temp = ~temp;                /* bottom of range ends in all 0 bits */
      x = ((nptx & temp) < qlox) ? qlox : (nptx & temp);
      y = ((npty & temp) < qloy) ? qloy : (npty & temp);
      intersect_lo = (((x >> (30-level)) & 2) | ((y >> (31-level)) & 1));

      /* loop over subquadrants in hilbert curve order to find lowest numbered
         quad that intersects with query region and is larger than start key */
      temp = (prune) ? startbits : 0;
      for (quadrant = temp; quadrant < 4; quadrant++)  {
         newnpt = *(dk[state] + quadrant);         /* get new 2 bit npt value */
         if (!((newnpt ^ intersect_hi) & (newnpt ^ intersect_lo)))
            break;   /* successfully found intersecting quadrant at this level */
         }
      if (quadrant < 4) {
         newstate = *(st[state] + quadrant);
         if  (prune && (quadrant > startbits))
            prune = 0;                            /* no more pruning required */
         }

      /* determine backtracking point, in case backtracking is needed */
      if (backtrack)
         for (i = quadrant+1; i < 4; i++)  {
            /* intersect test - subquad intersecting with search/query region */
            temp = *(dk[state] + i);               /* get new npt */
            if (!((temp ^ intersect_hi) & (temp ^ intersect_lo)))  {
               savelevel = level;
               savestate = *(st[state] + i);
               savenpt   = temp;
               savequad  = i;
               break;
               }
            }

      /* test if this round was successful (quadrant < 4), backtrack if not */
      if (quadrant == 4)  {
         if (savelevel == -1)
            return -1.0;             /* report that no solution is available */

         /* restore all conditions to backtrack point */
         level    = savelevel;
         newstate = savestate;
         newnpt   = savenpt;
         quadrant = savequad;
         prune     = 0;          /* no longer need to prune previous branches */
         backtrack = 0;          /* only 1 backtrack ever required, now done */

         /* discard results below backtrack level */
         keyx &= ~(IMAX >> savelevel);
         keyy &= ~(IMAX >> savelevel);
         
         nptx &= ~(IMAX >> savelevel);
         nptx &= ~(IMAX >> savelevel);
         }

      /* append current derived key and npoint to cumulative total */
      keyx |= ((quadrant & 2) << (30-level));
      keyy |= ((quadrant & 1) << (31-level));

      nptx |= ((newnpt   & 2) << (30-level));
      npty |= ((newnpt   & 1) << (31-level));

      state  = newstate;
      }

   /* now convert final derived key to floating point representation and exit */
   start[0] = start[1] = 0;
   for (i = 0; i < MAXLEVEL; i++)  {
      start[0] = (start[0] << 2) | (start[1] >> 30);
      start[1] = (start[1] << 2) | ((keyx >> (30-i)) & 2)
                                 | ((keyy >> (31-i)) & 1);
      }
   return ldexp ((double) start[0], -24) + ldexp((double) start[1], -56);
   }



/* returns the first Hilbert coordinate in [s,1] to enter the user's query box */
static double next_query_3d (ZZ *zz, double *lquerybox, double *hquerybox,
 double s)
   {
   int state, newstate = 0, savelevel, prune, backtrack;
   unsigned int temp, x, y, z;
   unsigned int savestate = 0, savequad = 0, savenpt = 0;
   unsigned int qlox, qloy, qloz, qhix, qhiy, qhiz;
   unsigned int nptx, npty, nptz, keyx, keyy, keyz, newnpt = 0;
   unsigned int start[2], startbits = 0;
   unsigned int intersect_hi, intersect_lo;
   double t;
   unsigned i, quadrant, level;
   static const unsigned int *dk[] =
      {data3d,      data3d +8,   data3d +16,  data3d +24,
       data3d +32,  data3d +40,  data3d +48,  data3d +56,
       data3d +64,  data3d +72,  data3d +80,  data3d +88,
       data3d +96,  data3d +104, data3d +112, data3d +120,
       data3d +128, data3d +136, data3d +144, data3d +152,
       data3d +160, data3d +168, data3d +176, data3d +184};

   static const unsigned int *st[] =
      {state3d,      state3d +8,   state3d +16,  state3d +24,
       state3d +32,  state3d +40,  state3d +48,  state3d +56,
       state3d +64,  state3d +72,  state3d +80,  state3d +88,
       state3d +96,  state3d +104, state3d +112, state3d +120,
       state3d +128, state3d +136, state3d +144, state3d +152,
       state3d +160, state3d +168, state3d +176, state3d +184};

   static const unsigned int MAXLEVEL = 18;  /* only 56 significant bits, 18 per dimension */

   /* convert floating query box corners to integers */
   qlox = (unsigned int) (lquerybox[0] * (double) IMAX);
   qloy = (unsigned int) (lquerybox[1] * (double) IMAX);
   qloz = (unsigned int) (lquerybox[2] * (double) IMAX);

   qhix = (unsigned int) (hquerybox[0] * (double) IMAX);
   qhiy = (unsigned int) (hquerybox[1] * (double) IMAX);
   qhiz = (unsigned int) (hquerybox[2] * (double) IMAX);

   /* convert floating minimum query hilbert coordinate to integer */
   start[1] = (unsigned int) (modf (s * (double) IMAX, &t) * (double) IMAX);
   start[0] = (unsigned int) t;

   /* initializations before main loop */
   state = 0;
   prune = 1;
   backtrack = 1;
   savelevel = -1;
   nptx = npty = nptz = 0;
   keyx = keyy = keyz = 0;

   /* now loop over each level of hilbert curve refinement */
   for (level = 0; level < MAXLEVEL; level++)  {
      if (prune)  {
         /* get next 3 bits of start for pruning, treat as shift register */
         startbits =  start[0] >> 29;                 /* extract top two bits */
         start[0]  = (start[0] << 3) | (start[1] >> 29);  /* shift hi word */
         start[1]  =  start[1] << 3;                      /* shift lo word */
         }

      /* compute current search space intersected with the query region */
      temp =  IMAX >> level;        /* top of range ends in all 1 bits */
      x = ((nptx | temp) > qhix) ? qhix : (nptx | temp);
      y = ((npty | temp) > qhiy) ? qhiy : (npty | temp);
      z = ((nptz | temp) > qhiz) ? qhiz : (nptz | temp);
      intersect_hi = ((x >> (29-level)) & 4) 
                   | ((y >> (30-level)) & 2)
                   | ((z >> (31-level)) & 1);

      temp = ~temp;             /* bottom of range ends in all 0 bits */
      x = ((nptx & temp) < qlox) ? qlox : (nptx & temp);
      y = ((npty & temp) < qloy) ? qloy : (npty & temp);
      z = ((nptz & temp) < qloz) ? qloz : (nptz & temp);
      intersect_lo = ((x >> (29-level)) & 4) 
                   | ((y >> (30-level)) & 2)
                   | ((z >> (31-level)) & 1);

      /* loop over subquadrants in hilbert curve order to find lowest numbered
         quad that intersects with query region and is larger than start key */
      temp = (prune) ? startbits : 0;
      for (quadrant = temp; quadrant < 8; quadrant++)  {
         newnpt = *(dk[state] + quadrant);         /* get new 3 bit npt value */
         if (!((newnpt ^ intersect_hi) & (newnpt ^ intersect_lo)))
            break;   /* successfully found intersecting quadrant at this level */
         }
      if (quadrant < 8)  {
         newstate = *(st[state] + quadrant);
         if  (prune && (quadrant > startbits))
            prune = 0;                            /* no more pruning required */
         }

      /* determine backtracking point, in case backtracking is needed */
      if (backtrack)
         for (i = quadrant+1; i < 8; i++)  {
            temp = *(dk[state] + i);         /* get new npt */

            /* intersect test - subquad intersecting with search/query region */
            if (!((temp ^ intersect_hi) & (temp ^ intersect_lo)))  {
               savelevel = level;
               savestate = *(st[state] + i);
               savenpt   = temp;
               savequad  = i;
               break;
               }
            }

      /* test if this round was successful (quadrant < 8), backtrack if not */
      if (quadrant == 8)  {
         if (savelevel == -1)
            return -1.0;             /* no solution available, report failure */

         /* restore all conditions to backtrack point */
         level    = savelevel;
         newstate = savestate;
         newnpt   = savenpt;
         quadrant = savequad;
         prune     = 0;          /* no longer need to prune previous branches */
         backtrack = 0;          /* only 1 backtrack ever required, now done */

         /* discard results below backtrack level */
         keyx &= ~(IMAX >> savelevel);
         keyy &= ~(IMAX >> savelevel);
         keyz &= ~(IMAX >> savelevel);

         nptx &= ~(IMAX >> savelevel);
         npty &= ~(IMAX >> savelevel);
         nptz &= ~(IMAX >> savelevel);
         }

      /* append current derived key and npoint to cumulative total */
      state  = newstate;

      keyx |= ((quadrant & 4) << (29-level));
      keyy |= ((quadrant & 2) << (30-level));
      keyz |= ((quadrant & 1) << (31-level));

      nptx |= ((newnpt & 4)   << (29-level));
      npty |= ((newnpt & 2)   << (30-level));
      nptz |= ((newnpt & 1)   << (31-level));
      }

   /* now convert final derived key to floating point representation and exit */
   start[0] = start[1] = 0;
   for (i = 0; i < MAXLEVEL; i++)  {
      start[0] = (start[0] << 3) | (start[1] >> 29);
      start[1] = (start[1] << 3) | ((keyx >> (29-i)) & 4)
                                 | ((keyy >> (30-i)) & 2) 
                                 | ((keyz >> (31-i)) & 1);
      }
   return ldexp ((double) start[0], -22) + ldexp((double) start[1], -54);
   }


/* Maintenance Note:  Per the review 04.15.03, the following section addresses
increasing the precision of this box assign.

The next_query_xx() routines can be extended to arbitrary precision. Currently
they reflect the precision limitation of using a double for the starting
Hilbert coordinate and the returned value:
   static double next_query_3d (ZZ *zz, double *lquerybox, double *hquerybox,
    double s)
Change this to:
   static int* next_query_3d (ZZ *zz, double *lquerybox, double *hquerybox,
    int *start).
Represent the Hilbert coordinate as an array of ints (2 for 2d, 3 for 3d) rather
than the double representation currently used.
Change the declaration:
   static const int MAXLEVEL = 18; (3d or 28 in 2d)
to
   static const int MAXLEVEL = 32;
Change the return statement:
   return ldexp ((double) start[0], -22) + ldexp((double) start[1], -54);
to
   return start;

Please see the related remarks in the file hsfc_hilbert.c.

The last change is to replace the criteria for the binary lookup of a hilbert
coordinate to find its partition to use the int arrays rather than double.
This requires a trivial change to hsfc.c to use int arrays for the hilbert
coordinate as well.
*/


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
