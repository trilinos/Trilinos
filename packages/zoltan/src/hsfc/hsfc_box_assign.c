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
#include <math.h>


/* For a detailed description of the following algorithm, please see the
   Developers Guide.  For instructions on use, please see the Users
   Guide.  */

static double next_query_2d (ZZ*, double *lquerybox, double *hquerybox, double);
static double next_query_3d (ZZ*, double *lquerybox, double *hquerybox, double);


/* returns list of processors and partitions falling within user's query box */
int Zoltan_HSFC_Box_Assign (
 ZZ *zz, double xlo, double ylo, double zlo,
         double xhi, double yhi, double zhi,
         int *procs, int *proc_count, int *parts, int *part_count)
   {
   double     xintl[3], xinth[3];         /* low and high bounded query box */
   double     start_time, end_time;       /* timing information */
   int       *part_array = NULL;
   int       *proc_array = NULL;
   HSFC_Data *d;
   int        n, i, loop;                 /* loop counters */
   int        first_proc, last_proc;
   double     fsfc, starting;
   Partition *p;
   const double FUZZY = 1.0e-9;
   int        err = ZOLTAN_OK;
   char *yo = "Zoltan_HSFC_Box_Assign";

   ZOLTAN_TRACE_ENTER (zz, yo);
   d = (HSFC_Data *) zz->LB.Data_Structure;
   if (d == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL,
          "No Decomposition Data available; use KEEP_CUTS parameter.");
   start_time = Zoltan_Time (zz->Timer);

   /* allocate memory to store results */
   *part_count = *proc_count = 0;
   part_array = (int *) ZOLTAN_MALLOC ((zz->LB.Num_Global_Parts + zz->Num_Proc)
    * sizeof (int));
   if (part_array == NULL)
      ZOLTAN_HSFC_ERROR (ZOLTAN_MEMERR, "Memory error.");
   proc_array = part_array + zz->LB.Num_Global_Parts;
   memset (part_array, 0, (zz->LB.Num_Global_Parts + zz->Num_Proc)*sizeof (int));

   /* One dimensional case is trival, do it and exit */
   if (d->ndimension == 1) {
      Zoltan_HSFC_Point_Assign (zz, &xlo, NULL, &n);
      Zoltan_HSFC_Point_Assign (zz, &xhi, NULL, &loop);
      for (i = n; i <= loop; i++)
         part_array[i] = 1;
      goto fini;
      }

   /* determine intersection of bounding box and dropped box */
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
   /* intersect unit cube.  FUZZY forces closest virtual overlap. */
   if (xinth[0] - xintl[0] < FUZZY)  {
       xintl[0] -= (FUZZY/2.0);
       xinth[0] += (FUZZY/2.0);
       }

   if (xinth[1] - xintl[1] < FUZZY)  {
       xintl[1] -= (FUZZY/2.0);
       xinth[1] += (FUZZY/2.0);
       }

   if (d->ndimension == 2)
      for (i = 0; i < zz->LB.Num_Global_Parts; i++)  {
         if (d->final_partition[i].l == d->final_partition[i].r)
            continue;            /* ignore empty partitions */

         /* find next query result in range [0.0, 1.0] */
         starting = (d->final_partition[i].l == 0.0)
          ? 0.0 : (d->final_partition[i].l + FUZZY);
         fsfc = next_query_2d (zz, xintl, xinth, starting);
         if (fsfc < 0.0)
            break;               /* done, no more partitions to be found */

         /* Find next nonempty partition entering or already in query region */
         p = (Partition *) bsearch (&fsfc, d->final_partition,
          zz->LB.Num_Global_Parts, sizeof (Partition), Zoltan_HSFC_compare);
         if (p == NULL)
            ZOLTAN_HSFC_ERROR (ZOLTAN_FATAL, "programming error, can't happen");
         if (d->final_partition[p->index].l != d->final_partition[p->index].r)
            part_array [p->index] = 1;         /* mark non empty partitions */
         i = p->index;
         }

   if (d->ndimension == 3)  {
      xintl[2] = (zlo - d->bbox_lo[2]) / d->bbox_extent[2];
      if      (xintl[2] < 0.0)   xintl[2] = 0.0;
      else if (xintl[2] > 1.0)   xintl[2] = 1.0;

      xinth[2] = (zhi - d->bbox_lo[2]) / d->bbox_extent[2];
      if      (xinth[2] < 0.0)   xinth[2] = 0.0;
      else if (xinth[2] > 1.0)   xinth[2] = 1.0;

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
         if (part_array[i] > 0)
            parts[(*part_count)++] = i;
      }

   if (procs)  {
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

   end_time = Zoltan_Time (zz->Timer);
   if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME && zz->Proc == 0)
      printf ("HSFC Box Assign time is %.6f seconds\n", end_time - start_time);

   return err;
   }



/* finds the next partition to enter the query space */
static double next_query_2d (ZZ *zz, double *lquerybox, double *hquerybox,
 double s)
   {
   int state, newstate, savelevel, prune, backtrack;
   unsigned int temp, x, y;
   unsigned int savestate, savequad, savenpt;  /* saved info for backtracking */
   unsigned int qlox, qloy, qhix, qhiy;        /* query box bounds */
   unsigned int nptx, npty, keyx, keyy, newnpt;
   unsigned int start[2], startbits;
   unsigned int intersect_hi, intersect_lo;
   double t;
   int level, quadrant, i;                    /* loop counters */
   static const unsigned *dk[] = {idata2d,  idata2d  +4, idata2d  +8, idata2d  +12};
   static const unsigned *st[] = {istate2d, istate2d +4, istate2d +8, istate2d +12};
   static const MAXLEVEL = 28;  /* only 56 significant bits, 28 per dimension */

   /* convert floating normalized, intersected query box corners to integers */
   qlox = (unsigned int) (lquerybox[0] * (double) IMAX);
   qloy = (unsigned int) (lquerybox[1] * (double) IMAX);
   qhix = (unsigned int) (hquerybox[0] * (double) IMAX);
   qhiy = (unsigned int) (hquerybox[1] * (double) IMAX);

   /* convert floating minimum query hilbert coordinate to integer */
   start[1] = (unsigned int) (modf (s * (double) IMAX, &t) * (double) IMAX);
   start[0] = (unsigned int) t;

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

      /* loop over subquadrants in hilbert curve order to find lowest that
         intersects with query region and is larger than start key */
      for (quadrant = 0; quadrant < 4; quadrant++)  {
         if  (prune && quadrant < startbits) continue;  /* subtree < startbits */
         if  (prune && quadrant > startbits) prune = 0; /* stop pruning */

         /* intersection test - subquad intersecting with search/query region */
         newnpt = *(dk[state] + quadrant);         /* get new 2 bit npt value */
         if (!((newnpt ^ intersect_hi) & (newnpt ^ intersect_lo)))
            break;   /* successfully found intersecting quadrant at this level */
         }
      newstate = *(st[state] + quadrant);

      /* determine backtracking point, in case backtracking is needed */
      if (backtrack == 1)
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

         /* discard results below backtrack level */
         nptx &= ~(IMAX >> savelevel);
         npty &= ~(IMAX >> savelevel);
         keyx &= ~(IMAX >> savelevel);
         keyy &= ~(IMAX >> savelevel);

         prune     = 0;  /* no longer need to prune previous branches */
         backtrack = 0;  /* only 1 backtrack ever required, now done */
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

/******************************************************************************/



static double next_query_3d (ZZ *zz, double *lquerybox, double *hquerybox,
 double s)
   {
   int state, newstate, savelevel, prune, backtrack;
   unsigned int temp, x, y, z;
   unsigned int savestate, savequad, savenpt;
   unsigned int qlox, qloy, qloz, qhix, qhiy, qhiz;
   unsigned int nptx, npty, nptz, keyx, keyy, keyz, newnpt;
   unsigned int start[2], startbits;
   unsigned int intersect_hi, intersect_lo;
   double t;
   int level, quadrant, i;
   static const unsigned int *dk[] =
      {idata3d,      idata3d +8,   idata3d +16,  idata3d +24,
       idata3d +32,  idata3d +40,  idata3d +48,  idata3d +56,
       idata3d +64,  idata3d +72,  idata3d +80,  idata3d +88,
       idata3d +96,  idata3d +104, idata3d +112, idata3d +120,
       idata3d +128, idata3d +136, idata3d +144, idata3d +152,
       idata3d +160, idata3d +168, idata3d +176, idata3d +184};

   static const unsigned int *st[] =
      {istate3d,      istate3d +8,   istate3d +16,  istate3d +24,
       istate3d +32,  istate3d +40,  istate3d +48,  istate3d +56,
       istate3d +64,  istate3d +72,  istate3d +80,  istate3d +88,
       istate3d +96,  istate3d +104, istate3d +112, istate3d +120,
       istate3d +128, istate3d +136, istate3d +144, istate3d +152,
       istate3d +160, istate3d +168, istate3d +176, istate3d +184};

   static const MAXLEVEL = 18;  /* only 56 significant bits, 18 per dimension */

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

      /* loop over subquadrants in hilbert curve order to find lowest that
         intersects with query region and is larger than start key */
      for (quadrant = 0; quadrant < 8; quadrant++)  {
         if  (prune && quadrant < startbits)   continue;   /* subtree < startbits */
         if  (prune && quadrant > startbits)   prune = 0;  /* pruning not required */

         /* intersection test - subquad intersecting with search/query region */
         newnpt = *(dk[state] + quadrant);         /* get new 3 bit npt value */
         if (!((newnpt ^ intersect_hi) & (newnpt ^ intersect_lo)))
            break;   /* successfully found quadrant at this level */
         }
      newstate = *(st[state] + quadrant);

      /* determine backtracking point, in case backtracking is needed */
      if (backtrack == 1)
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

         /* discard results below backtrack level */
         nptx &= ~(IMAX >> savelevel);
         npty &= ~(IMAX >> savelevel);
         nptz &= ~(IMAX >> savelevel);
         
         keyx &= ~(IMAX >> savelevel);
         keyy &= ~(IMAX >> savelevel);
         keyz &= ~(IMAX >> savelevel);

         prune     = 0;     /* no longer need to prune previous branches */
         backtrack = 0;     /* only 1 backtrack ever required, now done */
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

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
