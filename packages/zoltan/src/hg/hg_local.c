/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
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

#include "hypergraph.h"

typedef struct
   {
   double weight;
   float gain;
   int   vertex;
   int   source;
   int   destination;
   } Vdata;

static int isabove (int origin, int test, int p);
static int isbelow (int origin, int test, int p);
static int comparison (const void*, const void*);
static int comparison2 (const void*, const void*);



static ZOLTAN_HG_LOCAL_REF_FN local_no;
static ZOLTAN_HG_LOCAL_REF_FN local_fm2;
static ZOLTAN_HG_LOCAL_REF_FN local_fmkway;
static ZOLTAN_HG_LOCAL_REF_FN local_grkway;

/****************************************************************************/



ZOLTAN_HG_LOCAL_REF_FN *Zoltan_HG_Set_Local_Ref_Fn(char *str)
{
  if      (!strcasecmp(str, "fm2"))    return local_fm2;
  else if (!strcasecmp(str, "fmkway")) return local_fmkway;
  else if (!strcasecmp(str, "grkway")) return local_grkway;
  else if (!strcasecmp(str, "no"))     return local_no;
  else                                 return NULL;
}

/****************************************************************************/



int Zoltan_HG_Local(ZZ *zz, HGraph *hg, int p, Partition part, HGPartParams *hgp)
{
  return hgp->local_ref(zz, hg, p, part, hgp, hgp->bal_tol);
}

/****************************************************************************/


static int local_no (
  ZZ *zz,     /* Zoltan data structure */
  HGraph *hg,
  int p,
  Partition part,
  HGPartParams *hgp,
  float bal_tol
)
{
  return ZOLTAN_OK;
}

/****************************************************************************/



/*
static int gain_check (HGraph *hg, float *gain, int *part, int **cut)
{
float g;
int vertex, j, edge;

  for (vertex = 0; vertex < hg->nVtx; vertex++) {
     g = 0.0;
     for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
        edge = hg->vedge[j];
        if (cut[part[vertex]][edge] == 1)
           g += (hg->ewgt ? (hg->ewgt[edge]) : 1.0);
        else if (cut[1-part[vertex]][edge] == 0)
           g -= (hg->ewgt ? (hg->ewgt[edge]) : 1.0);
        }
     if (g != gain[vertex]) {
        printf("Wrong gain %f %f\n", g, gain[vertex]);
        return ZOLTAN_FATAL;
        }
     }
  return ZOLTAN_OK;
}
*/

/****************************************************************************/



int Zoltan_HG_move_vertex (HGraph *hg, int vertex, int sour, int dest,
 int *part, int **cut, double *gain, HEAP *heap)
{
int i, j, edge, v;

  gain[vertex] = 0.0;
  part[vertex] = dest;

  for (i = hg->vindex[vertex]; i < hg->vindex[vertex+1]; i++) {
     edge = hg->vedge[i];
     if (cut[sour][edge] == 1) {
        for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
           v = hg->hvertex[j];
           gain[v] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
           if (heap)
              Zoltan_HG_heap_change_value(&heap[part[v]], v, gain[v]);
           }
        }
     else if (cut[sour][edge] == 2) {
        for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
           v = hg->hvertex[j];
           if (part[v] == sour) {
              gain[v] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
              if (heap)
                 Zoltan_HG_heap_change_value(&heap[part[v]], v, gain[v]);
              break;
              }
           }
        }

     if (cut[dest][edge] == 0) {
        for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
           v = hg->hvertex[j];
           gain[v] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
           if (heap)
              Zoltan_HG_heap_change_value(&heap[part[v]], v, gain[v]);
           }
        }
     else if (cut[dest][edge] == 1) {
        for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
           v = hg->hvertex[j];
           if (v != vertex && part[v] == dest) {
              gain[v] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
              if (heap)
                 Zoltan_HG_heap_change_value(&heap[part[v]], v, gain[v]);
              break;
              }
           }
        }
     cut[sour][edge]--;
     cut[dest][edge]++;
     }
  return ZOLTAN_OK;
}



/****************************************************************************/

static int local_fm2 (
  ZZ *zz,
  HGraph *hg,
  int p,
  Partition part,
  HGPartParams *hgp,
  float bal_tol
)
{
int    i, j, vertex, edge, *cut[2], *locked = 0, *locked_list = 0, round = 0;
double total_weight, part_weight[2], max_weight[2];
double cutsize, best_cutsize, *gain = 0;
float  junk;
HEAP   heap[2];
int    steplimit;
char   *yo="local_fm2";

double ratio = hg->ratio, error, best_error;
int    best_imbalance, imbalance;

  if (p != 2) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "p!=2 not allowed for local_fm2.");
     return ZOLTAN_FATAL;
     }

  if (hg->nEdge == 0)
     return ZOLTAN_OK;

  /* Calculate the weights in each partition and total, then maxima */
  part_weight[0] = 0.0;
  part_weight[1] = 0.0;
  if (hg->vwgt)  {
     for (i = 0; i < hg->nVtx; i++)
        part_weight[part[i]] += hg->vwgt[i];
     total_weight = part_weight[0] + part_weight[1];
     }
  else  {
     total_weight = (double)(hg->nVtx);
     for (i = 0; i < hg->nVtx; i++)
        part_weight[part[i]] += 1.0;
     }
  max_weight[0] = total_weight * bal_tol *      ratio;
  max_weight[1] = total_weight * bal_tol * (1 - ratio);

  if (!(cut[0]      = (int*)   ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
   || !(locked      = (int*)   ZOLTAN_CALLOC    (hg->nVtx,  sizeof(int)))
   || !(locked_list = (int*)   ZOLTAN_CALLOC    (hg->nVtx,  sizeof(int)))
   || !(gain        = (double*)ZOLTAN_CALLOC    (hg->nVtx,  sizeof(double))) ) {
         Zoltan_Multifree(__FILE__,__LINE__, 4, &cut[0], &locked, &locked_list,
          &gain);
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
         }
  cut[1] = &(cut[0][hg->nEdge]);

  /* Initial calculation of the cut distribution and gain values */
  for (i = 0; i < hg->nEdge; i++)
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
        (cut[part[hg->hvertex[j]]][i])++;
  for (i = 0; i < hg->nVtx; i++)
     for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)
        {
        edge = hg->vedge[j];
        if (cut[part[i]][edge] == 1)
           gain[i] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
        else if (cut[1-part[i]][edge] == 0)
           gain[i] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
        }

  /* Initialize the heaps and fill them with the gain values */
  Zoltan_HG_heap_init(zz, &heap[0], hg->nVtx);
  Zoltan_HG_heap_init(zz, &heap[1], hg->nVtx);  
  for (i = 0; i < hg->nVtx; i++)
     Zoltan_HG_heap_input(&heap[part[i]], i, gain[i]);
  Zoltan_HG_heap_make(&heap[0]);
  Zoltan_HG_heap_make(&heap[1]);

  /* Initialize given partition as best partition */
  best_cutsize = cutsize = Zoltan_HG_hcut_size_total(hg, part);
  best_error = MAX (part_weight[0]-max_weight[0], part_weight[1]-max_weight[1]);
  best_imbalance = (part_weight[0]>max_weight[0])||(part_weight[1]>max_weight[1]);
  do {
    int step = 0, no_better_steps = 0, number_locked = 0, best_locked = 0;
    int sour, dest;
    double akt_cutsize=best_cutsize;

/*
int *rth = (int *) &locked_list[hg->nVtx-1];
*/

    round++;
    cutsize = best_cutsize;
    if (hgp->output_level > HG_DEBUG_LIST)
      printf("ROUND %d:\nSTEP VERTEX  PARTS MAX_WGT CHANGE CUTSIZE\n",round);

    steplimit = (hgp->noimprove_limit > 1.1) ? (int) hgp->noimprove_limit : (int) (hg->nVtx * hgp->noimprove_limit) ;
    /* steplimit = hg->nVtx/4;  Robsys previous choice */

    while (step < hg->nVtx && no_better_steps < steplimit) {
        step++;
        no_better_steps++;

        if (Zoltan_HG_heap_empty(&heap[0]))
           sour = 1;
        else if (Zoltan_HG_heap_empty(&heap[1]))
           sour = 0;
        else if (part_weight[0] > max_weight[0])
           sour = 0;
        else if (part_weight[1] > max_weight[1])
           sour = 1;
/*
else if (Zoltan_HG_heap_max_value(&heap[0]) == Zoltan_HG_heap_max_value(&heap[1])
 && hgp->tiestrategy == 2 && part_weight[0] > part_weight[1])
    sour = 0;
else if (Zoltan_HG_heap_max_value(&heap[0]) == Zoltan_HG_heap_max_value(&heap[1])
 && hgp->tiestrategy == 2 && part_weight[1] > part_weight[0])
   sour = 1;

else if (Zoltan_HG_heap_max_value(&heap[0]) == Zoltan_HG_heap_max_value(&heap[1])
 && hgp->tiestrategy == 1)
   sour = step % 2;
*/
        else if (Zoltan_HG_heap_max_value(&heap[0])
              >  Zoltan_HG_heap_max_value(&heap[1]))
           sour = 0;
        else
           sour = 1;
        dest = 1-sour;
        vertex = Zoltan_HG_heap_extract_max(&heap[sour], &junk);

/*
if ((hg->vwgt[vertex] + part_weight[dest] > max_weight[dest]) &&
(hgp->fmswitch == -2 || (hgp->fmswitch >= 0 && hg->info % 2 == hgp->fmswitch)))
   {
   *rth-- = vertex;
   continue;
   }
*/
        locked[vertex] = part[vertex] + 1;
        locked_list[number_locked++] = vertex;
        akt_cutsize -= gain[vertex];

        Zoltan_HG_move_vertex (hg, vertex, sour, dest, part, cut, gain, heap);
        part_weight[sour] -= (hg->vwgt ? hg->vwgt[vertex] : 1.0);
        part_weight[dest] += (hg->vwgt ? hg->vwgt[vertex] : 1.0);

        error = MAX (part_weight[0]-max_weight[0],part_weight[1]-max_weight[1]);
        imbalance = (part_weight[0]>max_weight[0])||(part_weight[1]>max_weight[1]);

#ifdef RTHRTH
if (step == 1)
  best_error = error;
#endif
        if ( ( best_imbalance && (error < best_error))
          || (!imbalance && (akt_cutsize < best_cutsize)))  {
            best_error   = error;
            best_imbalance = imbalance;
            best_locked  = number_locked;
            best_cutsize = akt_cutsize;
            no_better_steps = 0;
            }
        if (hgp->output_level > HG_DEBUG_LIST+1)
           printf ("%4d %6d %2d->%2d %7.2f %f %f\n", step, vertex, sour, dest,
            error, akt_cutsize - cutsize, akt_cutsize);
        }

     while (number_locked != best_locked) {
        vertex = locked_list[--number_locked];
        sour = part[vertex];
        dest = locked[vertex] - 1;

        Zoltan_HG_move_vertex (hg, vertex, sour, dest, part, cut, gain, heap);

        part_weight[sour] -= (hg->vwgt ? hg->vwgt[vertex] : 1.0);
        part_weight[dest] += (hg->vwgt ? hg->vwgt[vertex] : 1.0);
        Zoltan_HG_heap_input(&heap[dest], vertex, gain[vertex]);
        locked[vertex] = 0;
        }

     while (number_locked) {
        vertex = locked_list[--number_locked];
        locked[vertex] = 0;
        Zoltan_HG_heap_input(&heap[part[vertex]], vertex, gain[vertex]);
        }

/*
while (rth < &locked_list [hg->nVtx-1])
   {
   rth++;
   vertex = *rth;
   *rth = 0;
   Zoltan_HG_heap_input(&heap[part[vertex]], vertex, gain[vertex]);
   }
*/

     Zoltan_HG_heap_make(&(heap[0]));
     Zoltan_HG_heap_make(&(heap[1]));
     } while ((best_cutsize < cutsize) ||  (round < hgp->nlevelrepeat));

  /* gain_check (hg, gain, part, cut); */
  Zoltan_Multifree(__FILE__,__LINE__, 4, &cut[0], &locked, &locked_list, &gain);
  Zoltan_HG_heap_free(&heap[0]);
  Zoltan_HG_heap_free(&heap[1]);

  return ZOLTAN_OK;
}



/***************************************************************************/

static int local_fmkway (
  ZZ *zz,
  HGraph *hg,
  int p,
  Partition part,
  HGPartParams *hgp,
  float bal_tol
)
{
int    i, j, vertex, edge, *cut[2], *locked = 0, *locked_list = 0, round = 0;
double total_weight, part_weight[2], max_weight, max, best_max_weight;
double cutsize, best_cutsize, *gain = 0;
float  junk;
HEAP   heap[2];
int    steplimit;
char   *yo="local_fmkway";

  if (hg->nEdge == 0)
     return ZOLTAN_OK;

  /* Calculate the different weights */
  for (i = 0; i < p; i++)
     part_weight[i] = 0.0;
  if (hg->vwgt) {
     for (i = 0; i < hg->nVtx; i++)
        part_weight[part[i]] += hg->vwgt[i];
     total_weight = 0.0;
     for (i = 0; i < p; i++)
       total_weight += part_weight[i];
     }
  else {
     total_weight = (double)(hg->nVtx);
     for (i = 0; i < hg->nVtx; i++)
        part_weight[part[i]] += 1.0;
     }
  max_weight = (total_weight / (double)p) * bal_tol;

  if (!(cut[0]      = (int*)  ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
   || !(locked      = (int*)  ZOLTAN_CALLOC    (hg->nVtx,  sizeof(int)))
   || !(locked_list = (int*)  ZOLTAN_CALLOC    (hg->nVtx,  sizeof(int)))
   || !(gain        = (double*)ZOLTAN_CALLOC    (hg->nVtx,  sizeof(double))) ) {
         Zoltan_Multifree(__FILE__,__LINE__, 4, &cut[0], &locked, &locked_list,
          &gain);
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
         }
  cut[1] = &(cut[0][hg->nEdge]);

  /* Initial calculation of the cut distribution and gain values */
  for (i = 0; i < hg->nEdge; i++)
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
        (cut[part[hg->hvertex[j]]][i])++;
  for (i = 0; i < hg->nVtx; i++)
     for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)
        {
        edge = hg->vedge[j];
        if (cut[part[i]][edge] == 1)
           gain[i] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
        else if (cut[1-part[i]][edge] == 0)
           gain[i] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
        }

  /* Initialize the heaps and fill them with the gain values */
  for (i = 0; i < p; i++)
     Zoltan_HG_heap_init(zz, &heap[i], hg->nVtx);
  for (i = 0; i < hg->nVtx; i++)
     Zoltan_HG_heap_input(&heap[part[i]], i, gain[i]);
  for (i = 0; i < p; i++)
     Zoltan_HG_heap_make(&heap[i]);

  /* Initialize given partition as best partition */
  best_cutsize = cutsize = Zoltan_HG_hcut_size_total(hg, part);
  best_max_weight = part_weight[0];
  for (i = 1; i < p; i++)
     best_max_weight = MAX(best_max_weight, part_weight[i]);

  do {
     int step = 0, no_better_steps = 0, number_locked = 0, best_locked = 0;
     int sour, dest;
     double akt_cutsize=best_cutsize;

     round++;
     cutsize = best_cutsize;
     if (hgp->output_level > HG_DEBUG_LIST)
        printf("ROUND %d:\nSTEP VERTEX  PARTS MAX_WGT CHANGE CUTSIZE\n",round);

     steplimit = (hgp->noimprove_limit > 1.1) ? (int) hgp->noimprove_limit : (int) (hg->nVtx * hgp->noimprove_limit) ;

     while (step < hg->nVtx && no_better_steps < steplimit) {
        step++;
        no_better_steps++;

        if (Zoltan_HG_heap_empty(&heap[0]))
           sour = 1;
        else if (Zoltan_HG_heap_empty(&heap[1]))
           sour = 0;
        else if (part_weight[0] > max_weight)
           sour = 0;
        else if (part_weight[1] > max_weight)
           sour = 1;
        else if (Zoltan_HG_heap_max_value(&heap[0])
              >  Zoltan_HG_heap_max_value(&heap[1]))
           sour = 0;
        else
           sour = 1;
        dest = 1-sour;

        vertex = Zoltan_HG_heap_extract_max(&heap[sour], &junk);
        locked[vertex] = part[vertex] + 1;
        locked_list[number_locked++] = vertex;
        akt_cutsize -= gain[vertex];

        Zoltan_HG_move_vertex (hg, vertex, sour, dest, part, cut, gain, heap);
        part_weight[sour] -= (hg->vwgt ? hg->vwgt[vertex] : 1.0);
        part_weight[dest] += (hg->vwgt ? hg->vwgt[vertex] : 1.0);

        max = MAX(part_weight[0],part_weight[1]);
        if ((best_max_weight > max_weight && max < best_max_weight)
         || (max<=max_weight && akt_cutsize < best_cutsize)) {
            best_locked  = number_locked;
            best_cutsize = akt_cutsize;
            best_max_weight = max;
            if (hgp->output_level > HG_DEBUG_LIST)
               printf ("New Partition:%f\n", akt_cutsize);
            no_better_steps = 0;
            }
        if (hgp->output_level > HG_DEBUG_LIST+1)
           printf ("%4d %6d %2d->%2d %7.2f %f %f\n", step, vertex, sour, dest,
            max, akt_cutsize - cutsize, akt_cutsize);
        }

     while (number_locked != best_locked) {
        vertex = locked_list[--number_locked];
        sour = part[vertex];
        dest = locked[vertex] - 1;

        Zoltan_HG_move_vertex (hg, vertex, sour, dest, part, cut, gain, heap);

        part_weight[sour] -= (hg->vwgt ? hg->vwgt[vertex] : 1.0);
        part_weight[dest] += (hg->vwgt ? hg->vwgt[vertex] : 1.0);
        Zoltan_HG_heap_input(&heap[dest], vertex, gain[vertex]);
        locked[vertex] = 0;
        }

     while (number_locked) {
        vertex = locked_list[--number_locked];
        locked[vertex] = 0;
        Zoltan_HG_heap_input(&heap[part[vertex]], vertex, gain[vertex]);
        }

     for (i = 0; i < p; i++)
        Zoltan_HG_heap_make(&(heap[i]));

     } while (best_cutsize < cutsize  ||  round < hgp->nlevelrepeat);

  /* gain_check (hg, gain, part, cut); */
  Zoltan_Multifree(__FILE__,__LINE__, 4, &cut[0], &locked, &locked_list, &gain);
  Zoltan_HG_heap_free(&heap[0]);
  Zoltan_HG_heap_free(&heap[1]);

  return ZOLTAN_OK;
}

/***************************************************************************/

/* This algorithm is loosely based on "A Coarse-Grained Parallel Formulation */
/* of Multilevel k-way Graph Partitioning Algorithm", Karypis & Kumar, 1997. */
/* It is implimented in serial as a testbed for future parallel development  */



static int local_grkway (
  ZZ *zz,
  HGraph *hg,
  int p,
  Partition part,
  HGPartParams *hgp,
  float bal_tol
)
{
const int MAX_LOOP = 5;
int     i, loop, vertex, edge, ipart;  /* loop counters */
double *part_weight, total_weight, max_weight;
float  *gain, tgain;
int   **cuts, *store1, *ends;
Vdata **lists, *store2;
int     bestpart;
char   *yo="local_grkway";

double smallest;
int found, smallpart;

  /* allocate necessary storage for heaps and weight calculation */
  if  (!(part_weight = (double*) ZOLTAN_CALLOC (p,         sizeof (double)))
   ||  !(gain        = (float*)  ZOLTAN_CALLOC (p,         sizeof (float)))
   ||  !(cuts        = (int**)   ZOLTAN_CALLOC (hg->nEdge, sizeof (int)))
   ||  !(store1      = (int*)    ZOLTAN_CALLOC (hg->nEdge * p, sizeof (int)))
   ||  !(lists       = (Vdata**) ZOLTAN_CALLOC (p,         sizeof (Vdata)))
   ||  !(store2      = (Vdata*)  ZOLTAN_CALLOC (hg->nVtx * p, sizeof (Vdata)))
   ||  !(ends        = (int*)    ZOLTAN_CALLOC (p,         sizeof (int))))
     {
     Zoltan_Multifree(__FILE__,__LINE__, 7, &part_weight, &gain, &cuts,
      &store1, &lists, &store2, &ends);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }

  /* simulate 2 dimensional arrays whose dimensions are known at run time */
  for (edge = 0; edge < hg->nEdge; edge++) cuts[edge]   = store1 + edge  * p;
  for (ipart = 0; ipart < p; ipart++)      lists[ipart] = store2 + ipart * hg->nVtx;

  /* Calculate the total weights (local vertices weight) */
  /* most general case, need MPI global communication for total weight */
  total_weight = 0.0;
  memset (part_weight, 0, p * sizeof (double));
  if (hg->vwgt)
     for (vertex = 0; vertex < hg->nVtx; vertex++)  {
         total_weight              += hg->vwgt[vertex];
         part_weight[part[vertex]] += hg->vwgt[vertex];
         }
  else {
     total_weight = hg->nVtx;
     for (vertex = 0; vertex < hg->nVtx; vertex++)
        part_weight[part[vertex]] += 1.0;
     }
  max_weight = bal_tol * total_weight / p;

  /* determine if there are any overfilled partitions */
  smallest = part_weight[0];
  smallpart = 0;
  found = 0;
  memset (ends, 0, p * sizeof (int));
  for (ipart = 0; ipart < p; ipart++)  {
     if (part_weight[ipart] < smallest)  {
        smallest = part_weight[ipart];
        smallpart = ipart;
        }
     if (part_weight[ipart] > max_weight)
        found = 1;
     }

  /* overfilled partitions found, move vertices out of it */     
  if (found)  {
     for (vertex = 0; vertex < hg->nVtx; vertex++)  {
        ipart = part[vertex];

        lists[ipart][ends[ipart]].weight      = hg->vwgt[vertex];
        lists[ipart][ends[ipart]].vertex      = vertex;
        lists[ipart][ends[ipart]].source      = ipart;
        lists[ipart][ends[ipart]].destination = -1;
        lists[ipart][ends[ipart]].gain        = 0.0;
        ++ends[ipart];
        }
     for (ipart = 0; ipart < p; ipart++)
        qsort (lists[ipart], ends[ipart], sizeof (Vdata), comparison2);

     for (ipart = 0; ipart < p; ipart++)
        for (i = 0;  (part_weight[ipart] > max_weight) && (i < ends[ipart]); i++)  {
           part_weight[ipart] -= lists[ipart][i].weight;
           part[lists[ipart][i].vertex] = smallpart;
           }
     }

  /* algorithm loops to create interprocessor communication sections */
  for (loop = 0; loop < MAX_LOOP; loop++)   {
     int oddloop = loop & 1;         /* determines direction of legal moves */
     memset (ends, 0, p * sizeof (int));

     /* Calculate the total weights (local vertices weight) */
     /* most general case, need MPI global communication for total weight */
     total_weight = 0.0;
     memset (part_weight, 0, p * sizeof (double));
     if (hg->vwgt)
        for (vertex = 0; vertex < hg->nVtx; vertex++)  {
           total_weight              += hg->vwgt[vertex];
           part_weight[part[vertex]] += hg->vwgt[vertex];
           }
     else {
        total_weight = hg->nVtx;
        for (vertex = 0; vertex < hg->nVtx; vertex++)
           part_weight[part[vertex]] += 1.0;
        }
     max_weight = bal_tol * total_weight / p;

     /* Initial calculation of the cut distribution */
     memset (store1, 0, hg->nEdge * p);
     for (edge = 0; edge < hg->nEdge; edge++)
        for (i = hg->hindex[edge]; i < hg->hindex[edge+1]; i++)
           ++cuts[edge][part[hg->hvertex[i]]];

     for (vertex = 0; vertex < hg->nVtx; vertex++)  {
        /* calculate gains */
        memset (gain, 0, p * sizeof (float));
        for (i = hg->vindex[vertex]; i < hg->vindex[vertex+1]; i++)  {
           edge = hg->vedge[i];
           for (ipart = 0; ipart < p; ipart++)  {
              if (ipart == part[vertex])
                 continue;
              if (cuts[edge][ipart] != 0 && cuts[edge][part[vertex]] == 1)
                 gain[ipart] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
              if (cuts[edge][ipart] == 0)
                 gain[ipart] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
              }
           }

        /* save best move, if any, for each vertex */
        /* oddloop, isabove, isbelow control move direction each pass */
        bestpart = -1;                            /* arbitrary illegal value */
        tgain = -1.0e37;                          /* arbitrary small value */
        for (ipart = 0; ipart < p; ipart++)
           if (ipart != part[vertex] && gain[ipart] >= tgain
            && ((!oddloop && isabove (part[vertex], ipart, p))
            ||   (oddloop && isbelow (part[vertex], ipart, p))))  {
                     if (gain[ipart] > tgain)  {
                        bestpart = ipart;
                        tgain    = gain[ipart];
                        }
                     else if (part_weight[ipart] < part_weight[bestpart])  {
                        bestpart = ipart;
                        tgain    = gain[ipart];
                        }
                     }

        /* fill heaps with the best, legal gain value per vertex */
        if (bestpart != -1 && tgain >= 0.0)  {
           lists[bestpart][ends[bestpart]].weight      = hg->vwgt[vertex];
           lists[bestpart][ends[bestpart]].vertex      = vertex;
           lists[bestpart][ends[bestpart]].source      = part[vertex];
           lists[bestpart][ends[bestpart]].destination = bestpart;
           lists[bestpart][ends[bestpart]].gain        = tgain;
           ++ends[bestpart];
           }
        } /* end of loop over all vertices */

     for (ipart = 0; ipart < p; ipart++)
        qsort (lists[ipart], ends[ipart], sizeof (Vdata), comparison);

     /* make moves until while balance is OK upto a maximum fraction */
     for (ipart = 0; ipart < p; ipart++)
        for (i = 0; i < 0.30 * ends[ipart]; i++)  {
            vertex = lists[ipart][i].vertex;
            if (hg->vwgt[vertex] + part_weight[ipart] < max_weight)  {
               part[vertex] = ipart;
               part_weight[ipart] += hg->vwgt[vertex];
               }
            }

     /* communicate back moves made/rejected, update info (faked above) */
     /* update "ghost" vertices by either the previous comm is all to all, */
     /* or by a third comm by vertex owners to ghost owners */
     /* NOTE: this too is implicit in this serial version */

     }   /* end of loop over loop */
  Zoltan_Multifree(__FILE__,__LINE__, 7, &part_weight, &gain, &cuts, &store1,
   &store2, &lists, &ends);
  return ZOLTAN_OK;
}



static int isabove (int origin, int test, int p)
   {
   if (origin > p/2)  {
      if (test < origin && test > (origin - p/2))
         return 0;
      return 1;
      }
   else  {
      if (test > origin && test <= (origin + p/2))
         return 1;
      return 0;
      }
   }



static int isbelow (int origin, int test, int p)
   {
   if (origin > p/2)  {
      if (test < origin && test >= (origin - p/2))
         return 1;
      return 0;
      }
   else  {
      if (test > origin && test < (origin + p/2))
         return 0;
      return 1;
      }
   }

   

/* note: sort is normally asending order, need the opposite!  */
static int comparison (const void *a, const void *b)
   {
   if (((Vdata*)a)->gain < ((Vdata*)b)->gain)       return  1;
   if (((Vdata*)a)->gain > ((Vdata*)b)->gain)       return -1;

   if (((Vdata*)a)->weight < ((Vdata*)b)->weight)   return  1;
   if (((Vdata*)a)->weight > ((Vdata*)b)->weight)   return -1;

   if (((Vdata*)a)->vertex > ((Vdata*)b)->vertex)   return  1;
   if (((Vdata*)a)->vertex < ((Vdata*)b)->vertex)   return -1;

   return 0;
   }



static int comparison2 (const void *a, const void *b)
   {
   if (((Vdata*)a)->weight > ((Vdata*)b)->weight)   return  1;
   if (((Vdata*)a)->weight < ((Vdata*)b)->weight)   return -1;

   return 0;
   }



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
