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



static ZOLTAN_HG_MATCHING_FN matching_mxm;  /* maximal matching */
static ZOLTAN_HG_MATCHING_FN matching_rem;  /* random edge matching */
static ZOLTAN_HG_MATCHING_FN matching_rrm;  /* random, random edge matching */
static ZOLTAN_HG_MATCHING_FN matching_rm2;  /* random, light edge matching, random */
static ZOLTAN_HG_MATCHING_FN matching_rhm;  /* random, heavy edge matching */
static ZOLTAN_HG_MATCHING_FN matching_grm;  /* greedy edge matching, Hi to Lo */
static ZOLTAN_HG_MATCHING_FN matching_gm2;  /* greedy edge matching, Lo to Hi */
static ZOLTAN_HG_MATCHING_FN matching_gm3;  /* greedy edge matching, Hi to Lo weight/size ratio */
static ZOLTAN_HG_MATCHING_FN matching_gm4;  /* greedy edge matching, Lo to Hi weight/size ratio */
static ZOLTAN_HG_MATCHING_FN matching_lhm;  /* locally heaviest matching */
static ZOLTAN_HG_MATCHING_FN matching_pgm;  /* path growing matching */
static ZOLTAN_HG_MATCHING_FN matching_aug2; /* post matching optimizer */
static ZOLTAN_HG_MATCHING_FN matching_ipm;  /* inner product matching */

static double sim (HGraph*, int, int);
void draw_matrix(HGraph*);

/* static void check_upper_bound_of_matching_weight (Graph*, ZZ*, Matching); */
/* static int graph_connected_components (int, int*, int*, int);             */
/*****************************************************************************/



int Zoltan_HG_Set_Matching_Fn(HGPartParams *hgp)
{
int found = 1;

  if      (!strcasecmp(hgp->redm_str, "mxm"))  hgp->matching = matching_mxm;
  else if (!strcasecmp(hgp->redm_str, "rem"))  hgp->matching = matching_rem;
  else if (!strcasecmp(hgp->redm_str, "rrm"))  hgp->matching = matching_rrm;
  else if (!strcasecmp(hgp->redm_str, "rm2"))  hgp->matching = matching_rm2;
  else if (!strcasecmp(hgp->redm_str, "rhm"))  hgp->matching = matching_rhm;
  else if (!strcasecmp(hgp->redm_str, "grm"))  hgp->matching = matching_grm;
  else if (!strcasecmp(hgp->redm_str, "gm2"))  hgp->matching = matching_gm2;
  else if (!strcasecmp(hgp->redm_str, "gm3"))  hgp->matching = matching_gm3;
  else if (!strcasecmp(hgp->redm_str, "gm4"))  hgp->matching = matching_gm4;
  else if (!strcasecmp(hgp->redm_str, "lhm"))  hgp->matching = matching_lhm;
  else if (!strcasecmp(hgp->redm_str, "pgm"))  hgp->matching = matching_pgm;
  else if (!strcasecmp(hgp->redm_str, "ipm"))  hgp->matching = matching_ipm;
  else if (!strcasecmp(hgp->redm_str, "no"))   hgp->matching = NULL;
  else                            { found = 0; hgp->matching = NULL;}

  if (hgp->matching) {
     /* If reduction method is a matching, set the improvement and edge weight
        scaling functions accordingly. */

     /* Note: matching_aug1 is identical to matching_mxm -> it was eliminated */
     if     (!strcasecmp(hgp->redmo_str, "aug1"))hgp->matching_opt=matching_mxm;
     else if(!strcasecmp(hgp->redmo_str, "aug2"))hgp->matching_opt=matching_aug2;
     else                                        hgp->matching_opt=NULL;
     }
  return found;
}

/****************************************************************************/



/* This is the similarity measure between two vertices in a hypergraph.
   The similarity is equal to the scaled weight of the edge in the
   transformed graph. But with this function we calculate the edge
   weights on the fly without explicitly constructing the graph. */

static double sim (HGraph *hg, int a, int b)
{
int    i, j, edge, pins, end;
double  weight, sim = 0.0;

  /* First calculate the edge weight of the graph between a and b */
  for (i = hg->vindex[a]; i < hg->vindex[a+1]; i++) {
     edge = hg->vedge[i];
     end  = hg->hindex[edge+1];
     j    = hg->hindex[edge];
     while (j < end && hg->hvertex[j] != b)
        j++;
     if (j < end) {
        pins = end - hg->hindex[edge];
        weight = 2.0 / ((pins-1) * pins);
        if (hg->ewgt)
           weight *= hg->ewgt[edge];
        sim += weight;
        }
     }
  return sim;
}



/*****************************************************************************/

int Zoltan_HG_Matching (
  ZZ *zz,
  HGraph *hg,
  Matching match,
  HGPartParams *hgp,
  int *limit)
{
float *old_ewgt = NULL, *new_ewgt = NULL;
int   err;
char  *yo = "Zoltan_HG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hg->vwgt && hgp->ews) {
     if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto End;
        }
     Zoltan_HG_Scale_HGraph_Weight (zz, hg, new_ewgt, hgp->ews);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
     }

  /* Do the matching */
  if (hgp->matching) {
     err = hgp->matching (zz, hg, match, limit);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        goto End;
     }

  /* Optimization */
  if (hgp->matching_opt)
     err = hgp->matching_opt (zz, hg, match, limit);
     
End:

  /* Restore the old edge weights */
  if (hg->vwgt && hgp->ews)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}

/*****************************************************************************/



/* maximal matching, hypergraph version */
static int matching_mxm(ZZ *zz, HGraph *hg, Matching match, int *limit)
{
int  i, j, edge, vertex, *Hindex = NULL;
char *yo = "matching_mxm";

  if (!(Hindex = (int*) ZOLTAN_MALLOC (sizeof(int) * (hg->nEdge+1)) )) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i = 0; i <= hg->nEdge; i++)
     Hindex[i] = hg->hindex[i];

  for (i = 0; i < hg->nVtx  &&  *limit > 0; i++)
     for (j = hg->vindex[i]; match[i] == i && j < hg->vindex[i+1]; j++) {
        edge = hg->vedge[j];
        while (Hindex[edge] < hg->hindex[edge+1]) {
           vertex = hg->hvertex[Hindex[edge]++];
           if (vertex != i && match[vertex] == vertex) {
              match[i]      = vertex;
              match[vertex] = i;
              (*limit)--;
              break;   /* terminate while loop */
              }
           }
        }
        
  ZOLTAN_FREE ((void**) &Hindex);
  return ZOLTAN_OK;
}

/*****************************************************************************/



/* Random edge matching */
static int matching_rem(ZZ *zz, HGraph *hg, Matching match, int *limit)
{
int i, j, *edges = NULL, edge, random, vertex, *Hindex = NULL, partner;
char *yo = "matching_rem";

   if (!(Hindex = (int*) ZOLTAN_MALLOC ((hg->nEdge+1)* sizeof(int)) )
    || !(edges  = (int*) ZOLTAN_MALLOC  (hg->nEdge   * sizeof(int))) ) {
       Zoltan_Multifree (__FILE__, __LINE__, 2, &Hindex, &edges);
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
       return ZOLTAN_MEMERR;
       }
   for (i = 0; i <= hg->nEdge; i++)
      Hindex[i] = hg->hindex[i];
   for (i = 0; i < hg->nEdge; i++)
      edges[i] = i;

   for (i = hg->nEdge; i > 0  &&  *limit > 0; i--) {
      random = Zoltan_HG_Rand() % i;
      edge = edges[random];
      edges[random] = edges[i-1];
      for (j = Hindex[edge]++; j < hg->hindex[edge+1]; j = Hindex[edge]++) {
         vertex = hg->hvertex[j];
         if (match[vertex] == vertex) {
            while (match[vertex]==vertex
             && (j = Hindex[edge]++) < hg->hindex[edge+1]) {
                random = j + Zoltan_HG_Rand() % (hg->hindex[edge+1] - j);
                partner = hg->hvertex[random];
                hg->hvertex[random] = hg->hvertex[j];
                hg->hvertex[j]      = partner;
                if (match[partner] == partner) {
                   match[vertex]  = partner;
                   match[partner] = vertex;
                   (*limit)--;
                   }
                }
            break;
            }
         }
      }
   Zoltan_Multifree (__FILE__, __LINE__, 2, &Hindex, &edges);
   return ZOLTAN_OK;
   }

/*****************************************************************************/



/* random vertex, then random edge matching */
static int matching_rrm (ZZ *zz, HGraph *hg, Matching match, int *limit)
{
int i, j, k, edge, random, vertex, pstack, partner;
int *vertices = NULL, *stack = NULL;
char *yo = "matching_rrm";

   if (!(vertices = (int*) ZOLTAN_MALLOC (hg->nVtx   * sizeof(int)))
    || !(stack    = (int*) ZOLTAN_MALLOC (hg->nInput * sizeof(int))) ) {
       Zoltan_Multifree (__FILE__, __LINE__, 2, &stack, &vertices);
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
       return ZOLTAN_MEMERR;
       }
   for (i = 0; i < hg->nVtx;  i++)
      vertices[i] = i;

   for (i = hg->nVtx; i > 0  &&  *limit > 0; i--) {
      random = Zoltan_HG_Rand() % i;
      vertex = vertices[random];
      vertices[random] = vertices[i-1];
      if (match[vertex] == vertex) {
         pstack = 0;
         for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
            edge = hg->vedge[j];
            for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; k++) {
               partner = hg->hvertex[k];
               if (partner == match[partner] && partner != vertex)
                  stack[pstack++] = partner;
               }
            }
         if (pstack > 0) {
            random = stack[Zoltan_HG_Rand() % pstack];
            match[vertex] = random;
            match[random] = vertex;
            (*limit)--;
            }
         }
      }
   Zoltan_Multifree (__FILE__, __LINE__, 2, &stack, &vertices);
   return ZOLTAN_OK;
   }

/*****************************************************************************/



/* random vertex, then lightest hyperedge weight then random vertex matching */
/* this code is almost identical with rrm (above) */
static int matching_rm2 (ZZ *zz, HGraph *hg, Matching match, int *limit)
{
int i, j, k, edge, random, vertex, pstack, partner;
int *vertices = NULL, *stack = NULL;
char *yo = "matching_rrm";

   if (!(vertices = (int*) ZOLTAN_MALLOC (hg->nVtx   * sizeof(int)))
    || !(stack    = (int*) ZOLTAN_MALLOC (hg->nInput * sizeof(int))) ) {
       Zoltan_Multifree (__FILE__, __LINE__, 2, &stack, &vertices);
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
       return ZOLTAN_MEMERR;
       }
   for (i = 0; i < hg->nVtx;  i++)
      vertices[i] = i;

   for (i = hg->nVtx; i > 0  &&  *limit > 0; i--) {
      double weight, min_weight ;                  /* added to rrm */
      random = Zoltan_HG_Rand() % i;
      vertex = vertices[random];
      vertices[random] = vertices[i-1];
      pstack = 0;
      min_weight = 1.0e20 ;                    /* added to rrm */
      if (match[vertex] == vertex) {
         for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
            edge = hg->vedge[j];

            weight = hg->ewgt ? hg->ewgt[edge] : 1.0 ;
            for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; k++) {
               partner = hg->hvertex[k];

               /* the following if-else is different from rrm */
               if (partner == match[partner] && partner != vertex
                && min_weight < weight)
                  {
                  pstack=0 ;
                  min_weight = weight ;
                  stack[pstack++] = partner;
                  }
               else if (partner == match[partner] && partner != vertex
                && min_weight == weight)
                  stack[pstack++] = partner;
               }
            }
         if (pstack > 0) {
            random = stack[Zoltan_HG_Rand() % pstack];
            match[vertex] = random;
            match[random] = vertex;
            (*limit)--;
            }
         }
      }
   Zoltan_Multifree (__FILE__, __LINE__, 2, &stack, &vertices);
   return ZOLTAN_OK;
   }




/*****************************************************************************/



/* random heavy edge matching */
static int matching_rhm (ZZ *zz, HGraph *hg, Matching match, int *limit)
{
int i, j, k, edge, random, vertex, pstack, partner, pins;
int *vertices = NULL, *stack = NULL;
double weight, max_weight, sim, *sims = NULL;
char *yo = "matching_rhm";

   if (!(vertices = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
    || !(stack    = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
    || !(sims     = (double*) ZOLTAN_CALLOC (hg->nVtx,  sizeof(double))) ) {
       Zoltan_Multifree (__FILE__, __LINE__, 3, &stack, &vertices, &sims);
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
       return ZOLTAN_MEMERR;
       }
   for (i = 0; i < hg->nVtx;  i++)
      vertices[i] = i;

   for (i = hg->nVtx; i > 0  &&  *limit > 0; i--) {
      random = Zoltan_HG_Rand() % i;
      vertex = vertices[random];
      vertices[random] = vertices[i-1];
      if (match[vertex] == vertex) {
         for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
            edge = hg->vedge[j];
            pins = hg->hindex[edge+1] - hg->hindex[edge];
            weight = 2.0 / ((pins-1)*pins);
            if (hg->ewgt)
               weight *= hg->ewgt[edge];
            for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; k++)
               if ((partner = hg->hvertex[k]) != vertex
                && match[partner] == partner)
                   sims[partner] += weight;
            }

         pstack = 0;
         max_weight = 0.0;
         for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
            edge = hg->vedge[j];
            for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; k++) {
               if ((sim = sims[(partner = hg->hvertex[k])]) > 0.0) {
                   if (sim > max_weight) {
                       pstack = 1;
                       stack[0] = partner;
                       max_weight = sim;
                       }
                   else if (sim == max_weight)
                       stack[pstack++] = partner;
             /*    else   */
             /*       continue; */    /* ignore vertex with low sim() weight */
                   sims[partner] = 0.0;
                   }
               }
            }
         if (pstack > 0) {
            random = stack[Zoltan_HG_Rand() % pstack];
            match[vertex] = random;
            match[random] = vertex;
            (*limit)--;
            }
         }
      }
   Zoltan_Multifree (__FILE__, __LINE__, 3, &stack, &vertices, &sims);
   return ZOLTAN_OK;
   }

/*****************************************************************************/



/* greedy matching, hypergraph version */
static int matching_grm (ZZ *zz, HGraph *hg, Matching pack, int *limit)
{
int   i, j, *size = NULL, *sorted = NULL, vertex;
char *yo = "matching_grm";

/* Sort the hyperedges according to their weight and size */
  if (!(size   = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int)))
   || !(sorted = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &size, &sorted);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
  for (i = 0; i < hg->nEdge; i++)
     size[i] = -(hg->hindex[i+1] - hg->hindex[i]);
  for (i = 0; i < hg->nEdge; i++)
     sorted[i] = i;
  Zoltan_quicksort_pointer_dec_float_int (sorted,hg->ewgt,size,0,hg->nEdge -1);
  ZOLTAN_FREE ((void**) &size);

  /* Match hyperedges along decreasing weight */
  for (i = 0; i < hg->nEdge  &&  *limit > 0; i++)
     for (j = hg->hindex[sorted[i]]; j < hg->hindex[sorted[i]+1]; j++)
        if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
           vertex = hg->hvertex[j];
           for (j++; j < hg->hindex[sorted[i]+1]  &&  *limit > 0; j++)
              if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                 pack[vertex]         = hg->hvertex[j];
                 pack[hg->hvertex[j]] = vertex;
                 (*limit)--;
                 break;
                 }
           break;
           }
  ZOLTAN_FREE ((void**) &sorted);
  return ZOLTAN_OK;
  }



/*****************************************************************************/

/* greedy matching, hypergraph version, matches light to heavy! */
static int matching_gm2 (ZZ *zz, HGraph *hg, Matching pack, int *limit)
{
int   i, j, *size = NULL, *sorted = NULL, vertex;
char *yo = "matching_gm2";

/* Sort the hyperedges according to their weight and size */
  if (!(size   = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int)))
   || !(sorted = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &size, &sorted);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
  for (i = 0; i < hg->nEdge; i++)
     size[i] = -(hg->hindex[i+1] - hg->hindex[i]);
  for (i = 0; i < hg->nEdge; i++)
     sorted[i] = i;
  Zoltan_quicksort_pointer_dec_float_int (sorted,hg->ewgt,size,0,hg->nEdge -1);
  ZOLTAN_FREE ((void**) &size);

  /* Match hyperedges along increasing weight */
  for (i = hg->nEdge-1; i >= 0  &&  *limit > 0; i--)
     for (j = hg->hindex[sorted[i]]; j < hg->hindex[sorted[i]+1]; j++)
        if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
           vertex = hg->hvertex[j];
           for (j++; j < hg->hindex[sorted[i]+1]  &&  *limit > 0; j++)
              if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                  pack[vertex]          = hg->hvertex[j];
                  pack[hg->hvertex[j]]  = vertex;
                  (*limit)--;
                  break;
                  }
           break;
           }
  ZOLTAN_FREE ((void**) &sorted);
  return ZOLTAN_OK;
  }

/*****************************************************************************/


/* greedy matching, hypergraph version, matching based on
   hyperdge weight divided by size. */

static int matching_gm3 (ZZ *zz, HGraph *hg, Matching pack, int *limit)
{
int   i, j, *sorted = NULL, vertex;
float *ratio = NULL;
char *yo = "matching_gm3";

/* Sort the hyperedges according to their weight/size ratio */
  if (!(ratio   = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))
   || !(sorted = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &ratio, &sorted);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
  for (i = 0; i < hg->nEdge; i++)
     ratio[i] = (hg->ewgt?hg->ewgt[i]:1.0)/(hg->hindex[i+1] - hg->hindex[i]);
  for (i = 0; i < hg->nEdge; i++)
     sorted[i] = i;
  Zoltan_quicksort_pointer_dec_float (sorted,ratio,0,hg->nEdge -1);
  ZOLTAN_FREE ((void**) &ratio);

  /* Match hyperedges along decreasing weight/size ratio  */
  for (i = 0; i < hg->nEdge  &&  *limit > 0; i++)
     for (j = hg->hindex[sorted[i]]; j < hg->hindex[sorted[i]+1]; j++)
        if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
           vertex = hg->hvertex[j];
           for (j++; j < hg->hindex[sorted[i]+1]  &&  *limit > 0; j++)
              if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                  pack[vertex]          = hg->hvertex[j];
                  pack[hg->hvertex[j]]  = vertex;
                  (*limit)--;
                  break;
                  }
           break;
           }
  ZOLTAN_FREE ((void**) &sorted);
  return ZOLTAN_OK;
  }

/*****************************************************************************/


/* greedy matching, hypergraph version, matching based on
   hyperdge weight divided by size. Reverse order of gm3. */

static int matching_gm4 (ZZ *zz, HGraph *hg, Matching pack, int *limit)
{
int   i, j, *sorted = NULL, vertex;
float *ratio = NULL;
char *yo = "matching_gm4";

/* Sort the hyperedges according to their weight/size ratio */
  if (!(ratio   = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))
   || !(sorted = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &ratio, &sorted);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
  for (i = 0; i < hg->nEdge; i++)
     ratio[i] = (hg->ewgt?hg->ewgt[i]:1.0)/(hg->hindex[i+1] - hg->hindex[i]);
  for (i = 0; i < hg->nEdge; i++)
     sorted[i] = i;
  Zoltan_quicksort_pointer_dec_float (sorted,ratio,0,hg->nEdge -1);
  ZOLTAN_FREE ((void**) &ratio);

  /* Match hyperedges along increasing weight/size ratio  */
  for (i = hg->nEdge-1; i >= 0  &&  *limit > 0; i--)
     for (j = hg->hindex[sorted[i]]; j < hg->hindex[sorted[i]+1]; j++)
        if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
           vertex = hg->hvertex[j];
           for (j++; j < hg->hindex[sorted[i]+1]  &&  *limit > 0; j++)
              if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                  pack[vertex]          = hg->hvertex[j];
                  pack[hg->hvertex[j]]  = vertex;
                  (*limit)--;
                  break;
                  }
           break;
           }
  ZOLTAN_FREE ((void**) &sorted);
  return ZOLTAN_OK;
  }

/*****************************************************************************/

#undef LAM_ORIG
/* helper function for locally heavy matching (matching_lhm) below */
static int lhm_match (int a, int b, float ewgt_ab, int *nindex, int *neigh,
 float *ewgt, int *Nindex, int *limit, int *match)
{
int Nindex_a_old = Nindex[a], Nindex_b_old = Nindex[b], c, c_deg,
 a_deg = nindex[a+1] - nindex[a], b_deg = nindex[b+1] - nindex[b];
float c_ewgt;

  while (match[a] == a && match[b] == b &&
   (Nindex[a] < nindex[a+1] || Nindex[b] < nindex[b+1]) && (*limit) > 0 ) {
      if (Nindex[a] < nindex[a+1]) {
         c = neigh[Nindex[a]];
         c_ewgt = ewgt[Nindex[a]];
         c_deg = nindex[c+1] - nindex[c];
         (Nindex[a])++;
#ifdef LAM_ORIG
         if (c_ewgt > ewgt_ab)
#else
         if ((b_deg > 1 && ((c_deg > 1 && (c_ewgt > ewgt_ab || (c_ewgt == ewgt_ab
          && c_deg < b_deg))) || (c_deg == 1 && 2.0 * c_ewgt >= ewgt_ab)))
          || (b_deg == 1 && ((c_deg  > 1 && c_ewgt > 2.0 * ewgt_ab)
          || (c_deg == 1 &&   c_ewgt > ewgt_ab)) ) )
#endif
         lhm_match(a, c, c_ewgt, nindex, neigh, ewgt, Nindex, limit, match);
         }
      if (match[b] == b  &&  Nindex[b] < nindex[b+1]  && (*limit) > 0) {
         c = neigh[Nindex[b]];
         c_ewgt = ewgt[Nindex[b]];
         c_deg = nindex[c+1] - nindex[c];
         (Nindex[b])++;
#ifdef LAM_ORIG
         if (c_ewgt > ewgt_ab)
#else
         if ((a_deg > 1 && ((c_deg > 1 && (c_ewgt > ewgt_ab || (c_ewgt == ewgt_ab
          && c_deg < a_deg))) || (c_deg == 1 && 2.0 * c_ewgt >= ewgt_ab)))
          || (a_deg == 1 && ((c_deg > 1 && c_ewgt > 2.0 * ewgt_ab)
          || (c_deg == 1 && c_ewgt > ewgt_ab)) ))
#endif
          lhm_match(b, c, c_ewgt, nindex, neigh, ewgt, Nindex, limit, match);
          }
      }

  if (match[a] == a  &&  match[b] == b  &&  *limit > 0) {
     match[a] = b;
     match[b] = a;
     (*limit)--;
     }
  else if (match[a] == a  &&  match[b] != b)
     Nindex[a] = Nindex_a_old;
  else if (match[a] != a  &&  match[b] == b)
     Nindex[b] = Nindex_b_old;

  return ZOLTAN_OK;
}


/* Locally Heaviest matching. 
   It looks for locally heaviest edges by backtracking. For a current
   edge it has to check the weight of all intersecting edges. Edges
   may be checked several times, but it is amortized not more than
   time O(k*|I|) and guarantees an approximation of 1/k. */
   

/* locally heavy matching, graph version                                      */
/* This is an implimentation of LAM created by Robert Preis, Linear Time      */
/* 1/2-Approximation Algorithm for Maximum Weighted Matching in General Graphs*/
/* Symposium on Theoretical Aspects of Computer Science, STACS 99, C. Meinel, */
/* S. Tison (eds.), Springer, LNCS 1563, 1999, 259-269                        */
static int matching_lhm (ZZ *zz, HGraph *hg, Matching match, int *limit)
{
int  i, j, *Nindex = NULL, err;
Graph g;
char *yo = "matching_lhm";

     Zoltan_HG_Graph_Init (&g);
     err = Zoltan_HG_HGraph_to_Graph (zz, hg, &g);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
        Zoltan_HG_Graph_Free(&g);
        return err;
        }

  if (!(Nindex = (int*) ZOLTAN_MALLOC ((hg->nVtx+1) * sizeof(int))))
     {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  memcpy(Nindex, g.nindex, (g.nVtx+1) * sizeof(int));

  for (i = 0; i < hg->nVtx  &&  *limit > 0; i++)
     for (j = g.nindex[i]; match[i] == i && j < g.nindex[i+1]; j++)
        if (match[g.neigh[j]] == g.neigh[j])
           lhm_match(i, g.neigh[j], g.ewgt[j], g.nindex, g.neigh, g.ewgt,
            Nindex, limit, match);

  ZOLTAN_FREE ((void**) &Nindex);
  Zoltan_HG_Graph_Free(&g);
  return ZOLTAN_OK;
}

/*****************************************************************************/



/* path growing matching, hypergraph version */
static int matching_pgm (ZZ *zz, HGraph *hg, Matching match, int *limit)
{
int i, j, k, side = 0, edge, vertex, *Match[2], limits[2], neighbor,
 next_vertex, pins;
double w[2]={0.0,0.0}, weight, max_weight, *sims;
char  *yo = "matching_pgm";

  limits[0] = limits[1] = *limit;
  Match[0] = match;
  if (!(Match[1] = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
   || !(sims     = (double*) ZOLTAN_CALLOC (hg->nVtx,  sizeof(double))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &Match[1], &sims);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
  for (i = 0; i < hg->nVtx; i++)
     Match[1][i] = i;
     
  for (i = 0; i < hg->nVtx  &&  limits[side] > 0; i++)
     if (Match[0][i] == i && Match[1][i] == i) {
        vertex = i;
        while (vertex > 0 && limits[side] > 0) {
           max_weight = 0.0;
           next_vertex = -1;

           for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
              edge = hg->vedge[j];
              pins = hg->hindex[edge+1] - hg->hindex[edge];
              weight = 2.0 / ((pins-1)*pins);
              if (hg->ewgt)
                 weight *= hg->ewgt[edge];
              for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; k++) {
                 neighbor = hg->hvertex[k];
                 if (neighbor != vertex && Match[0][neighbor] == neighbor
                  && Match[1][neighbor]==neighbor)
                     sims[neighbor] += weight;
                 }
              }

        for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
           edge = hg->vedge[j];
           for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; k++) {
              neighbor = hg->hvertex[k];
              if (sims[neighbor] > 0.0) {
                 if (sims[neighbor] > max_weight) {
                    max_weight = sims[neighbor];
                    next_vertex = neighbor;
                    }
                 sims[neighbor] = 0.0;
                 }
              }
           }

        if (next_vertex >= 0) {
           Match[side][vertex] = next_vertex;
           Match[side][next_vertex] = vertex;
           limits[side]--;
           w[side] += max_weight;
           side = 1-side;
           }
        vertex = next_vertex;
        }
     }

  if (w[0] < w[1]) {
     for (i = 0; i < hg->nVtx; i++)
        match[i] = Match[1][i];
     *limit = limits[1];
     }
  else
     *limit = limits[0];

  Zoltan_Multifree (__FILE__, __LINE__, 2, &Match[1], &sims);
  return ZOLTAN_OK;
}

/*****************************************************************************/

/* inner product matching                                                    */ 
/* based on Rob Bisseling's implementation in Mondriaan                      */
/* for each vertex, we match with the unmatched vertex which has the most    */
/* hyperedges in common with it (ie, the pair with greatest inner product).  */

static int matching_ipm(ZZ *zz, HGraph *hg, Matching match, int *limit)
{
    int   i, j, k, v1, v2, edge, ip, maxip, maxindex;
    int   *checked, *vector;
    char  *yo = "matching_ipm";

    if (!(checked = (int*) ZOLTAN_MALLOC (hg->nVtx  * sizeof(int))) 
     || !(vector  = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int)))) {
        Zoltan_Multifree (__FILE__, __LINE__, 2, &checked, &vector);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }
    
    //j = 0;
    //
    //printf("Debugging info:\nvedges:\n----------\n");
    //for(i = 0; i < hg->vindex[hg->nVtx]; ++i) {
    //    printf("%2d ", hg->vedge[i]);
    //    if(hg->vindex[j+1] == i + 1) {
    //        j++;
    //        printf("\n");
    //    }
    //}
    //printf("\nvindex:\t");
    //for(i = 0; i < hg->nVtx + 1; ++i)
    //    printf("%2d ", hg->vindex[i]);
    //printf("\n\nhvertex\n-------\n");
    //j = 0;
    //for(i = 0; i < hg->hindex[hg->nEdge]; ++i) {
    //    printf("%2d ", hg->hvertex[i]);
    //    if(hg->hindex[j+1] == i + 1) {
    //        j++;
    //        printf("\n");
    //    }
    //}
    //printf("\nhindex:\t");
    //for(i = 0; i < hg->nEdge + 1; ++i)
    //    printf("%2d ", hg->hindex[i]);
    //
    //printf("\n");
    //draw_matrix(hg);
    //printf("\n");
        
    /*
     * Next optimization: Vertex switching
     * When we check a vertex, we move it to the end of the data structure,
     * rather than marking a checked array.  We then avoid the overhead of
     * misses when we check the checked array.
     *
     * To do this, introduce another array, indices into the vindex array.
     * Initially each element is its index, but is swapped with the last
     * element whenever its vertex is checked.  We thereby avoid actual
     * modification of the HGraph structure.  A similar procedure is used
     * in the randomized algorithms above.
     */
     
    for(i = 0; i < hg->nVtx; i++)
        checked[i] = -1;
    for(i = 0; i < hg->nEdge; i++)
        vector[i] = -1;
    
    /* for every vertex */
    for (v1 = 0; v1 < hg->nVtx  &&  *limit > 0; v1++) {
        if(match[v1] != v1)
            continue;
        
        maxip = 0;
        maxindex = -1;
        for(i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++)
            vector[hg->vedge[i]] = v1;
            
        /* for every hyperedge containing the vertex */
        for(i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
            edge = hg->vedge[i];
                
            /* for every other vertex in the hyperedge */
            for(j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
                v2 = hg->hvertex[j];
                
                /* ignore matched and previously checked vertices */
                if(match[v2] != v2 || checked[v2] == v1 || v1 == v2)
                    continue;
            
                /* compute inner product of vertex v1 with vertex v2 */
                ip = 0;
                for(k = hg->vindex[v2]; k < hg->vindex[v2+1]; k++) {
                    if(vector[hg->vedge[k]] == v1)
                        ip++;
                }
                checked[v2] = v1;
                
                /* keep track of best match seen */
                if(maxip < ip) {
                    maxip = ip;
                    maxindex = v2;
                }
                //printf("IP of %d with %d is %d\n", v1, v2, ip);
            }
        }

        //printf("Done with %d, best match is %d with product %d\n",
        //        v1, maxindex, maxip);

        /* match v1 with v2 having greatest inner product */
        if(maxindex != -1) {
            match[v1] = maxindex;
            match[maxindex] = v1;
            (*limit)--;
        } 
    }

    //printf("Final Matching:\n");
    //for(i = 0; i < hg->nVtx; i++)
    //    printf("%2d ",i);
    //printf("\n");
    //for(i = 0; i < hg->nVtx; i++)
    //    printf("%2d ",match[i]);
    //printf("\n");

    Zoltan_Multifree (__FILE__, __LINE__, 2, &checked, &vector);
    return ZOLTAN_OK;
}

/*****************************************************************************/



/* optimization passed on augmenting path of length 2 */
static int matching_aug2 (ZZ *zz, HGraph *hg, Matching match, int *limit)
{
int   i, j, k, *stack, free_p=0, edge, vertex, best_2;
double gain, gain_2;
char  *yo = "matching_aug2";

  if (!(hg->ewgt))
     return matching_mxm (zz, hg, match, limit);

  if (!(stack = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i = 0; i < hg->nVtx; i++)
     if (match[i] == i)
        stack[free_p++] = i;

  while (free_p  &&  *limit > 0) {
     i = stack[--free_p];
     if (match[i] == i) {
        gain_2 = EPS;
        best_2 = -1;

        for (j = hg->vindex[i]; match[i] == i && j < hg->vindex[i+1]; j++) {
           edge = hg->vedge[j];
           for (k = hg->hindex[edge]; match[i] == i && k < hg->hindex[edge+1];
            k++)
              if ((vertex = hg->hvertex[k]) != i) {
                 if (match[vertex] == vertex) {
                    match[i]      = vertex;
                    match[vertex] = i;
                    (*limit)--;
                    }
                 else if ((gain = sim(hg, i, vertex)
                  - sim(hg, vertex, match[vertex])) > (1.0+EPS) * gain_2) {
                      gain_2 = gain;
                      best_2 = vertex;
                      }
                  }
           }

      if (match[i] == i && best_2 >= 0)  {
         stack[free_p++]      = match[best_2];
         match[match[best_2]] = match[best_2];
         match[best_2]        = i;
         match[i]             = best_2;
         }
      }
  }
  ZOLTAN_FREE ((void**) &stack);
  return ZOLTAN_OK;
}

/*****************************************************************************/



/* Old version of aug3 which runs on the graph. I am not sure if it
   is worth while to transfer it to only use the hypergraph and
   the sim-function. If there will be time, one could do that. */
/*
static int matching_aug3 (ZZ *zz, HGraph *hg, Graph *g, Matching match,
 int *limit)
{
int i, j, k, *stack, free_p = 0, best_2 = -1, best_near = -1,
int best_middle = -1, best_distant = -1, neigh1, neigh2, neigh3;
double w_1, w_2, w_3, gain_2, gain_3;
char   *yo = "matching_aug3";

  if (!(stack = (int*) ZOLTAN_MALLOC (g->nVtx * sizeof(int)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i = 0; i < g->nVtx; i++)
     if (match[i] == i)
        stack[free_p++] = i;

  while (free_p && (*limit) > 0) {
     i = stack[--free_p];
     if (match[i] == i) {
        gain_2 = 0.0;
        gain_3 = -FLT_MIN;
        for (j = g->nindex[i]; match[i] == i && j < g->nindex[i+1]; j++) {
           neigh1 = g->neigh[j];
           w_1 = (double)(g->ewgt ? g->ewgt[j] : 1.0);
           if (match[neigh1] == neigh1) {
              match[i] = neigh1;
              match[neigh1] = i;
              (*limit)--;
              }
           else {
              neigh2 = match[neigh1];
              k = g->nindex[neigh1];
              while (neigh2 != g->neigh[k])
                 k++;
              w_2 = (double)(g->ewgt ? g-> ewgt[k] : 1.0);
              if (w_1 - w_2 > gain_2) {
                 gain_2 = w_1 - w_2;
                 best_2 = neigh1;
                 }
           for (k = g->nindex[neigh2]; k < g->nindex[neigh2+1]; k++) {
              neigh3 = g->neigh[k];
              if (match[neigh3] == neigh3 && neigh3 != i) {
                 w_3 = (double)(g->ewgt ? g->ewgt[k] : 1.0);
                 if (w_1 - w_2 + w_3 > gain_3) {
                    best_near = neigh1;
                    best_middle = neigh2;
                    best_distant = neigh3;
                    gain_3 = w_1 - w_2 + w_3;
                    }
                 }
              }
           }
        }

      if (match[i] == i) {
         if (gain_3 > -FLT_MIN) {
            match[i] = best_near;
            match[best_near] = i;
            match[best_middle] = best_distant;
            match[best_distant] = best_middle;
            (*limit)--;
            }
         else if (gain_2 > 0.0) {
            match[i] = best_2;
            stack[free_p++] = match[best_2];
            match[match[best_2]] = match[best_2];
            match[best_2] = i;
            }
         }
      }
  }
  ZOLTAN_FREE ((void **) &stack);
  return ZOLTAN_OK;
}
*/

/*****************************************************************************/



/*
static void check_upper_bound_of_matching_weight(Graph *g,ZZ *zz,Matching match)
{ int i, j;
  double highest, matching_weight=0.0, upper_bound=0.0;

  for (i = 0; i < g->nVtx; i++) {
     highest = 0.0;
     for (j = g->nindex[i]; j < g->nindex[i+1]; j++) {
        highest = MAX(highest, g->ewgt[j]);
        if (g->neigh[j] == match[i]) {
           if (g->ewgt[j] == FLT_MAX)
              matching_weight = FLT_MAX;
           else
              matching_weight += (g->ewgt[j]);
           }
        }
     if (highest < FLT_MAX)
        upper_bound += highest;
     else
        upper_bound = FLT_MAX;
     }
  printf("matching/upper_bound ratio:%.3f\n", matching_weight / upper_bound);
}



static int graph_connected_components (int n, int *ep, int *edge, int Out)
{ int i=0, *locked, *queue, queue_head, queue_tail,
  int components=0, first_unlocked=0;

  if (!(queue  = (int*) ZOLTAN_MALLOC (n*sizeof(int)))
   || !(locked = (int*) ZOLTAN_CALLOC (n,sizeof(int))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &locked, &queue);
      return ZOLTAN_MEMERR;
      }

  if (Out > 1)
     printf ("Component Sizes     : ");
  while (i < n) {
     int k, size=1, vertex;
     queue_head = 0;
     queue_tail = 1;
     components++;
     while (locked[first_unlocked])
        first_unlocked++;
     queue[0] = first_unlocked;
     locked[first_unlocked] = 1;
     while (queue_head < queue_tail) {
        vertex = queue[queue_head++];
        for (k = ep[vertex]; k < ep[vertex+1]; k++)
           if (!(locked[edge[k]])) {
              size++;
              locked[edge[k]] = 1;
              queue[queue_tail++] = edge[k];
              }
        }
     if (Out > 1)
        printf ("%d ", size);
     i += size;
     }
  if (Out > 1)
     printf ("\n");
  if (Out > 0)
     printf("Components          : %d\n",components);
  Zoltan_Multifree (__FILE__, __LINE__, 2, &locked, &queue);
  return components;
}
*/



/*****************************************************************************/

/* 
 * draws the sparse matrix representation of a hypergraph, with vertices as
 * rows and columns as edges.
 */
void draw_matrix(HGraph* hg)
{
    int position, i, j;
    
    for(i = 0; i < hg->nVtx; ++i) {
        j = position = 0;
        while(j < hg->nEdge) {
            if(hg->vedge[hg->vindex[i] + position] == j) {
                position++;
                printf("x ");
            }
            else
                printf("0 ");
            j++;
        }
        printf("\n");
    }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

