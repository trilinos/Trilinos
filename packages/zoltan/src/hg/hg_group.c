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



static ZOLTAN_HG_GROUPING_FN grouping_mxg;  /* maximal grouping */
static ZOLTAN_HG_GROUPING_FN grouping_reg;  /* random edge grouping */
static ZOLTAN_HG_GROUPING_FN grouping_rrg;  /* random, random grouping */
static ZOLTAN_HG_GROUPING_FN grouping_rhg;  /* random, heavy grouping */
static ZOLTAN_HG_GROUPING_FN grouping_grg;  /* greedy grouping */
static ZOLTAN_HG_GROUPING_FN grouping_aug2; /* augmenting path; length 2 */

/****************************************************************************/



int Zoltan_HG_Set_Grouping_Fn(HGPartParams *hgp)
{
int found = 1;

  if      (!strcasecmp(hgp->redm_str,"mxg")) hgp->grouping = grouping_mxg;
  else if (!strcasecmp(hgp->redm_str,"reg")) hgp->grouping = grouping_reg;
  else if (!strcasecmp(hgp->redm_str,"rrg")) hgp->grouping = grouping_rrg;
  else if (!strcasecmp(hgp->redm_str,"rhg")) hgp->grouping = grouping_rhg;
  else if (!strcasecmp(hgp->redm_str,"grg")) hgp->grouping = grouping_grg;
  else if (!strcasecmp(hgp->redm_str, "no")) hgp->grouping = NULL;
  else                          { found = 0; hgp->grouping = NULL; }

  if (hgp->grouping) {
     /* If reduction method is a grouping, set the improvement and edge weight
      * scaling functions accordingly */

     /* register optimization function */
     if     (!strcasecmp(hgp->redmo_str,"aug1"))hgp->grouping_opt=grouping_mxg;
     else if(!strcasecmp(hgp->redmo_str,"aug2"))hgp->grouping_opt=grouping_aug2;
     else                                       hgp->grouping_opt=NULL;
     }

  return found;
}

/****************************************************************************/



int Zoltan_HG_Grouping (ZZ *zz, HGraph *hg, Packing pack, HGPartParams *hgp,
 int *limit)
{
int   err = ZOLTAN_OK;
float *old_ewgt = NULL, *new_ewgt = NULL;
char  *yo = "Zoltan_HG_Grouping";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hg->vwgt && hgp->ews) {
     if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
        }
     Zoltan_HG_Scale_HGraph_Weight (zz, hg, new_ewgt, hgp->ews);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
     }

  /* Do the grouping */
  if (hgp->grouping) {
     err = hgp->grouping(zz, hg, pack, limit);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        goto End;
     }

  /* Optimization */
  if (hgp->grouping_opt != NULL)
     err = hgp->grouping_opt (zz, hg, pack, limit);

End:
  /* Restore the old edge weights */
  if (hg->vwgt && hgp->ews)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}

/****************************************************************************/



/* Maximal grouping. Just goes through all edges and groups whatever is
   available. Time O(|I|) */
static int grouping_mxg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
   {
int i, j, vertex, first_vertex;

   for (i = 0; i < hg->nEdge  &&  *limit > 0; i++)
      for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++) {
         vertex = hg->hvertex[j];
         if (pack[vertex] == vertex) {
            first_vertex   = vertex;
            for (j++; j < hg->hindex[i+1]  &&  *limit > 0; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                  vertex = pack[vertex]  = hg->hvertex[j];
                  (*limit)--;
                  }
            pack[vertex] = first_vertex;
            }
         }
   return ZOLTAN_OK;
   }

/****************************************************************************/



/* Random Egde Grouping. Randomly traverses through all edges and
   groups whatever is available in that edge. Time O(|I|)  */
static int grouping_reg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
int i, j, *edges = NULL, edge, random, vertex, first_vertex;
char *yo = "grouping_reg";

   if (!(edges = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0; i < hg->nEdge; i++)
      edges[i] = i;

   for (i = hg->nEdge; i > 0  &&  *limit > 0; i--) {
      random = Zoltan_HG_Rand() % i;
      edge = edges[random];
      edges[random] = edges[i-1];

      j = hg->hindex[edge];
      while (j < hg->hindex[edge+1]) {
         vertex = hg->hvertex[j];
         if (pack[vertex] == vertex) {
            first_vertex = vertex;
            j++;
            while (j < hg->hindex[edge+1]  &&  *limit > 0) {
               if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                  vertex = pack[vertex] = hg->hvertex[j];
                  (*limit)--;
                  }
               j++;
               }
            pack[vertex] = first_vertex;
            }
         j++;
         }
      }
   ZOLTAN_FREE ((void**) &edges);
   return ZOLTAN_OK;
   }

/****************************************************************************/



/* Random Random Grouping. Randomly traverses through the vertices
   and for each vertex it randomly takes one covering edge and groups
   whatever is available. Time O(|I|) */
static int grouping_rrg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
int i, j, edge, random, *vertices = NULL, vertex, first_vertex, count;
char *yo = "grouping_rrg";

   if (!(vertices  = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0; i < hg->nVtx;  i++)
      vertices[i] = pack[i] = i;

   for (i = hg->nVtx; i > 0  &&  *limit > 0; i--) {
      random = Zoltan_HG_Rand() % i;
      vertex = vertices[random];
      vertices[random] = vertices[i-1];
      if (pack[vertex] != vertex)
         continue;          /* vertex already packed, move on */

      count = hg->vindex[vertex+1] - hg->vindex[vertex];  /* count edges */
      if (count == 0)
         continue;

      edge = hg->vedge[hg->vindex[vertex] + (Zoltan_HG_Rand() % count)];  /* random edge */
      for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
            first_vertex = vertex = hg->hvertex[j];
            for (j++; j < hg->hindex[edge+1]  &&  *limit > 0; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                  vertex = pack[vertex] = hg->hvertex[j];
                  (*limit)--;
                  }
            pack[vertex] = first_vertex;
            break;
            }
      }
   ZOLTAN_FREE ((void **) &vertices);
   return ZOLTAN_OK;
   }

/****************************************************************************/



/* Random Heavy Grouping. Randomly traverses through the vertices
   and for each vertex it takes the heaviest covering edge and groups
   whatever is available. Time O(|I|) */
static int grouping_rhg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
int   i, j, *vertices = NULL, *del_edges = NULL, vertex, first_vertex, edge,
 number, best_edge, best_size, best_neighbors, random, size;
float best_ewgt;
char  *yo = "grouping_rhg";

   if (!(vertices  = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
    || !(del_edges = (int*) ZOLTAN_CALLOC (hg->nEdge, sizeof(int)))) {
       Zoltan_Multifree (__FILE__, __LINE__, 2, &vertices, &del_edges);
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
       return ZOLTAN_MEMERR;
       }
   for (i = 0; i < hg->nVtx; i++)
      vertices[i] = i;

   for (i = hg->nVtx; i > 0  &&  *limit > 0; i--) {
      number = Zoltan_HG_Rand() % i;
      vertex = vertices[number];
      vertices[number] = vertices[i-1];
      if (pack[vertex] != vertex)
         continue;            /* vertex is already matched, move on */

      best_neighbors = 0;
      best_edge = best_size = -1;
      best_ewgt = -1.0;
      for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
         int size;
         edge = hg->vedge[j];
         size = hg->hindex[edge+1] - hg->hindex[edge];
         if (del_edges[edge] == 0 && ((!(hg->ewgt) && size == best_size)
          || (hg->ewgt && hg->ewgt[edge] == best_ewgt && size == best_size)))
              best_neighbors++;
         else if (del_edges[edge] == 0 && ((!(hg->ewgt) && size < best_size)
          || (hg->ewgt && (hg->ewgt[edge] >  best_ewgt
          || (hg->ewgt[edge] == best_ewgt && size < best_size))))) {
              best_neighbors = 1;
              best_edge = edge;
              best_ewgt = hg->ewgt[edge];
              best_size = hg->hindex[edge+1] - hg->hindex[edge];
              }
         }
      if (best_neighbors == 0)
         continue;                       /* no suitable edge found */

      if (best_neighbors > 1) {
         random = Zoltan_HG_Rand() % best_neighbors;
         for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
            edge = hg->vedge[j];
            size = hg->hindex[edge+1] - hg->hindex[edge];
            if (del_edges[edge]==0 && ((!(hg->ewgt)       && size == best_size)
              || (hg->ewgt && hg->ewgt[edge] == best_ewgt && size == best_size))
              && --best_neighbors == random) {
                 best_edge = edge;
                 break;
                 }
            }
         }

      for (j = hg->hindex[best_edge]; j < hg->hindex[best_edge+1]; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
            first_vertex = vertex = hg->hvertex[j];
            for (j++; j < hg->hindex[best_edge+1] && (*limit) > 0; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                  vertex = pack[vertex] = hg->hvertex[j];
                  (*limit)--;
                  }
            pack[vertex] = first_vertex;
            break;
            }
      del_edges[best_edge] = 1;
      }
   Zoltan_Multifree (__FILE__, __LINE__, 2, &vertices, &del_edges);
   return ZOLTAN_OK;
}

/****************************************************************************/



/* Greedy Grouping. It sorts the edges due to their weight and size
   and groups whatever is available. Time O(|I|+|E|*log(|E|)) */
static int grouping_grg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
int   i, j, *size = NULL, *sorted = NULL, first_vertex, vertex;
char *yo = "grouping_grg";

  /* Sort the hyperedges according to their weight and size */
  if (!(size   = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int)))
   || !(sorted = (int*) ZOLTAN_MALLOC (hg->nEdge * sizeof(int)))) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &size, &sorted);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
  for (i = 0; i < hg->nEdge; i++) /* negate sizes */
      size[i] = -(hg->hindex[i+1] - hg->hindex[i]);
  for (i = 0; i < hg->nEdge; i++)
      sorted[i] = i;
  Zoltan_quicksort_pointer_dec_float_int(sorted, hg->ewgt, size,0,hg->nEdge-1);
  ZOLTAN_FREE ((void**) &size);

  /* Match hyperedges along decreasing weight */
  for (i = 0; i < hg->nEdge  &&  *limit > 0; i++)
     for (j = hg->hindex[sorted[i]]; j < hg->hindex[sorted[i]+1]; j++)
        if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
           first_vertex = vertex = hg->hvertex[j];
           for (j++; j < hg->hindex[sorted[i]+1]  &&  *limit > 0; j++)
              if (pack[hg->hvertex[j]] == hg->hvertex[j]) {
                 vertex = pack[vertex] = hg->hvertex[j];
                 (*limit)--;
                 }
           pack[vertex] = first_vertex;
           break;
           }
  ZOLTAN_FREE ((void**) &sorted);
  return ZOLTAN_OK;
  }

/****************************************************************************/



static int grouping_aug2 (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
/* Placeholder for grouping_aug2. */
  return grouping_mxg (zz, hg, pack, limit);
}

/****************************************************************************/



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
