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
static ZOLTAN_HG_GROUPING_FN grouping_aug1; /* augmenting path; length 1 */
static ZOLTAN_HG_GROUPING_FN grouping_aug2; /* augmenting path; length 2 */
static ZOLTAN_HG_GROUPING_FN grouping_aug3; /* augmenting path; length 3 */

/****************************************************************************/

int Zoltan_HG_Set_Grouping_Fn(HGPartParams *hgp)
{
  static int srand_set = 0;

  if (srand_set == 0)
     {
     srand_set = 1 ;
     srand ((unsigned long) RANDOM_SEED) ;
     }


  if      (!strcasecmp(hgp->redm_str, "mxg"))  hgp->grouping = grouping_mxg;
  else if (!strcasecmp(hgp->redm_str, "reg"))  hgp->grouping = grouping_reg;
  else if (!strcasecmp(hgp->redm_str, "rrg"))  hgp->grouping = grouping_rrg;
  else if (!strcasecmp(hgp->redm_str, "rhg"))  hgp->grouping = grouping_rhg;
  else if (!strcasecmp(hgp->redm_str, "grg"))  hgp->grouping = grouping_grg;
  else                                         hgp->grouping = NULL;

  if (hgp->grouping) {
  /*
   * If reduction method is a grouping, set the improvement and
   * edge weight scaling functions accordingly.
   */
    if      (!strcasecmp(hgp->redmo_str,"aug1")) hgp->grouping_opt = grouping_aug1;
    else if (!strcasecmp(hgp->redmo_str,"aug2")) hgp->grouping_opt = grouping_aug2;
    else if (!strcasecmp(hgp->redmo_str,"aug3")) hgp->grouping_opt = grouping_aug3;
    else                                         hgp->grouping_opt = NULL;
  }

  return hgp->grouping ? 1 : 0 ;
}

/****************************************************************************/

int Zoltan_HG_Grouping (ZZ *zz, HGraph *hg, Packing pack, HGPartParams *hgp, int *limit)
{ int   ierr = ZOLTAN_OK;
  float *old_ewgt=NULL, *new_ewgt;
  char  *yo = "Zoltan_HG_Grouping";

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (hg->vwgt && hgp->ews)
  { if (!(new_ewgt = (float *) ZOLTAN_MALLOC (sizeof (float) * hg->nEdge)))
    { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    Zoltan_HG_Scale_HGraph_Weight (zz, hg, new_ewgt);
    old_ewgt = hg->ewgt;
    hg->ewgt = new_ewgt;
  }

  ierr = hgp->grouping(zz,hg,pack,limit);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    goto End;
  }

  /* Optimization */
  if (hgp->grouping_opt != NULL)
    ierr = hgp->grouping_opt (zz,hg,pack,limit);

End:
  if (hg->vwgt && hgp->ews)
  { hg->ewgt = old_ewgt;
    ZOLTAN_FREE ((void **) &new_ewgt);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/****************************************************************************/

static int grouping_mxg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
   {
   int i, j, vertex, first_vertex ;

   for (i = 0 ; i < hg->nVtx ; i++)
      pack[i] = i ;

   for (i = 0 ; i < hg->nEdge && (*limit)>0 ; i++)
      for (j = hg->hindex[i] ; j < hg->hindex[i+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[i+1] && (*limit)>0 ; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  {
                  vertex = pack[vertex] = hg->hvertex[j] ;
                  (*limit)--;
                  }
            pack[vertex] = first_vertex ;
            break ;       /* not required, might improve speed */
            }
   return ZOLTAN_OK ;
   }

/****************************************************************************/

static int grouping_reg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
   {
   int i, j, *edges=NULL, edge, random, vertex, first_vertex ;
   char *yo = "grouping_reg" ;

   if (!(edges = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge)))
      {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nEdge ; i++)
      edges[i] = i ;
   for (i = 0 ; i < hg->nVtx ;  i++)
      pack[i]  = i ;

   for (i = hg->nEdge ; i > 0 && (*limit)>0 ; i--)
      {
      random = rand() % i ;
      edge = edges[random] ;
      edges[random] = edges[i-1] ;

      for (j = hg->hindex[edge] ; j < hg->hindex[edge+1]; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[edge+1] && (*limit)>0 ; j++)
                if (pack[hg->hvertex[j]] == hg->hvertex[j])
                   {
                   vertex = pack[vertex] = hg->hvertex[j] ;
                   (*limit)--;
                   }
            pack[vertex] = first_vertex ;
            break ;
            }
      }
   ZOLTAN_FREE ((void **) &edges) ;
   return ZOLTAN_OK ;
   }

/****************************************************************************/

static int grouping_rrg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
   {
   int i, j, edge, random, *vertices=NULL, vertex, first_vertex, count ;
   char *yo = "grouping_rrg" ;

   if (!(vertices  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
      {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nVtx ;  i++)
      vertices[i] = pack[i] = i ;

   for (i = hg->nVtx ; i > 0 && (*limit)>0; i--)
      {
      random = rand() % i ;
      vertex = vertices[random] ;
      vertices[random] = vertices[i-1] ;
      if (pack[vertex] != vertex)
         continue ;          /* vertex already packed, move on */

      count = hg->vindex[vertex+1] - hg->vindex[vertex] ;  /* count edges */
      if (count == 0)
         continue ;

      edge = hg->vedge[hg->vindex[vertex] + (rand() % count)] ;  /* random edge */
      for (j = hg->hindex[edge] ; j < hg->hindex[edge+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[edge+1] && (*limit)>0; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  {
                  vertex = pack[vertex] = hg->hvertex[j] ;
                  (*limit)--;
                  }
            pack[vertex] = first_vertex ;
            break ;
            }
      }
   ZOLTAN_FREE ((void **) &vertices) ;
   return ZOLTAN_OK ;
   }

/****************************************************************************/

static int grouping_rhg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
   int   i, j, *vertices=NULL, *del_edges=NULL, vertex, first_vertex, edge,
         number, best_edge, best_size, best_neighbors, random, size;
   float best_ewgt;
   char  *yo = "grouping_heg" ;

   if (!(vertices  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)) ||
       !(del_edges = (int *) ZOLTAN_CALLOC (hg->nEdge,sizeof (int))) )
      {
      ZOLTAN_FREE ((void **) &vertices) ;
      ZOLTAN_FREE ((void **) &del_edges) ;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nVtx ; i++)
      pack[i] = vertices[i] = i ;

   for (i = hg->nVtx ; i > 0 && (*limit)>0 ; i--)
      {
      number = rand() % i ;
      vertex = vertices[number] ;
      vertices[number] = vertices[i-1] ;
      if (pack[vertex] != vertex)
         continue ;            /* vertex is already matched, move on */

      best_neighbors = 0 ;
      best_edge = best_size = -1 ;
      best_ewgt = -1.0 ;
      for (j = hg->vindex[vertex] ; j < hg->vindex[vertex+1] ; j++)
         {
         int size ;
         edge = hg->vedge[j];
         size = hg->hindex[edge+1] - hg->hindex[edge] ;
         if (del_edges[edge]==0 && ((!(hg->ewgt) && size == best_size)
          || (hg->ewgt && hg->ewgt[edge] == best_ewgt && size == best_size)))
            {
            best_neighbors++;
            }
         else if (del_edges[edge]==0 && ((!(hg->ewgt) && size < best_size)
          || (hg->ewgt && (hg->ewgt[edge] >  best_ewgt
                       || (hg->ewgt[edge] == best_ewgt && size < best_size)))))
              {
              best_neighbors = 1 ;
              best_edge = edge ;
              best_ewgt = hg->ewgt[edge] ;
              best_size = hg->hindex[edge+1]-hg->hindex[edge];
              }
         }
      if (best_neighbors == 0)
         continue ;                       /* no suitable edge found */

      if (best_neighbors > 1)
         {
         random = rand() % best_neighbors;
         for (j = hg->vindex[vertex] ; j < hg->vindex[vertex+1] ; j++)
            {
            edge = hg->vedge[j] ;
            size = hg->hindex[edge+1] - hg->hindex[edge] ;
            if (del_edges[edge]==0 && ((!(hg->ewgt) && size == best_size)
              || (hg->ewgt && hg->ewgt[edge] == best_ewgt && size == best_size))
              && --best_neighbors==random)
                 {
                 best_edge = edge ;
                 break ;
                 }
            }
         }

      for (j = hg->hindex[best_edge] ; j < hg->hindex[best_edge+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[best_edge+1] && (*limit)>0 ; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  {
                  vertex = pack[vertex] = hg->hvertex[j] ;
                  (*limit)--;
                  }
            pack[vertex] = first_vertex ;
            break ;
            }
      del_edges[best_edge] = 1 ;
      }
   ZOLTAN_FREE ((void **) &vertices) ;
   ZOLTAN_FREE ((void **) &del_edges) ;
   return ZOLTAN_OK ;
}

/****************************************************************************/

static int grouping_grg (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
  int   i, j, *size=NULL, *sorted=NULL, first_vertex, vertex ;
  char *yo = "grouping_grg" ;

  for (i=0; i<hg->nVtx; i++)
    pack[i] = i;

/* Sort the hyperedges according to their weight and size */
  if (!(size   = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge))  ||
      !(sorted = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge))   )
     {
     ZOLTAN_FREE ((void **) &size) ;
     ZOLTAN_FREE ((void **) &sorted) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i=0; i<hg->nEdge; i++)
     size[i] = -(hg->hindex[i+1]-hg->hindex[i]);
  for (i=0; i<hg->nEdge; i++)
     sorted[i] = i;
  quicksort_pointer_dec_float_int(sorted,hg->ewgt,size,0,hg->nEdge-1);
  ZOLTAN_FREE ((void **) &size);

  /* Match hyperedges along decreasing weight */
  for (i=0; i<hg->nEdge && (*limit)>0; i++)
      for (j = hg->hindex[sorted[i]] ; j < hg->hindex[sorted[i]+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[sorted[i]+1] && (*limit)>0; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  {
                  vertex = pack[vertex] = hg->hvertex[j] ;
                  (*limit)--;
                  }
            pack[vertex] = first_vertex ;
            break ;
            }
  ZOLTAN_FREE ((void **) &sorted);
  return ZOLTAN_OK;
  }

/****************************************************************************/

static int grouping_aug1 (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
/* Placeholder for grouping_aug1. */
  return ZOLTAN_OK;
}

/****************************************************************************/

static int grouping_aug2 (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
/* Placeholder for grouping_aug2. */
  return ZOLTAN_OK;
}

/****************************************************************************/

static int grouping_aug3 (ZZ *zz, HGraph *hg, Packing pack, int *limit)
{
/* Placeholder for grouping_aug3. */
  return ZOLTAN_OK;
}

/****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
