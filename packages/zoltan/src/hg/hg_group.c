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

extern int srand_set;   /* flag indicating random number generator is seeded */

static ZOLTAN_HG_GROUPING_FN grouping_mxg;  /* maximal packing */
static ZOLTAN_HG_GROUPING_FN grouping_reg;  /* random edge packing */
static ZOLTAN_HG_GROUPING_FN grouping_heg;  /* heavy edge packing */
static ZOLTAN_HG_GROUPING_FN grouping_grp;  /* greedy packing */
/* static ZOLTAN_HG_GROUPING_FN grouping_lhp; */ /* locally heaviest packing */
/* static ZOLTAN_HG_GROUPING_FN grouping_pgp; */ /* path growing packing */
static ZOLTAN_HG_GROUPING_FN grouping_deg;  /* decreasing packing */
static ZOLTAN_HG_GROUPING_FN grouping_rrg;  /* random vertex, random edge packing */

/****************************************************************************/

ZOLTAN_HG_GROUPING_FN *Zoltan_HG_Set_Grouping_Fn(char *str)
{

  if      (strcasecmp(str, "mxg") == 0)  return grouping_mxg;
  else if (strcasecmp(str, "reg") == 0)  return grouping_reg;
  else if (strcasecmp(str, "heg") == 0)  return grouping_heg;
  else if (strcasecmp(str, "grp") == 0)  return grouping_grp;
  else if (strcasecmp(str, "deg") == 0)  return grouping_deg;
  else if (strcasecmp(str, "rrg") == 0)  return grouping_rrg;

/*
  else if (strcasecmp(str, "lhp") == 0)  return grouping_lhp;
  else if (strcasecmp(str, "pgp") == 0)  return grouping_pgp;
*/

  else                                   return NULL;
}

/****************************************************************************/

int Zoltan_HG_Grouping (ZZ *zz, HGraph *hg, Packing pack, HGParams *hgp)
{ int   i, j;
  int   ierr = ZOLTAN_OK;
  float *old_ewgt=NULL, weight, sum1, sum2;
  char  *yo = "Zoltan_HG_Grouping";

  old_ewgt = hg->ewgt;
  hg->ewgt = (float *) ZOLTAN_MALLOC (sizeof (float) * hg->nEdge);
  if (hg->ewgt == NULL)
     {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }

  for (i=0; i<hg->nEdge; i++)
  { sum1 = sum2 = 0.0;
    if (hg->vwgt)
    { for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      { weight = hg->vwgt[hg->hvertex[j]];
        sum1 += weight;
        sum2 += weight*weight;
    } }
    else
      sum1 = sum2 = (float)(hg->hindex[i+1]-hg->hindex[i]);
    sum1 = (sum1*sum1-sum2)/2.0;
    if (sum1 == 0.0)
      hg->ewgt[i] = FLT_MAX;
    else
      hg->ewgt[i] = (old_ewgt?old_ewgt[i]:1.0)/sum1;
  }

  ierr = hgp->packing(zz,hg,pack);
  ZOLTAN_FREE ((void **) &hg->ewgt);
  hg->ewgt = old_ewgt;
  return ierr;
}

/****************************************************************************/

static int grouping_mxg (ZZ *zz, HGraph *hg, Packing pack)
   {
   int i, j, vertex, first_vertex ;
   char *yo = "grouping_mxg" ;

   for (i = 0 ; i < hg->nVtx ; i++)
      pack[i] = i ;

   for (i = 0 ; i < hg->nEdge ; i++)
      for (j = hg->hindex[i] ; j < hg->hindex[i+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[i+1] ; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  vertex = pack[vertex] = hg->hvertex[j] ;
            pack[vertex] = first_vertex ;
            break ;       /* not required, might improve speed */
            }
   return ZOLTAN_OK ;
   }

/****************************************************************************/

static int grouping_reg (ZZ *zz, HGraph *hg, Packing pack)
   {
   int i, j, *edges=NULL, edge, random, vertex, first_vertex ;
   char *yo = "grouping_reg" ;

   if (!srand_set)
      {
      srand_set = 1 ;
      srand ((unsigned long) RANDOM_SEED) ;
      }

   edges = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge) ;
   if (edges == NULL)
      {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nEdge ; i++)
      edges[i] = i ;
   for (i = 0 ; i < hg->nVtx ;  i++)
      pack[i]  = i ;

   for (i = hg->nEdge ; i > 0 ; i--)
      {
      random = rand() % i ;
      edge = edges[random] ;
      edges[random] = edges[i-1] ;

      for (j = hg->hindex[edge] ; j < hg->hindex[edge+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[edge+1] ; j++)
                if (pack[hg->hvertex[j]] == hg->hvertex[j])
                   vertex = pack[vertex] = hg->hvertex[j] ;
            pack[vertex] = first_vertex ;
            break ;
            }
      }
   ZOLTAN_FREE ((void **) &edges) ;
   return ZOLTAN_OK ;
   }

/****************************************************************************/

static int grouping_rrg (ZZ *zz, HGraph *hg, Packing pack)
   {
   int i, j, k, edge, random, *vertices=NULL, vertex ;
   int *del_edges=NULL, count ;
   char *yo = "grouping_rrg" ;

   if (!srand_set)
      {
      srand_set = 1 ;
      srand ((unsigned long) RANDOM_SEED) ;
      }

   vertices  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx) ;
   del_edges = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge) ;
   if (vertices == NULL || del_edges == NULL)
      {
      ZOLTAN_FREE ((void **) &vertices) ;
      ZOLTAN_FREE ((void **) &del_edges) ;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nVtx ;  i++)
      vertices[i] = pack[i] = i ;
   for (i = 0 ; i < hg->nEdge ; i++)
      del_edges[i] = 0 ;

   for (i = hg->nVtx ; i > 0 ; i--)
      {
      random = rand() % i ;
      vertex = vertices[random] ;
      vertices[random] = vertices[i-1] ;
      if (pack[vertex] != vertex)
         continue ;          /* vertex already packed, move on */

      count = 0 ;            /* count will be number of viable edges */
      for (k = hg->vindex[vertex] ; k < hg->vindex[vertex+1] ; k++)
         {
         edge = hg->vedge[k] ;
         if (del_edges[edge] == 1)
            continue ;       /* edge has been deleted for use already */

         for (j = hg->hindex[edge] ; j < hg->hindex[edge+1] ; j++)
            if (pack[hg->hvertex[j]] != hg->hvertex[j])
               break ;
         if (j == hg->hindex[edge+1])
            count++ ;                 /* found a possible edge */
         else
            del_edges[edge] = 1 ;     /* edge unusable, delete it from consideration */
         }
      if (count == 0)
         continue ;                  /* vertex has no free edges available */

      random = (count > 1) ? rand() % count : 0 ;      /* randomly select from available edges */
      for (k = hg->vindex[vertex] ; k < hg->vindex[vertex+1] ; k++)
         {
         edge = hg->vedge[k] ;
         if (del_edges[edge] == 0 && --count == random)
            {
            for (j = hg->hindex[edge] ; j < hg->hindex[edge+1]-1 ; j++)
               pack[hg->hvertex[j]] = hg->hvertex[j+1] ;
            pack[hg->hvertex[j]] = hg->hvertex[hg->hindex[edge]] ;
            break ;              /* packed edge, escape to outer loop */
            }
         del_edges[edge] = 1 ;
         }
      }
   ZOLTAN_FREE ((void **) &vertices) ;
   ZOLTAN_FREE ((void **) &del_edges) ;
   return ZOLTAN_OK ;
   }

/****************************************************************************/

static int grouping_heg (ZZ *zz, HGraph *hg, Packing pack)
{
   int   i, j, k, *vertices=NULL, *del_edges=NULL, vertex, edge,
         number, best_edge, best_size;
   float best_ewgt;
   char  *yo = "grouping_heg" ;

   if (!srand_set)
     {
     srand_set = 1 ;
     srand ((unsigned long) RANDOM_SEED) ;
     }

   vertices  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx) ;
   del_edges = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge) ;
   if (vertices == NULL  ||  del_edges == NULL)
      {
      ZOLTAN_FREE ((void **) &vertices) ;
      ZOLTAN_FREE ((void **) &del_edges) ;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nVtx ; i++)
      pack[i] = vertices[i] = i ;
   for (i=0; i<hg->nEdge; i++)
      del_edges[i] = 0;

   for (i = hg->nVtx ; i > 0 ; i--)
      {
      number = rand() % i ;
      vertex = vertices[number] ;
      vertices[number] = vertices[i-1] ;
      if (pack[vertex] != vertex)
         continue ;            /* vertex is already matched, move on */

      best_edge = best_size = -1;
      best_ewgt = -1.0 ;
      for (j = hg->vindex[vertex] ; j < hg->vindex[vertex+1] ; j++)
      { edge = hg->vedge[j];
        if (del_edges[edge]==0 &&
             (hg->ewgt && (hg->ewgt[edge]>best_ewgt || (hg->ewgt[edge]==best_ewgt && hg->hindex[edge+1]-hg->hindex[edge]<best_size))) ||
             (hg->ewgt==NULL && hg->hindex[edge+1]-hg->hindex[edge]<best_size))
            {
            best_edge = edge ;
            best_ewgt = hg->ewgt[best_edge] ;
            best_size = hg->hindex[best_edge+1]-hg->hindex[best_edge];
            }
      }
      if (best_edge == -1)
         continue ;

      for (j = hg->hindex[best_edge] ; j < hg->hindex[best_edge+1] ; j++)
        for (k=hg->vindex[hg->hvertex[j]]; k<hg->vindex[hg->hvertex[j]+1]; k++)
          del_edges[hg->vedge[k]] = 1;

      for (j = hg->hindex[best_edge] ; j < hg->hindex[best_edge+1]-1 ; j++)
         pack[hg->hvertex[j]] = hg->hvertex[j+1] ;
      pack[hg->hvertex[j]] = hg->hvertex[hg->hindex[best_edge]] ;
      }
   ZOLTAN_FREE ((void **) &vertices) ;
   ZOLTAN_FREE ((void **) &del_edges) ;
   return ZOLTAN_OK ;
}

/****************************************************************************/

static void quickpart_dec_float_int (float *val1, int *val2, int *sorted,
int start, int end, int *equal, int *smaller)
{ int	i, next, key2, key2_next;
  float key1, key1_next;

  i = (end+start)/2;
  key1 = val1?val1[sorted[i]]:1.0;
  key2 = val2?val2[sorted[i]]:1;

  (*equal) = (*smaller) = start;
  for (i=start; i<=end; i++)
  { next = sorted[i];
    key1_next = val1?val1[next]:1.0;
    key2_next = val2?val2[next]:1;
    if (key1_next>key1 || (key1_next==key1 && key2_next>key2))
    { sorted[i] = sorted[(*smaller)];
      sorted[(*smaller)++] = sorted[(*equal)];
      sorted[(*equal)++] = next;
    }
    else if (key1_next==key1 && key2_next==key2)
    { sorted[i] = sorted[(*smaller)];
      sorted[(*smaller)++] = next;
  } }
}


static void quicksort_dec_float_int (float* val1, int *val2, int* sorted,
int start, int end)
{ int  equal, smaller;

  if (start < end)
  { quickpart_dec_float_int (val1,val2,sorted,start,end,&equal,&smaller);
    quicksort_dec_float_int (val1,val2,sorted,start,equal-1);
    quicksort_dec_float_int (val1,val2,sorted,smaller,end);
  }
}

/****************************************************************************/

static int grouping_grp (ZZ *zz, HGraph *hg, Packing pack)
{ int   i, j, *size=NULL, *sorted=NULL;
  char *yo = "grouping_grp" ;

  for (i=0; i<hg->nVtx; i++)
    pack[i] = i;

/* Sort the hyperedges according to their weight and size */
  size   = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  sorted = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  if (size == NULL || sorted == NULL)
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
  quicksort_dec_float_int(hg->ewgt,size,sorted,0,hg->nEdge-1);
  ZOLTAN_FREE ((void **) &size);

/* Match hyperedges along decreasing weight */
  for (i=0; i<hg->nEdge; i++)
  { for (j=hg->hindex[sorted[i]]; j<hg->hindex[sorted[i]+1]; j++)
      if (pack[hg->hvertex[j]] != hg->hvertex[j])
        break;
    if (j == hg->hindex[sorted[i]+1])
    { for (j=hg->hindex[sorted[i]]; j<hg->hindex[sorted[i]+1]-1; j++)
        pack[hg->hvertex[j]] = hg->hvertex[j+1];
      pack[hg->hvertex[j]] = hg->hvertex[hg->hindex[sorted[i]]];
  } }
  ZOLTAN_FREE ((void **) &sorted);
  return ZOLTAN_OK;
}

/****************************************************************************/

static int lhp_pack (ZZ *zz, HGraph *hg, int edge, int *del_edge, int *Vindex,
                     Packing pack)
{ int	i, j, vertex, *Vindex_old, next_edge, done=0;
  char *yo = "lhp_pack" ;

  Vindex_old = (int *) ZOLTAN_MALLOC (sizeof (int)
   * (hg->hindex[edge+1]-hg->hindex[edge]));
  if (Vindex_old == NULL)
     {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i=0; i<hg->hindex[edge+1]-hg->hindex[edge]; i++)
    Vindex_old[i] = Vindex[hg->hvertex[i+hg->hindex[edge]]];

  while (hg->ewgt && del_edge[edge]==0 && done==0)
  { done = 1;
    for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++)
    { vertex = hg->hvertex[i];
      if (Vindex[vertex] < hg->vindex[vertex+1])
      { next_edge = hg->vedge[Vindex[vertex]];
        (Vindex[vertex])++;
        if (hg->ewgt[next_edge] > hg->ewgt[edge])
          lhp_pack(zz,hg,next_edge,del_edge,Vindex,pack);
        done=0;
  } } }

  if (del_edge[edge] == 0)
  { for (i=hg->hindex[edge]; i<hg->hindex[edge+1]-1; i++)
      pack[hg->hvertex[i]] = hg->hvertex[i+1];
    pack[hg->hvertex[i]] = hg->hvertex[hg->hindex[edge]];
    for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++)
    { vertex = hg->hvertex[i];
      for (j=hg->vindex[vertex]; j<hg->vindex[vertex+1]; j++)
        del_edge[hg->vedge[j]] = 1;
    }
  }
  else
    for (i=0; i<hg->hindex[edge+1]-hg->hindex[edge]; i++)
    { vertex = hg->hvertex[i+hg->hindex[edge]];
      if (pack[vertex] == vertex)
        Vindex[vertex] = Vindex_old[i];
    }

  ZOLTAN_FREE((void **) &Vindex_old);
  return ZOLTAN_OK;
}

/****************************************************************************/

static int grouping_lhp (ZZ *zz, HGraph *hg, Packing pack)
{ int	i, *Vindex=NULL, *del_edge=NULL;
  char *yo = "grouping_lhp" ;

  for (i=0; i<hg->nVtx; i++)
    pack[i] = i;

  del_edge = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  Vindex   = (int *) ZOLTAN_MALLOC (sizeof (int) * (hg->nVtx+1));
  if (del_edge == NULL || Vindex == NULL)
     {
     ZOLTAN_FREE ((void **) &del_edge) ;
     ZOLTAN_FREE ((void **) &Vindex) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i=0; i<hg->nEdge; i++)
     del_edge[i] = 0 ;
  memcpy(Vindex,hg->vindex,(hg->nVtx+1)*sizeof(int));

  for (i=0; i<hg->nEdge; i++)
    if (del_edge[i] == 0)
      lhp_pack(zz,hg,i,del_edge,Vindex,pack);

  ZOLTAN_FREE ((void **) &del_edge);
  ZOLTAN_FREE ((void **) &Vindex);
  return ZOLTAN_OK;
}

/****************************************************************************/

static int grouping_pgp (ZZ *zz, HGraph *hg, Packing pack)
{ int   i, j, k, vertex, edge, *pack1=pack, *pack2=NULL, *Pack=pack,
        *taken_edge=NULL, *taken_vertex=NULL, cur_edge, best_edge;
  float	best_weight, w1=0.0, w2=0.0;
  char *yo = "grouping_pgp" ;

  pack2        = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx);
  taken_edge   = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  taken_vertex = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx);
  if (pack2 == NULL || taken_edge == NULL || taken_vertex == NULL)
     {
     ZOLTAN_FREE ((void **) &pack2) ;
     ZOLTAN_FREE ((void **) &taken_edge) ;
     ZOLTAN_FREE ((void **) &taken_vertex) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i=0; i<hg->nEdge; i++)
     taken_edge[i]=0 ;
  for (i=0; i< hg->nVtx; i++)
     taken_vertex[i]=0 ;
  for (i=0; i<hg->nVtx; i++)
     pack1[i] = pack2[i] = i;

  for (i=0; i<hg->nEdge; i++)
    if (taken_edge[i] == 0)
    { cur_edge = i;
      taken_edge[cur_edge] = 1;
      while (cur_edge >= 0)
      { for (j=hg->hindex[cur_edge]; j<hg->hindex[cur_edge+1]-1; j++)
          Pack[hg->hvertex[j]] = hg->hvertex[j+1];
        Pack[hg->hvertex[j]] = hg->hvertex[hg->hindex[cur_edge]];
        if (Pack == pack1)
        { w1 += (hg->ewgt?hg->ewgt[cur_edge]:1.0);
          Pack = pack2;
        }
        else
        { w2 += (hg->ewgt?hg->ewgt[cur_edge]:1.0);
          Pack = pack1;
        }

        best_weight = 0.0;
        best_edge = -1;
        for (j=hg->hindex[cur_edge]; j<hg->hindex[cur_edge+1]; j++)
          if (taken_vertex[(vertex=hg->hvertex[j])] == 0)
          { for (k=hg->vindex[vertex]; k<hg->vindex[vertex+1]; k++)
              if (taken_edge[(edge=hg->vedge[k])]==0)
              { if ((hg->ewgt?hg->ewgt[edge]:1.0)>best_weight)
                { best_edge = edge;
                  best_weight = (hg->ewgt?hg->ewgt[best_edge]:1.0);
                }
                taken_edge[edge] = 1;
              }
            taken_vertex[vertex] = 1;
          }
        cur_edge = best_edge;
    } }

  if (w1 < w2)
    memcpy(pack,pack2,(hg->nVtx) * sizeof(int));

  ZOLTAN_FREE ((void **) &pack2);
  ZOLTAN_FREE ((void **) &taken_edge);
  ZOLTAN_FREE ((void **) &taken_vertex);
  return ZOLTAN_OK;
}

/****************************************************************************/

static int grouping_deg (ZZ *zz, HGraph *hg, Packing pack)
{ int   i, j, *size=NULL, *sorted=NULL, vertex, first_vertex;
  char *yo = "grouping_deg" ;

  for (i=0; i<hg->nVtx; i++)
    pack[i] = i;

/* Sort the hyperedges according to their weight and size */
  size   = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  sorted = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  if (size == NULL || sorted == NULL)
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
  quicksort_dec_float_int(hg->ewgt,size,sorted,0,hg->nEdge-1);
  ZOLTAN_FREE ((void **) &size);

/* Pack hyperedges along decreasing weight */
  for (i=0; i<hg->nEdge; i++)
  { j = hg->hindex[sorted[i]];
    while (j<hg->hindex[sorted[i]+1] && pack[hg->hvertex[j]]!=hg->hvertex[j])
      j++;
    if (j<hg->hindex[sorted[i]+1])
    { first_vertex = vertex = hg->hvertex[j];
      j++;
      while (j<hg->hindex[sorted[i]+1])
      { if (pack[hg->hvertex[j]] == hg->hvertex[j])
        { pack[vertex] = hg->hvertex[j];
          vertex = hg->hvertex[j];
        }
        j++;
      }
      pack[vertex] = first_vertex;
  } }
  ZOLTAN_FREE ((void **) &sorted);
  return ZOLTAN_OK;
}

/****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
