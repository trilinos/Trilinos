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
static ZOLTAN_HG_MATCHING_FN matching_rhm;  /* random, heavy edge matching */
static ZOLTAN_HG_MATCHING_FN matching_grm;  /* greedy edge matching */
static ZOLTAN_HG_MATCHING_FN matching_lhm;  /* locally heaviest matching */
static ZOLTAN_HG_MATCHING_FN matching_pgm;  /* path growing matching */
static ZOLTAN_HG_MATCHING_FN matching_aug1; /* post matching optimizer */
static ZOLTAN_HG_MATCHING_FN matching_aug2; /* post matching optimizer */
static ZOLTAN_HG_MATCHING_FN matching_aug3; /* post matching optimizer */

/*****************************************************************************/

int Zoltan_HG_Set_Matching_Fn(HGPartParams *hgp)
{
  if      (!strcasecmp(hgp->redm_str, "mxm"))  hgp->matching = matching_mxm;
  else if (!strcasecmp(hgp->redm_str, "rem"))  hgp->matching = matching_rem;
  else if (!strcasecmp(hgp->redm_str, "rrm"))  hgp->matching = matching_rrm;
  else if (!strcasecmp(hgp->redm_str, "rhm"))  hgp->matching = matching_rhm;
  else if (!strcasecmp(hgp->redm_str, "grm"))  hgp->matching = matching_grm;
  else if (!strcasecmp(hgp->redm_str, "lhm"))  hgp->matching = matching_lhm;
  else if (!strcasecmp(hgp->redm_str, "pgm"))  hgp->matching = matching_pgm;
  else                                         hgp->matching = NULL;

  if (hgp->matching)
  { /* If reduction method is a matching, set the improvement and
     * edge weight scaling functions accordingly.
     */

    if      (!strcasecmp(hgp->redmo_str, "aug1")) hgp->matching_opt = matching_aug1;
    else if (!strcasecmp(hgp->redmo_str, "aug2")) hgp->matching_opt = matching_aug2;
    else if (!strcasecmp(hgp->redmo_str, "aug3")) hgp->matching_opt = matching_aug3;
    else                                          hgp->matching_opt = NULL;
  }
  return  hgp->matching ? 1 : 0 ;
}

/*****************************************************************************/

/*
static float sim (HGraph *hg, int a, int b)
{ int   i, j, edge, pins;
  float sim=0.0, weight;

  for (i=hg->vindex[a]; i<hg->vindex[a+1]; i++)
  { edge = hg->vedge[i];
    j = hg->hindex[edge];
    while (j<hg->hindex[edge+1] && hg->hvertex[j]!=b)
      j++;
    if (j < hg->hindex[edge+1])
    { pins = hg->hindex[edge+1]-hg->hindex[edge];
      weight = 2.0/((pins-1)*pins);
      if (hg->ewgt)
        weight *= hg->ewgt[edge];
      sim += weight;
    }
  }
  if (hg->vwgt)
  { if (hg->vwgt[a]<=0.0 || hg->vwgt[b]<=0.0)
      sim = FLT_MAX;
    else
      sim = sim/hg->vwgt[a]/hg->vwgt[b];
  }
  return sim;
}

static void sim_check (HGraph *hg, Graph *g)
{ int i, j;
  for (i=0; i<g->nVtx; i++)
  { for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
     if ((fabs(g->ewgt[j] - sim(hg,i,g->neigh[j]))) > EPS*(g->ewgt[j]))
       printf("%d %d %.20f %.20f\n",i,g->neigh[j],g->ewgt[j],sim(hg,i,g->neigh[j])
);
  }
}
*/

/*****************************************************************************/

int Zoltan_HG_Matching (
  ZZ *zz,
  HGraph *hg,
  Matching match,
  HGPartParams *hgp,
  int *limit)
{ int   i, ierr, need_graph ;
  char  *yo = "Zoltan_HG_Matching";
  float *old_ewgt=NULL, *new_ewgt;
  Graph g;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* temporary code while graphs are replaced by hypergraphs */
  need_graph = 0 ;
  if (hgp->matching == matching_rem
   || hgp->matching == matching_rrm
   || hgp->matching == matching_rhm
   || hgp->matching == matching_grm
   || hgp->matching == matching_lhm
   || hgp->matching == matching_pgm
   || hgp->matching_opt == matching_aug2
   || hgp->matching_opt == matching_aug3)
  {
      need_graph = 1 ;
      Zoltan_HG_Graph_Init (&g) ;
      ierr = Zoltan_HG_HGraph_to_Graph (zz, hg, &g) ;
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
         {
         Zoltan_HG_Graph_Free (&g) ;
         return ierr ;
         }
/*
    graph_connected_components (g.nVtx,g.nindex,g.neigh,zz->Debug_Level);
*/
    if (g.vwgt && hgp->ews>0)
    { if (!(new_ewgt = (float *) ZOLTAN_MALLOC (g.nEdge*sizeof(float))))
      { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      Zoltan_HG_Scale_Graph_Weight (zz, &g, new_ewgt,hgp->ews);
      old_ewgt = g.ewgt;
      g.ewgt = new_ewgt;
    }
/*
    sim_check(hg,&g);
*/
  }

  for (i=0; i<hg->nVtx; i++)
    match[i] = i;

  /* Call matching routine specified by parameters. */
  ierr = hgp->matching(zz,hg,&g,match,limit);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    goto End;
  }

  /* Optimization */
  if (hgp->matching_opt)
    ierr = hgp->matching_opt (zz,hg,&g,match,limit);

  /* compare matching weight to an upper bound */
  if (g.ewgt && zz->Debug_Level>=2)
  { int i, j;
    double highest, matching_weight=0.0, upper_bound=0.0;
    
    for (i=0; i<g.nVtx; i++)
    { highest = 0.0;
      for (j=g.nindex[i]; j<g.nindex[i+1]; j++)
      { highest = MAX(highest,g.ewgt[j]);
        if (g.neigh[j] == match[i])
          if (g.ewgt[j] == FLT_MAX)
            matching_weight = FLT_MAX;
          else
            matching_weight += (g.ewgt[j]);
      }
      if (highest < FLT_MAX)
        upper_bound += highest;
      else
        upper_bound = FLT_MAX;
    }
    printf("matching/upper_bound ratio:%.3f\n",matching_weight/upper_bound);
  }

End:
/* more temporary code while graphs are being replaced by hypergraphs */
if (need_graph == 1)
{
  if (g.vwgt && hgp->ews)
  { g.ewgt = old_ewgt;
    ZOLTAN_FREE ((void **) &new_ewgt);
  }
  Zoltan_HG_Graph_Free(&g) ;
}

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/

static int matching_mxm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int  i, j, edge, vertex, *Hindex;
  char *yo = "matching_mxm";

  if (!(Hindex    = (int *) ZOLTAN_MALLOC (sizeof (int) * (hg->nEdge+1)) ))
  { ZOLTAN_FREE ((void **) &Hindex) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR ;
  }
  for (i=0; i<=hg->nEdge; i++)
    Hindex[i] = hg->hindex[i];
  
  for (i=0; i<hg->nVtx && (*limit)>0; i++)
    for (j=hg->vindex[i]; match[i]==i && j<hg->vindex[i+1]; j++)
    { edge = hg->vedge[j];
      while (Hindex[edge]<hg->hindex[edge+1])
      { vertex = hg->hvertex[Hindex[edge]++];
        if (vertex!=i && match[vertex]==vertex)
        { match[i] = vertex;
          match[vertex] = i;
          (*limit)--;
          break ;   /* terminate while loop */
    } } }

  ZOLTAN_FREE ((void **)&Hindex);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_rem (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int i, j, k=0, *v1, *v2, random ;
  char *yo = "matching_rem" ;

  if (!(v1    = (int *) ZOLTAN_MALLOC (g->nEdge*sizeof(int))))
  { ZOLTAN_FREE ((void **) &v1) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR ;
  }
  v2 = &(v1[g->nEdge/2]);
  for (i = 0 ; i < g->nVtx ; i++)
    for (j = g->nindex[i] ; j < g->nindex[i+1] ; j++)
      if (i < g->neigh[j])
      { v1[k] = i ;
        v2[k++] = g->neigh[j] ;
      }
  while (k>0 && (*limit)>0)
  { i = v1[random=Zoltan_HG_Rand()%k];
    j = v2[random];
    v1[random] = v1[--k] ;
    v2[random] = v2[k] ;
    if (match[i] == i && match[j] == j)
    { match[i] = j ;
      match[j] = i ;
      (*limit)-- ;
    }
  }

  ZOLTAN_FREE ((void **) &v1) ;
  return ZOLTAN_OK ;
}


/*****************************************************************************/

static int matching_rrm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int   i, j, *vertices, vertex, number, free_neighbors, random;
  char *yo = "matching_rrm" ;

  if (!(vertices = (int *) ZOLTAN_MALLOC (g->nVtx*sizeof(int))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nVtx; i++)
    vertices[i] = i;
  for (i=g->nVtx; i>0 && (*limit)>0; i--)
  { vertex = vertices[number=Zoltan_HG_Rand()%i];
    vertices[number] = vertices[i-1];
    if (match[vertex] == vertex)
    { free_neighbors = 0;
      for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
        if (match[g->neigh[j]]==g->neigh[j])
          free_neighbors++;
      if (free_neighbors > 0)
      { random = Zoltan_HG_Rand()%free_neighbors;
        for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
          if (match[g->neigh[j]]==g->neigh[j] && --free_neighbors==random)
          { match[vertex] = g->neigh[j];
            match[g->neigh[j]] = vertex;
            (*limit)--;
            break;
  } } }   }
  ZOLTAN_FREE ((void  **) &vertices);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_rhm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int   i, j, *vertices, vertex, number, best_neighbors, random;
  float best_ewgt;
  char *yo = "matching_rhm" ;

  if (!g->ewgt)
     return matching_rrm(zz,hg,g,match,limit);

  if (!(vertices = (int *) ZOLTAN_MALLOC (g->nVtx*sizeof(int))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nVtx; i++)
    vertices[i] = i;
  for (i=g->nVtx; i>0 && (*limit)>0; i--)
  { vertex = vertices[number=Zoltan_HG_Rand()%i];
    vertices[number] = vertices[i-1];
    if (match[vertex] == vertex)
    { best_neighbors = 0;
      best_ewgt = 0.0;
      for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
        if (match[g->neigh[j]]==g->neigh[j])
        { if (g->ewgt[j] > best_ewgt)
          { best_neighbors = 1;
            best_ewgt = g->ewgt[j];
          }
          else if (g->ewgt[j] == best_ewgt)
            best_neighbors++;
        }
      if (best_neighbors > 0)
      { random = Zoltan_HG_Rand()%best_neighbors;
        for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
          if (match[g->neigh[j]]==g->neigh[j] && g->ewgt[j]==best_ewgt
              && --best_neighbors==random)
          { match[vertex] = g->neigh[j];
            match[g->neigh[j]] = vertex;
            (*limit)--;
            break;
  } } }   }
  ZOLTAN_FREE ((void **) &vertices);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_grm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int   i, j, *order, *vertex, vertex1, vertex2;
  char *yo = "matching_grm" ;

  if (!g->ewgt)
    return matching_mxm(zz,hg,g,match,limit);

  if (!(order  = (int *) ZOLTAN_MALLOC (g->nEdge*sizeof(int))) ||
      !(vertex = (int *) ZOLTAN_MALLOC (g->nEdge*sizeof(int)))  )
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_FREE ((void **) &vertex) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nEdge; i++)
    order[i] = i;

  for (i=0; i<g->nVtx; i++)
    for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
      vertex[j] = i;

  quicksort_pointer_dec_float (order,g->ewgt,0,g->nEdge-1);

  for (i=0; i<g->nEdge && (*limit)>0; i++)
  { vertex1 = vertex[order[i]];
    vertex2 = g->neigh[order[i]];
    if (vertex1!=vertex2 && match[vertex1]==vertex1 && match[vertex2]==vertex2)
    { match[vertex1] = vertex2;
      match[vertex2] = vertex1;
      (*limit)--;
  } }
  ZOLTAN_FREE ((void **) &order);
  ZOLTAN_FREE ((void **) &vertex);
  return ZOLTAN_OK;
}

/*****************************************************************************/

#undef LAM_ORIG

static int lhm_match (int a, int b, float ewgt_ab, int *nindex, int *neigh,
 float *ewgt, int *Nindex, int *limit, int *match)
{ int	Nindex_a_old=Nindex[a], Nindex_b_old=Nindex[b], c, c_deg,
   a_deg=nindex[a+1]-nindex[a], b_deg=nindex[b+1]-nindex[b];
  float	c_ewgt;

  while (match[a]==a && match[b]==b &&
         (Nindex[a]<nindex[a+1]||Nindex[b]<nindex[b+1]) && (*limit)>0 )
  { if (Nindex[a] < nindex[a+1])
    { c = neigh[Nindex[a]];
      c_ewgt = ewgt[Nindex[a]];
      c_deg = nindex[c+1]-nindex[c];
      (Nindex[a])++;
#ifdef LAM_ORIG
      if (c_ewgt>ewgt_ab)
#else
      if (( b_deg>1 && ((c_deg>1&&(c_ewgt>ewgt_ab||(c_ewgt==ewgt_ab&&c_deg<b_deg)))
       || (c_deg==1&&2.0*c_ewgt>=ewgt_ab))) ||
	  ( b_deg==1 && ((c_deg>1&&c_ewgt>2.0*ewgt_ab)||(c_deg==1&&c_ewgt>ewgt_ab)) ) )
#endif
        lhm_match(a,c,c_ewgt,nindex,neigh,ewgt,Nindex,limit,match);
    }
    if (match[b] == b  &&  Nindex[b] < nindex[b+1]  && (*limit)>0)
    { c = neigh[Nindex[b]];
      c_ewgt = ewgt[Nindex[b]];
      c_deg = nindex[c+1]-nindex[c];
      (Nindex[b])++;
#ifdef LAM_ORIG
      if (c_ewgt>ewgt_ab)
#else
      if (( a_deg>1 && ((c_deg>1&&(c_ewgt>ewgt_ab||(c_ewgt==ewgt_ab&&c_deg<a_deg)))
       || (c_deg==1&&2.0*c_ewgt>=ewgt_ab))) ||
	  ( a_deg==1 && ((c_deg>1&&c_ewgt>2.0*ewgt_ab)||(c_deg==1&&c_ewgt>ewgt_ab)) ))
#endif
	lhm_match(b,c,c_ewgt,nindex,neigh,ewgt,Nindex,limit,match);
  } }

  if (match[a]==a && match[b]==b && (*limit)>0)
  { match[a] = b;
    match[b] = a;
    (*limit)--;
  }
  else if (match[a]==a && match[b]!=b)
    Nindex[a] = Nindex_a_old;
  else if (match[a]!=a && match[b]==b)
    Nindex[b] = Nindex_b_old;

  return ZOLTAN_OK;
}

static int matching_lhm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int	i, j, *Nindex;
  char *yo = "matching_lhm" ;

  if (!g->ewgt)
     return matching_mxm(zz,hg,g,match,limit);

  if (!(Nindex = (int *) ZOLTAN_MALLOC ((g->nVtx+1)*sizeof(int))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  memcpy(Nindex,g->nindex,(g->nVtx+1)*sizeof(int));

  for (i=0; i<g->nVtx && (*limit)>0; i++)
    for (j=g->nindex[i]; match[i]==i && j<g->nindex[i+1]; j++)
      if (match[g->neigh[j]] == g->neigh[j])
        lhm_match(i,g->neigh[j],g->ewgt[j],g->nindex,g->neigh,g->ewgt,Nindex,limit,match);
  ZOLTAN_FREE ((void **) &Nindex);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_pgm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int	i, j, side=0, vertex, *Match[2], limits[2], neighbor, next_vertex;
  float	w[2]={0.0,0.0}, weight;
  char *yo = "matching_pgm" ;

  if (!g->ewgt)
    return matching_mxm(zz,hg,g,match,limit);

  limits[0] = limits[1] = (*limit);
  Match[0] = match;
  if (!(Match[1] = (int *) ZOLTAN_MALLOC (g->nVtx*sizeof(int))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nVtx; i++)
    Match[1][i] = i;

  for (i=0; i<g->nVtx && limits[side]>0; i++)
    if (Match[0][i]==i && Match[1][i]==i)
    { vertex = i;
      while (vertex>0 && limits[side]>0)
      { weight = 0.0;
        next_vertex = -1;
        for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
        { neighbor = g->neigh[j];
          if (Match[0][neighbor]==neighbor && Match[1][neighbor]==neighbor && g->ewgt[j]>weight)
            { weight = g->ewgt[j];
              next_vertex = neighbor;
        }   }
        if (next_vertex >= 0)
        { Match[side][vertex] = next_vertex;
          Match[side][next_vertex] = vertex;
          limits[side]--;
          w[side] += weight;
          side = 1-side;
        }
        vertex = next_vertex;
      }
    }

  if (w[0] < w[1])
  { for (i=0; i<g->nVtx; i++)
      match[i] = Match[1][i];
    (*limit) = limits[1];
  }
  else
    (*limit) = limits[0];

  ZOLTAN_FREE ((void **) &Match[1]);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_aug1 (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ 
  return matching_mxm (zz,hg,g,match,limit);
}

/*****************************************************************************/

static int matching_aug2 (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int    i, j, k, *stack, free_p=0, best_2=-1, neigh1, neigh2;
  float  gain_2;
  char   *yo = "matching_aug2" ;

  if (!(g->ewgt))
    return matching_aug1 (zz,hg,g,match,limit);

  if (!(stack = (int*) ZOLTAN_MALLOC (g->nVtx*sizeof(int))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nVtx; i++)
    if (match[i] == i)
      stack[free_p++] = i;

  while (free_p && (*limit)>0)
  { i = stack[--free_p];
    if (match[i]==i)
    { gain_2 = 0.0;
      for (j=g->nindex[i]; match[i]==i && j<g->nindex[i+1]; j++)
      { neigh1 = g->neigh[j];
        if (match[neigh1]==neigh1)
        { match[i] = neigh1;
          match[neigh1] = i;
          (*limit)--;
        }
        else if (g->ewgt[j] > gain_2)
        { neigh2 = match[neigh1];
          k=g->nindex[neigh1];
          while (neigh2 != g->neigh[k])
            k++;
          if (g->ewgt[j]-g->ewgt[k] > gain_2)
          { gain_2 = g->ewgt[j]-g->ewgt[k];
            best_2 = neigh1;
      } } }

      if (match[i]==i && gain_2>0.0)
      { match[i] = best_2;
        stack[free_p++] = match[best_2];
        match[match[best_2]] = match[best_2];
        match[best_2] = i;
  } } }

  ZOLTAN_FREE ((void **) &stack);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_aug3 (ZZ *zz, HGraph *hg, Graph *g, Matching match, int *limit)
{ int   i, j, k, *stack, free_p=0, best_2=-1, best_near=-1,
        best_middle=-1, best_distant=-1, neigh1, neigh2, neigh3;
  float gain, gain_2, gain_3;
  char  *yo = "matching_aug3" ;

  if (!(stack = (int*) ZOLTAN_MALLOC (g->nVtx*sizeof(int))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nVtx; i++)
    if (match[i] == i)
      stack[free_p++] = i;

  while (free_p && (*limit)>0)
  { i = stack[--free_p];
    if (match[i]==i)
    { gain_2 = 0.0;
      gain_3 = -FLT_MIN ;
      for (j=g->nindex[i]; match[i]==i && j<g->nindex[i+1]; j++)
      { neigh1 = g->neigh[j];
        if (neigh1!=i)
        { if (match[neigh1]==neigh1)
          { match[i] = neigh1;
  	    match[neigh1] = i;
            (*limit)--;
          }
          else
          { neigh2 = match[neigh1];
            if (g->ewgt)
            { k=g->nindex[neigh1];
              while (neigh2 != g->neigh[k])
                k++;
              gain = g->ewgt[j]-g->ewgt[k];
              if (gain > gain_2)
              { gain_2 = gain;
                best_2 = neigh1;
              }
              for (k=g->nindex[neigh2]; k<g->nindex[neigh2+1]; k++)
              { neigh3 = g->neigh[k];
	        if (match[neigh3]==neigh3 && neigh3!=i && gain+g->ewgt[k]>gain_3)
                { best_near = neigh1;
                  best_middle = neigh2;
                  best_distant = neigh3;
                  gain_3 = gain+g->ewgt[k];
            } } }
            else
            { for (k=g->nindex[neigh2]; k<g->nindex[neigh2+1] && gain_3<1.0; k++)
              { neigh3 = g->neigh[k];
                if (match[neigh3]==neigh3 && neigh3!=i)
                { best_near = neigh1;
                  best_middle = neigh2;
                  best_distant = neigh3;
                  gain_3 = 1.0;
      } } } } } }
      if (match[i] == i)
      { if (gain_3 > -FLT_MIN)
        { match[i] = best_near;
	  match[best_near] = i;
	  match[best_middle] = best_distant;
	  match[best_distant] = best_middle;
          (*limit)--;
        }
        else if (gain_2 > 0.0)
        { match[i] = best_2;
          stack[free_p++] = match[best_2];
          match[match[best_2]] = match[best_2];
          match[best_2] = i;
  } } } }

  ZOLTAN_FREE ((void **) &stack);
  return ZOLTAN_OK;
}

/*****************************************************************************/

/*
int graph_connected_components (int n, int *ep, int *edge, int Out)
{ int   i=0, *locked, *queue, queue_head, queue_tail,
        components=0, first_unlocked=0;

  if (!(queue  = (int*) ZOLTAN_MALLOC (n*sizeof(int))) ||
      !(locked = (int*) ZOLTAN_CALLOC (n,sizeof(int)))   )
    return ZOLTAN_MEMERR;

  if (Out > 1)
    printf ("Component Sizes     : ");
  while (i<n)
  { int k, size=1, vertex;
    queue_head=0;
    queue_tail=1;
    components++;
    while (locked[first_unlocked])
      first_unlocked++;
    queue[0] = first_unlocked;
    locked[first_unlocked] = 1;
    while (queue_head < queue_tail)
    { vertex = queue[queue_head++];
      for (k=ep[vertex]; k<ep[vertex+1]; k++)
        if (!(locked[edge[k]]))
        { size++;
          locked[edge[k]] = 1;
          queue[queue_tail++] = edge[k];
    }   }
    if (Out > 1)
      printf ("%d ", size);
    i += size;
  }
  if (Out > 1)
    printf ("\n");
  if (Out > 0)
    printf("Components          : %d\n",components);
  ZOLTAN_FREE((void **) &locked);
  ZOLTAN_FREE((void **) &queue);
  return components;
}
*/

/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
