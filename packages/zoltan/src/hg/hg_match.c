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
static ZOLTAN_HG_MATCHING_FN matching_aug2; /* post matching optimizer */
static ZOLTAN_HG_MATCHING_FN matching_aug3; /* post matching optimizer */

/*****************************************************************************/

int Zoltan_HG_Set_Matching_Fn(HGPartParams *hgp)
{
  static int srand_set = 0;

  if (srand_set == 0)
  { srand_set = 1 ;
    srand ((unsigned long) RANDOM_SEED) ;
  }

  if      (!strcasecmp(hgp->redm_str, "mxm"))  hgp->matching = matching_mxm;
  else if (!strcasecmp(hgp->redm_str, "rem"))  hgp->matching = matching_rem;
  else if (!strcasecmp(hgp->redm_str, "rrm"))  hgp->matching = matching_rrm;
  else if (!strcasecmp(hgp->redm_str, "rhm"))  hgp->matching = matching_rhm;
  else if (!strcasecmp(hgp->redm_str, "grm"))  hgp->matching = matching_grm;
  else if (!strcasecmp(hgp->redm_str, "lhm"))  hgp->matching = matching_lhm;
  else if (!strcasecmp(hgp->redm_str, "pgm"))  hgp->matching = matching_pgm;
  else                                         hgp->matching = NULL;

  if (hgp->matching) {

    /*
     * If reduction method is a matching, set the improvement and
     * edge weight scaling functions accordingly.
     */

    if      (!strcasecmp(hgp->rli_str, "aug2")) hgp->matching_rli = matching_aug2;
    else if (!strcasecmp(hgp->rli_str, "aug3")) hgp->matching_rli = matching_aug3;
    else                                        hgp->matching_rli = NULL;
  }
  return  hgp->matching ? 1 : 0 ;
}

/*****************************************************************************/

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
  return sim;
}

/*****************************************************************************/

int Zoltan_HG_Matching (
  ZZ *zz,
  HGraph *hg,
  Graph *g,
  Matching match,
  HGPartParams *hgp,
  int limit)
{ int ierr;
  char *yo = "Zoltan_HG_Matching";
  float *old_ewgt=NULL, *new_ewgt;

  ZOLTAN_TRACE_ENTER(zz, yo);

/*
{ int i, j;
  for (i=0; i<g->nVtx; i++)
    for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
      if (g->ewgt[j] != sim(hg,i,g->neigh[j]))
        printf("%d %d %f %f\n",i,g->neigh[j],g->ewgt[j],sim(hg,i,g->neigh[j]));
}
*/
 
  if (g->vwgt && hgp->ews)
  { if (!(new_ewgt = (float *) ZOLTAN_MALLOC (sizeof (float) * g->nEdge)))
    { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    Zoltan_HG_Scale_Graph_Weight (zz, g, new_ewgt);
    old_ewgt = g->ewgt;
    g->ewgt = new_ewgt;
  }
  
  /* Call matching routine specified by parameters. */
  ierr = hgp->matching(zz,hg,g,match,limit);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    goto End;
  }

  /* Optimization */
  if (hgp->matching_rli)
    ierr = hgp->matching_rli (zz,hg,g,match,limit);

End:
  if (g->vwgt && hgp->ews)
  { g->ewgt = old_ewgt;
    ZOLTAN_FREE ((void **) &new_ewgt);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/

static int matching_mxm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
{ int   i, j ;

  for (i=0; i<g->nVtx; i++)
    match[i] = i;

  for (i=0; i<g->nVtx && limit>0; i++)
    for (j=g->nindex[i]; j<g->nindex[i+1] && match[i]==i; j++)
      if (match[g->neigh[j]] == g->neigh[j])
      { match[i] = g->neigh[j];
        match[g->neigh[j]] = i;
        limit--;
      }
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_rem (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
   {
   int i,j,k=0, ii, jj ;
   int *v1, *v2 ;
   int random ;
   char *yo = "matching_rem" ;

   v1 = (int *) ZOLTAN_MALLOC (sizeof (int) * g->nEdge) ;
   v2 = (int *) ZOLTAN_MALLOC (sizeof (int) * g->nEdge) ;
   if (v1 == NULL || v2 == NULL)
      {
      ZOLTAN_FREE ((void **) &v1) ;
      ZOLTAN_FREE ((void **) &v2) ;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR ;
      }

   for (i = 0 ; i < g->nVtx ; i++)
      match[i] = i ;

   for (i = 0 ; i < g->nVtx ; i++)
      for (j = g->nindex[i] ; j < g->nindex[i+1] ; j++)
          if (i < g->neigh[j])
             {
             v1[k] = i ;
             v2[k] = g->neigh[j] ;
             k++ ;
             }

   for (i = k ; i > 0 && limit > 0 ; i--)
      {
      random = rand() % i ;
      ii = v1[random] ;
      jj = v2[random] ;
      v1[random] = v1[i-1] ;
      v2[random] = v2[i-1] ;

      if (match[ii] == ii && match[jj] == jj)
         {
         match[ii] = jj ;
         match[jj] = ii ;
         limit-- ;
         }
      }

   ZOLTAN_FREE (&v1) ;
   ZOLTAN_FREE (&v2) ;
   return ZOLTAN_OK ;
   }


/*****************************************************************************/

static int matching_rrm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
{ int   i, j, *vertices, vertex, number, free_neighbors, random;
  char *yo = "matching_rrm" ;

  if (!(vertices = (int *) ZOLTAN_MALLOC (sizeof (int) * g->nVtx)))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nVtx; i++)
    match[i] = vertices[i] = i;
  for (i=g->nVtx; i>0 && limit>0; i--)
  { vertex = vertices[number=rand()%i];
    vertices[number] = vertices[i-1];
    if (match[vertex] == vertex)
    { free_neighbors = 0;
      for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
        if (match[g->neigh[j]]==g->neigh[j])
          free_neighbors++;
      if (free_neighbors > 0)
      { random = rand()%free_neighbors;
        for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
          if (match[g->neigh[j]]==g->neigh[j] && --free_neighbors==random)
          { match[vertex] = g->neigh[j];
            match[g->neigh[j]] = vertex;
            limit--;
            break;
  } } }   }
  ZOLTAN_FREE ((void  **) &vertices);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_rhm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
{ int   i, j, *vertices, vertex, number, best_neighbors, random;
  float best_ewgt;
  char *yo = "matching_rhm" ;

  if (!g->ewgt)
     return matching_rrm(zz,hg,g,match,limit);

  if (!(vertices = (int *) ZOLTAN_MALLOC (sizeof (int) * g->nVtx)))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nVtx; i++)
    match[i] = vertices[i] = i;
  for (i=g->nVtx; i>0 && limit>0; i--)
  { vertex = vertices[number=rand()%i];
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
      { random = rand()%best_neighbors;
        for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
          if (match[g->neigh[j]]==g->neigh[j] && g->ewgt[j]==best_ewgt
              && --best_neighbors==random)
          { match[vertex] = g->neigh[j];
            match[g->neigh[j]] = vertex;
            limit--;
            break;
  } } }   }
  ZOLTAN_FREE ((void **) &vertices);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_grm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
{ int   i, j, *sorted, *vertex, vertex1, vertex2;
  char *yo = "matching_grm" ;

  if (!g->ewgt)
     return matching_mxm(zz,hg,g,match,limit);

  for (i=0; i<g->nVtx; i++)
    match[i] = i;
  if (!(sorted = (int *) ZOLTAN_MALLOC (sizeof (int) * g->nindex[g->nVtx])) ||
      !(vertex = (int *) ZOLTAN_MALLOC (sizeof (int) * g->nindex[g->nVtx]))  )
  { ZOLTAN_FREE ((void **) &sorted) ;
    ZOLTAN_FREE ((void **) &vertex) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nindex[g->nVtx]; i++)
    sorted[i] = i;

  for (i=0; i<g->nVtx; i++)
    for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
      vertex[j] = i;

  quicksort_pointer_dec_float (sorted,g->ewgt,0,g->nindex[g->nVtx]-1);

  for (i=0; i<g->nindex[g->nVtx] && limit>0; i++)
  { vertex1 = vertex[sorted[i]];
    vertex2 = g->neigh[sorted[i]];
    if (vertex1!=vertex2 && match[vertex1]==vertex1 && match[vertex2]==vertex2)
    { match[vertex1] = vertex2;
      match[vertex2] = vertex1;
      limit--;
  } }
  ZOLTAN_FREE ((void **) &sorted);
  ZOLTAN_FREE ((void **) &vertex);
  return ZOLTAN_OK;
}

/*****************************************************************************/

#undef LAM_ORIG

static int lhm_match (int a, int b, float ewgt_ab, int *nindex, int *neigh,
 float *ewgt, int *Nindex, int limit, int *match, int *size)
{ int	Nindex_a_old=Nindex[a], Nindex_b_old=Nindex[b], c, c_deg,
   a_deg=nindex[a+1]-nindex[a], b_deg=nindex[b+1]-nindex[b];
  float	c_ewgt;

  while (match[a]==a && match[b]==b &&
         (Nindex[a]<nindex[a+1]||Nindex[b]<nindex[b+1]) && (*size)<limit )
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
        lhm_match(a,c,c_ewgt,nindex,neigh,ewgt,Nindex,limit,match,size);
    }
    if (match[b] == b  &&  Nindex[b] < nindex[b+1])
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
	lhm_match(b,c,c_ewgt,nindex,neigh,ewgt,Nindex,limit,match,size);
  } }

  if (match[a]==a && match[b]==b && (*size)<limit)
  { match[a] = b;
    match[b] = a;
    (*size)++;
  }
  else if (match[a]==a && match[b]!=b)
    Nindex[a] = Nindex_a_old;
  else if (match[a]!=a && match[b]==b)
    Nindex[b] = Nindex_b_old;

  return ZOLTAN_OK;
}

static int matching_lhm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
{ int	i, j, size=0, *Nindex;
  char *yo = "matching_lhm" ;

  if (!g->ewgt)
     return matching_mxm(zz,hg,g,match,limit);

  for (i=0; i<g->nVtx; i++)
    match[i] = i;

  if (!(Nindex = (int *) ZOLTAN_MALLOC (sizeof (int) *  g->nVtx+1)))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  memcpy(Nindex,g->nindex,(g->nVtx+1)*sizeof(int));

  for (i=0; i<g->nVtx && size<limit; i++)
    for (j=g->nindex[i]; match[i]==i && j<g->nindex[i+1]; j++)
      if (match[g->neigh[j]] == g->neigh[j])
        lhm_match(i,g->neigh[j],g->ewgt[j],g->nindex,g->neigh,g->ewgt,Nindex,limit,match,&size);
  ZOLTAN_FREE ((void **) &Nindex);
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_pgm (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
{ int	i, j, vertex, *match1=match, *match2, *M, size1=0, size2=0,
        neighbor, next_vertex;
  float	w1=0.0, w2=0.0, weight;
  char *yo = "matching_pgm" ;

  if (!g->ewgt)
    return matching_mxm(zz,hg,g,match,limit);

  if (!(match2 = (int *) ZOLTAN_MALLOC (sizeof (int) * g->nVtx)))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<g->nVtx; i++)
    match1[i] = match2[i] = i;
  M = match1;

  for (i=0; i<g->nVtx && size2<limit; i++)
    if (match1[i]==i && match2[i]==i)
    { vertex = i;
      do
      { weight = 0.0;
        next_vertex = -1;
        for (j=g->nindex[vertex]; j<g->nindex[vertex+1]; j++)
        { neighbor = g->neigh[j];
          if (match1[neighbor]==neighbor && match2[neighbor]==neighbor && g->ewgt[j]>weight)
            { weight = g->ewgt[j];
              next_vertex = neighbor;
        }   }
        if (next_vertex >= 0)
        { M[vertex] = next_vertex;
          M[next_vertex] = vertex;
          if (M == match1)
          { size1++;
            w1 += weight;
            M = match2;
          }
          else
          { size2++;
            w2 += weight;
            M = match1;
          }
          vertex = next_vertex;
        }
      } while (next_vertex>=0 && size2<limit);
    }

  if (w1 < w2)
    for (i=0; i<g->nVtx; i++)
      match[i] = match2[i];

  ZOLTAN_FREE ((void **) &match2);
  return ZOLTAN_OK;
}

/*****************************************************************************/
static int matching_aug2 (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
{
/* Placeholder for matching_aug2. */
  return ZOLTAN_OK;
}

/*****************************************************************************/

static int matching_aug3 (ZZ *zz, HGraph *hg, Graph *g, Matching match, int limit)
{ int		i, j, k, size=0, *stack, free_p=0, best_2=-1, best_near=-1,
		best_middle=-1, best_distant=-1, neigh1, neigh2, neigh3;
  float		w_1, w_2, w_3, gain_2, gain_3;
  char *yo = "matching_aug3" ;

  if (!(stack = (int*) ZOLTAN_MALLOC (sizeof (int) * g->nVtx)))
     {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i=0; i<g->nVtx; i++)
    if (match[i] == i)
      stack[free_p++] = i;
    else
      size++;
  size /= 2;

  while (free_p && size<limit)
  { i = stack[--free_p];
    if (match[i]==i)
    { gain_2 = gain_3 = -1;
      for (j=g->nindex[i]; match[i]==i && j<g->nindex[i+1]; j++)
      { neigh1 = g->neigh[j];
        if (neigh1!=i)
        { w_1 = (g->ewgt?g->ewgt[j]:1.0);
          if (match[neigh1]==neigh1)
          { match[i] = neigh1;
  	    match[neigh1] = i;
	    size++;
          }
          else
          { neigh2 = match[neigh1];
	    for (k=g->nindex[neigh1]; neigh2!=g->neigh[k]; k++) ; /* ; is O.K.!!*/
	    w_2 = (g->ewgt?g->ewgt[k]:1.0);
	    if (w_1-w_2 > gain_2)
            { gain_2 = w_1-w_2;
              best_2 = neigh1;
            }

            for (k=g->nindex[neigh2]; k<g->nindex[neigh2+1]; k++)
            { neigh3 = g->neigh[k];
	      if (neigh3!=neigh2 && match[neigh3]==neigh3 && neigh3!=i)
	      { w_3 = (g->ewgt?g->ewgt[k]:1.0);
                if (w_1-w_2+w_3>gain_3)
	        { best_near = neigh1;
	          best_middle = neigh2;
	          best_distant = neigh3;
	          gain_3 = w_1-w_2+w_3;
      } } } } } }

      if (match[i] == i)
      { if (gain_3 >= EPS)
        { match[i] = best_near;
	  match[best_near] = i;
	  match[best_middle] = best_distant;
	  match[best_distant] = best_middle;
	  size++;
        }
        else if (gain_2>EPS)
        { match[i] = best_2;
          stack[free_p++] = match[best_2];
          match[match[best_2]] = match[best_2];
          match[best_2] = i;
  } } } }
  ZOLTAN_FREE ((void **) &stack);
  return ZOLTAN_OK;
}

/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
