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



static int *HSIZE;
static int *SUM;
static int compare (const int *a, const int *b)
{ 
  if (SUM[(*a)] < SUM[(*b)])
    return -1;
  else if (SUM[(*a)] > SUM[(*b)])
    return 1;
  else
  { if (HSIZE[(*a)] < HSIZE[(*b)])
      return -1;
    else if (HSIZE[(*a)] > HSIZE[(*b)])
      return 1;
  }
  return 0;
}


static int compare_int (const int *a, const int *b)
{ if ((*a) < (*b))
    return -1;
  else if ((*a) > (*b))
    return 1;
  return 0;
}


static int *HINDEX;
static int *HVERTEX;
static int compare_HSIZE (const int *a, const int *b)
{ int k=0;

  for (k=0; k<HSIZE[(*a)]; k++)
  { if (HVERTEX[HINDEX[(*a)]+k] < HVERTEX[HINDEX[(*b)]+k])
      return -1;
    else if (HVERTEX[HINDEX[(*a)]+k] > HVERTEX[HINDEX[(*b)]+k])
      return 1;
  }
  return 0;
}



/* Time O(|P|*log(|P|)), due to sorting */
int Zoltan_HG_Coarsening (
  ZZ *zz,
  HGraph *hg, 
  int *pack, 
  HGraph *c_hg, 
  int *LevelMap)
{ 
  int	i, j, k, l, old, vertex, new_vertex, *used_vertices=NULL,
   *sorted=NULL, deleted_he, deleted_pins, *c_hindex=NULL, *c_hvertex=NULL;
  float *c_ewgt=NULL;
  char *yo = "hgraph_coarsening" ;

  Zoltan_HG_HGraph_Init(c_hg);
  c_hg->info = hg->info+1;

/* Calculate the number of vertices */
  c_hg->nVtx = 0;
  for (i=0; i<hg->nVtx; i++)
    if (pack[i] >= 0)
    { (c_hg->nVtx)++;
      vertex = i;
      while (pack[vertex] >= 0)
      { old = vertex;
        vertex = pack[old];
        pack[old] = -pack[old]-1;
    } }
  c_hg->vwgt = (float*) ZOLTAN_MALLOC (sizeof (float) * c_hg->nVtx) ;
  if (c_hg->vwgt == NULL)
     {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i = 0 ; i < c_hg->nVtx ; i++)
      c_hg->vwgt[i] = 0.0 ;

/* Construct the LevelMap */
  new_vertex = 0;
  for (i=0; i<hg->nVtx; i++)
    if (pack[i] < 0)
    { vertex = i;
      while (pack[vertex] < 0)
      { LevelMap[vertex] = new_vertex;
        c_hg->vwgt[new_vertex] += (hg->vwgt?hg->vwgt[vertex]:1.0);
        pack[vertex] = -pack[vertex]-1;
        vertex = pack[vertex];
      }
      new_vertex++;
    }

/* Coarsen the hyperedges */
  used_vertices = (int *)   ZOLTAN_MALLOC (sizeof (int)   * c_hg->nVtx) ;
  c_ewgt        = (float *) ZOLTAN_MALLOC (sizeof (float) *   hg->nEdge) ;
  c_hindex      = (int *)   ZOLTAN_MALLOC (sizeof (int)   *  (hg->nEdge+1)) ;
  c_hvertex     = (int *)   ZOLTAN_MALLOC (sizeof (int)   *   hg->nPin) ;
  if (used_vertices==NULL || c_ewgt==NULL || c_hindex==NULL || c_hvertex==NULL)
     {
     ZOLTAN_FREE ((void **) &used_vertices) ;
     ZOLTAN_FREE ((void **) &c_hindex) ;
     ZOLTAN_FREE ((void **) &c_hvertex) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i = 0 ; i < c_hg->nVtx ; i++)
     used_vertices[i] = 0 ;

  c_hindex[0] = c_hg->nEdge = c_hg->nPin = 0;
  for (i=0; i<hg->nEdge; i++)
  { for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
    { new_vertex = LevelMap[hg->hvertex[j]];
      if (used_vertices[new_vertex] <= i)
      { used_vertices[new_vertex] = i+1;
        c_hvertex[(c_hg->nPin)++] = new_vertex;
    } }
    if (c_hg->nPin > c_hindex[c_hg->nEdge]+1)
    { c_ewgt[c_hg->nEdge++] = (hg->ewgt?hg->ewgt[i]:1.0);
      c_hindex[c_hg->nEdge] = c_hg->nPin;
    }
    else
      c_hg->nPin = c_hindex[c_hg->nEdge];
  }
  ZOLTAN_FREE((void **) &used_vertices);

/* Move weight of identical hyperedges to one of them */
  sorted = (int *) ZOLTAN_MALLOC (sizeof (int) * c_hg->nEdge) ;
  HSIZE  = (int *) ZOLTAN_MALLOC (sizeof (int) * c_hg->nEdge) ;
  SUM    = (int *) ZOLTAN_MALLOC (sizeof (int) * c_hg->nEdge) ;
  if (sorted == NULL || HSIZE == NULL || SUM == NULL)
     {
     ZOLTAN_FREE ((void **) &c_hindex) ;
     ZOLTAN_FREE ((void **) &c_hvertex) ;
     ZOLTAN_FREE ((void **) &sorted) ;
     ZOLTAN_FREE ((void **) &HSIZE) ;
     ZOLTAN_FREE ((void **) &SUM) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i = 0; i < c_hg->nEdge; i++)
  { sorted[i] = i;
    HSIZE[i] = c_hindex[i+1]-c_hindex[i];
    SUM[i] = 0 ;
    for (j=c_hindex[i]; j<c_hindex[i+1]; j++)
      SUM[i] += c_hvertex[j];
  }

/* sort the edges according to their size and their sum of vertexIDs */
  qsort(sorted,c_hg->nEdge,sizeof(int),(int(*)(const void*,const void*))compare);

  HINDEX = c_hindex;
  HVERTEX = c_hvertex;
  i = 0;
  while (i<c_hg->nEdge)
  { j = i+1;
    while (j<c_hg->nEdge && SUM[sorted[j]]==SUM[sorted[i]]
     && HSIZE[sorted[j]]==HSIZE[sorted[i]])
       j++;
    if (j > i+1)
    { for (k=i; k<j; k++)
        /* sort the vertex list of these edges */
        qsort (&(c_hvertex[c_hindex[sorted[k]]]),
         c_hindex[sorted[k]+1] -c_hindex[sorted[k]], sizeof(int),
        (int(*)(const void*, const void*))compare_int);

      /* sorted edges according to their sorted vertex lists */
      qsort (&(sorted[i]), j-i, sizeof(int),
       (int(*)(const void*, const void*))compare_HSIZE);
      /* check if the vertex lists of two hyperedges are identical */
      for (k=i; k<j-1; k++)
      { for (l=0; l<HSIZE[sorted[i]]
         && c_hvertex[c_hindex[sorted[k]]+l]==c_hvertex[c_hindex[sorted[k+1]]+l];
         l++) ;
        if (l==HSIZE[sorted[i]])
        { c_ewgt[sorted[k+1]] += c_ewgt[sorted[k]];
          c_ewgt[sorted[k]] = 0.0;
      } }
      i = j;
    }
    else
      i++;
  }
  ZOLTAN_FREE((void **) &sorted);
  ZOLTAN_FREE((void **) &HSIZE);
  ZOLTAN_FREE((void **) &SUM);

/* delete hyperedges without weight */
  deleted_he = deleted_pins = 0;
  for (i=0; i<c_hg->nEdge; i++)
  { if (c_ewgt[i] <= EPS)
    { deleted_he ++;
      deleted_pins += (c_hindex[i+1]-c_hindex[i]);
    }
    else if (deleted_pins > 0)
    { for (j=c_hindex[i]; j<c_hindex[i+1]; j++)
        c_hvertex[j-deleted_pins] = c_hvertex[j];
      c_ewgt[i-deleted_he] = c_ewgt[i];
      c_hindex[i+1-deleted_he] = c_hindex[i+1]-deleted_pins;
  } }
  c_hg->nEdge -= deleted_he;
  c_hg->nPin = c_hindex[c_hg->nEdge];

/* Reallocate the arrays to their exact size */
  c_hg->ewgt = (float *) ZOLTAN_MALLOC (sizeof (float) * c_hg->nEdge);
  if (c_hg->ewgt == NULL)
     {
     ZOLTAN_FREE ((void **) &c_hindex) ;
     ZOLTAN_FREE ((void **) &c_hvertex) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  memcpy(c_hg->ewgt,c_ewgt,c_hg->nEdge*sizeof(float));
  ZOLTAN_FREE ((void **) &c_ewgt);

  c_hg->hindex = (int *) ZOLTAN_MALLOC (sizeof (int) * (c_hg->nEdge+1));
  if (c_hg->hindex == NULL)
     {
     ZOLTAN_FREE ((void **) &c_hindex) ;
     ZOLTAN_FREE ((void **) &c_hvertex) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  memcpy(c_hg->hindex,c_hindex,(c_hg->nEdge+1)*sizeof(int));
  ZOLTAN_FREE ((void **) &c_hindex);

  c_hg->hvertex = (int *) ZOLTAN_MALLOC (sizeof (int) * c_hg->nPin);
  if (c_hg->hvertex == NULL)
     {
     ZOLTAN_FREE ((void **) &c_hvertex) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  memcpy(c_hg->hvertex,c_hvertex,c_hg->nPin*sizeof(int));
  ZOLTAN_FREE ((void **) &c_hvertex);
  return Zoltan_HG_Create_Mirror(zz,c_hg);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
