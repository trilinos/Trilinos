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

/****************************************************************************/

/* Procedure to coarsen a hypergraph based on a packing. All vertices
   of one pack are clustered to one vertex. Identical hyperedges are
   collapsed to a single hyperedge with combined weight. The array
   LevelMap is the mapping of the old vertices to the new vertices. It
   will be used to pass a partition of the coarse graph back to the
   original graph. Time O(|I|*log(|I|)), due to sorting.
*/
int Zoltan_HG_Coarsening (
  ZZ *zz,
  HGraph *hg, 
  int *pack, 
  HGraph *c_hg, 
  int *LevelMap)
{ 
  int  i, j, k, l, old, vertex, new_vertex, *hsize, *sum, *used_vertices,
       *sorted, deleted_he, deleted_pins, *c_hindex, *c_hvertex;
  float *c_ewgt;
  char *yo = "Zoltan_HG_Coarsening" ;

  Zoltan_HG_HGraph_Init(c_hg);
  c_hg->info = hg->info+1;

/* Calculate the number of coarse vertices. Change pack[vertex] to
   -pack[vertex]-1.
*/
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
  if (!(c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx,sizeof(float))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

/* Construct the LevelMap. pack[vertex] is changed back to original value.
   Also sum up the vertex weights.
*/
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

/* Coarsen the hyperedges and avoid edges with only one vertex */
  if (!(used_vertices = (int *)  ZOLTAN_CALLOC(c_hg->nVtx,sizeof(int)))   ||
      !(c_ewgt        = (float *)ZOLTAN_MALLOC(hg->nEdge*sizeof(float)))  ||
      !(c_hindex      = (int *)  ZOLTAN_MALLOC((hg->nEdge+1)*sizeof(int)))||
      !(c_hvertex     = (int *)  ZOLTAN_MALLOC(hg->nInput*sizeof(int)))      )
  { ZOLTAN_FREE ((void **) &used_vertices) ;
    ZOLTAN_FREE ((void **) &c_hindex) ;
    ZOLTAN_FREE ((void **) &c_hvertex) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  c_hindex[0] = c_hindex[1] = c_hg->nEdge = c_hg->nInput = 0;
  for (i=0; i<hg->nEdge; i++)
  { for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
    { new_vertex = LevelMap[hg->hvertex[j]];
      if (used_vertices[new_vertex] <= i)
      { used_vertices[new_vertex] = i+1;
        c_hvertex[(c_hg->nInput)++] = new_vertex;
    } }
    if (c_hg->nInput > c_hindex[c_hg->nEdge]+1)   /* a new hyperedge */
    { c_ewgt[c_hg->nEdge++] = (hg->ewgt?hg->ewgt[i]:1.0);
      c_hindex[c_hg->nEdge] = c_hg->nInput;
    }
    else /* no new hyperedge, because it only covers one vertex */
      c_hg->nInput = c_hindex[c_hg->nEdge];
  }
  ZOLTAN_FREE((void **) &used_vertices);

/* Done if there are no remaining edges */
  if (c_hg->nEdge == 0)
  { c_hg->ewgt = NULL;
    if (!(c_hg->hindex = (int *) ZOLTAN_CALLOC (1,sizeof(int))))
    { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }
    c_hg->hvertex = NULL;
    ZOLTAN_FREE ((void **) &c_ewgt);
    ZOLTAN_FREE ((void **) &c_hindex) ;
    ZOLTAN_FREE ((void **) &c_hvertex) ;
    return Zoltan_HG_Create_Mirror(zz,c_hg);
  }

/* Move weight of identical hyperedges to one of them */
  if (!(sorted = (int *) ZOLTAN_MALLOC (c_hg->nEdge*sizeof(int)))  ||
      !(hsize  = (int *) ZOLTAN_MALLOC (c_hg->nEdge*sizeof (int))) ||
      !(sum    = (int *) ZOLTAN_CALLOC (c_hg->nEdge,sizeof (int)))  )
  { ZOLTAN_FREE ((void **) &c_hindex) ;
    ZOLTAN_FREE ((void **) &c_hvertex) ;
    ZOLTAN_FREE ((void **) &sorted) ;
    ZOLTAN_FREE ((void **) &hsize) ;
    ZOLTAN_FREE ((void **) &sum) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i = 0; i < c_hg->nEdge; i++)
  { sorted[i] = i;
    hsize[i] = c_hindex[i+1]-c_hindex[i];
    for (j=c_hindex[i]; j<c_hindex[i+1]; j++)
      sum[i] += c_hvertex[j];
  }
  /* sort the edges according to their size and their sum of vertexIDs */
  quicksort_pointer_inc_int_int (sorted, sum, hsize, 0, c_hg->nEdge-1);
  i = 0;
  while (i<c_hg->nEdge)
  { j = i+1;
    while (j<c_hg->nEdge && sum[sorted[j]]==sum[sorted[i]]
     && hsize[sorted[j]]==hsize[sorted[i]])
       j++;
    if (j > i+1)
    { for (k=i; k<j; k++)
        /* sort the vertex list of these edges */
        quicksort_list_inc_int (&(c_hvertex[c_hindex[sorted[k]]]),0,c_hindex[sorted[k]+1]-c_hindex[sorted[k]]-1);
      /* sort edges according to their sorted vertex lists */
      quicksort_pointer_inc_int_mult (sorted,i,j-1,c_hindex,c_hvertex);
      /* check if the vertex lists of two hyperedges are identical */
      for (k=i; k<j-1; k++)
      { l = 0;
        while (l<hsize[sorted[i]] && c_hvertex[c_hindex[sorted[k]]+l]==c_hvertex[c_hindex[sorted[k+1]]+l])
          l++;
        if (l==hsize[sorted[i]])
        { c_ewgt[sorted[k+1]] += c_ewgt[sorted[k]];
          c_ewgt[sorted[k]] = 0.0;
      } }
      i = j;
    }
    else
      i++;
  }
  ZOLTAN_FREE((void **) &sorted);
  ZOLTAN_FREE((void **) &hsize);
  ZOLTAN_FREE((void **) &sum);

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
  c_hg->nInput = c_hindex[c_hg->nEdge];

/* Reallocate the arrays to their exact size */
  if (!(c_hg->ewgt = (float *) ZOLTAN_MALLOC (c_hg->nEdge*sizeof(float))))
  { ZOLTAN_FREE ((void **) &c_hindex) ;
    ZOLTAN_FREE ((void **) &c_hvertex) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  memcpy(c_hg->ewgt,c_ewgt,c_hg->nEdge*sizeof(float));
  ZOLTAN_FREE ((void **) &c_ewgt);
  if (!(c_hg->hindex = (int *) ZOLTAN_MALLOC ((c_hg->nEdge+1)*sizeof(int))))
  { ZOLTAN_FREE ((void **) &c_hindex) ;
    ZOLTAN_FREE ((void **) &c_hvertex) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  memcpy(c_hg->hindex,c_hindex,(c_hg->nEdge+1)*sizeof(int));
  ZOLTAN_FREE ((void **) &c_hindex);
  if (!(c_hg->hvertex = (int *) ZOLTAN_MALLOC (c_hg->nInput*sizeof(int))))
  { ZOLTAN_FREE ((void **) &c_hvertex) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  memcpy(c_hg->hvertex,c_hvertex,c_hg->nInput*sizeof(int));
  ZOLTAN_FREE ((void **) &c_hvertex);

  return Zoltan_HG_Create_Mirror(zz,c_hg);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
