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

static ZOLTAN_HG_LOCAL_REF_FN local_no;
static ZOLTAN_HG_LOCAL_REF_FN local_fm;

/****************************************************************************/

ZOLTAN_HG_LOCAL_REF_FN *Zoltan_HG_Set_Local_Ref_Fn(char *str)
{
  if      (!strcasecmp(str, "fm")) return local_fm;
  else if (!strcasecmp(str, "no")) return local_no;
  else                             return NULL;
}

/****************************************************************************/

int Zoltan_HG_Local (ZZ *zz, HGraph *hg, int p, Partition part, HGPartParams *hgp)
{
  return hgp->local_ref(zz,hg,p,part,hgp->bal_tol);
}

/****************************************************************************/

static int local_no (
  ZZ *zz,
  HGraph *hg,
  int p,
  Partition part,
  float bal_tol
)
{ return ZOLTAN_OK;
}

/****************************************************************************/

/*
static int gain_check (HGraph *hg, float *gain, int *part, int **cut)
{ float g;
  int vertex, j, edge;

  for (vertex=0; vertex<hg->nVtx; vertex++)
  { g = 0.0;
    for (j=hg->vindex[vertex]; j<hg->vindex[vertex+1]; j++)
    { edge = hg->vedge[j];
      if (cut[part[vertex]][edge] == 1)
        g += (hg->ewgt?(hg->ewgt[edge]):1.0);
      else if (cut[1-part[vertex]][edge] == 0)
        g -= (hg->ewgt?(hg->ewgt[edge]):1.0);
    }
    if (g != gain[vertex]){
      printf("Wrong gain %f %f\n",g,gain[vertex]);
      return ZOLTAN_FATAL;
    }
  }
  return ZOLTAN_OK;
}
*/

/****************************************************************************/

int move_vertex (HGraph *hg, int vertex, int sour, int dest, int *part,
	int **cut, float *gain, HEAP *heap)
{ int i, j, edge, v;

  gain[vertex] = 0.0;
  part[vertex] = dest;
  
  for (i=hg->vindex[vertex]; i<hg->vindex[vertex+1]; i++)
  { edge = hg->vedge[i];
    if (cut[sour][edge] == 1)
    { for (j=hg->hindex[edge]; j<hg->hindex[edge+1]; j++)
      { v = hg->hvertex[j];
        gain[v] -= (hg->ewgt?hg->ewgt[edge]:1.0);
        if (heap) heap_change_value(&(heap[part[v]]),v,gain[v]);
    } }
    else if (cut[sour][edge] == 2)
    { for (j=hg->hindex[edge]; j<hg->hindex[edge+1]; j++)
      { v = hg->hvertex[j];
        if (part[v] == sour)
        { gain[v] += (hg->ewgt?hg->ewgt[edge]:1.0);
          if (heap) heap_change_value(&(heap[part[v]]),v,gain[v]);
          break;
    } } }
    if (cut[dest][edge] == 0)
    { for (j=hg->hindex[edge]; j<hg->hindex[edge+1]; j++)
      { v = hg->hvertex[j];
        gain[v] += (hg->ewgt?hg->ewgt[edge]:1.0);
        if (heap) heap_change_value(&(heap[part[v]]),v,gain[v]);
    } }
    else if (cut[dest][edge] == 1)
    { for (j=hg->hindex[edge]; j<hg->hindex[edge+1]; j++)
      { v = hg->hvertex[j];
        if (v!=vertex && part[v]==dest)
        { gain[v] -= (hg->ewgt?hg->ewgt[edge]:1.0);
          if (heap) heap_change_value(&(heap[part[v]]),v,gain[v]);
          break;
    } } }
    cut[sour][edge]--;
    cut[dest][edge]++;
  }
  return ZOLTAN_OK;
}

/****************************************************************************/

static int local_fm (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part,
  float bal_tol
)
{ int    i, j, vertex, edge, *cut[2], *locked, *locked_list, round=0;
  float  total_weight, max_weight, max, best_max_weight, *gain,
         part_weight[2], cutsize, best_cutsize;
  HEAP   heap[2];
  char   *yo="local_fm";

  if (p != 2)
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "p!=2 not yet implemented for local_fm.");
    return ZOLTAN_FATAL;
  }

  if (hg->nEdge == 0)
    return ZOLTAN_OK;

  /* Calculate the different weights */
  for (i=0; i<p; i++)
    part_weight[i] = 0.0;
  if (hg->vwgt)
  { total_weight = 0.0;
    for (i=0; i<hg->nVtx; i++)
    { total_weight += hg->vwgt[i];
      part_weight[part[i]] += hg->vwgt[i];
    }
  }
  else
  { total_weight = (float)(hg->nVtx);
    for (i=0; i<hg->nVtx; i++)
      part_weight[part[i]] += 1.0;

  }
  max_weight = (total_weight/(float)p)*bal_tol;
  
  if (!(cut[0]      = (int *)  ZOLTAN_CALLOC(2*hg->nEdge,sizeof(int))) ||
      !(locked      = (int *)  ZOLTAN_CALLOC(hg->nVtx,sizeof(int)))    ||
      !(locked_list = (int *)  ZOLTAN_CALLOC(hg->nVtx,sizeof(int)))    ||
      !(gain        = (float *)ZOLTAN_CALLOC(hg->nVtx,sizeof(float)))   )
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  cut[1] = &(cut[0][hg->nEdge]);

  /* Initial calculation of the cut distribution and gain values */
  for (i=0; i<hg->nEdge; i++)
    for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      (cut[part[hg->hvertex[j]]][i])++;
  for (i=0; i<hg->nVtx; i++)
    for (j=hg->vindex[i]; j<hg->vindex[i+1]; j++)
    { edge = hg->vedge[j];
      if (cut[part[i]][edge] == 1)
        gain[i] += (hg->ewgt?(hg->ewgt[edge]):1.0);
      else if (cut[1-part[i]][edge] == 0)
        gain[i] -= (hg->ewgt?(hg->ewgt[edge]):1.0);
    }

  /* Initialize the heaps and fill them with the gain values */
  for (i=0; i<p; i++)
    heap_init(zz, &(heap[i]),hg->nVtx);
  for (i=0; i<hg->nVtx; i++)
    heap_input(&(heap[part[i]]),i,gain[i]);
  for (i=0; i<p; i++)
    heap_make(&(heap[i]));

  /* Initialize given partition as best partition */
  best_cutsize = cutsize = hcut_size_total(hg,part);
  best_max_weight = part_weight[0];
  for (i=1; i<p; i++)
    best_max_weight = MAX(best_max_weight,part_weight[i]);

  do
  { int   step=0, no_better_steps=0, number_locked=0, best_locked=0,
          sour, dest;
    float akt_cutsize=best_cutsize;

    round++;
    cutsize = best_cutsize;
    if (zz->Debug_Level > ZOLTAN_DEBUG_LIST)
      printf("ROUND %d:\nSTEP VERTEX  PARTS MAX_WGT CHANGE CUTSIZE\n",round);

    while (step<hg->nVtx && no_better_steps<hg->nVtx/4)
    { step++;
      no_better_steps ++;

      if (heap_empty(&(heap[0])))
        sour = 1;
      else if (heap_empty(&(heap[1])))
        sour = 0;
      else if (part_weight[0] > max_weight)
        sour = 0;
      else if (part_weight[1] > max_weight)
        sour = 1;
      else if (heap_max_value(&(heap[0])) > heap_max_value(&(heap[1])))
        sour = 0;
      else
        sour = 1;
      dest = 1-sour;

      vertex = heap_extract_max(&(heap[sour]));
      locked[vertex] = part[vertex]+1;
      locked_list[number_locked++] = vertex;
      akt_cutsize -= gain[vertex];

      move_vertex (hg,vertex,sour,dest,part,cut,gain,heap);
      part_weight[sour] -= (hg->vwgt?hg->vwgt[vertex]:1.0);
      part_weight[dest] += (hg->vwgt?hg->vwgt[vertex]:1.0);

      max = MAX(part_weight[0],part_weight[1]);
      if ((best_max_weight>max_weight && max<best_max_weight) ||
          (max<=max_weight && akt_cutsize<best_cutsize))
      { best_locked = number_locked;
        best_cutsize = akt_cutsize;
        best_max_weight = max;
        if (zz->Debug_Level > ZOLTAN_DEBUG_LIST)
          printf ("New Partition:%f\n",akt_cutsize);
        no_better_steps = 0;
      }
      if (zz->Debug_Level > ZOLTAN_DEBUG_LIST+1)
        printf ("%4d %6d %2d->%2d %7.2f %f %f\n",step,vertex,sour,dest,
                max,akt_cutsize-cutsize,akt_cutsize);
    }

    while (number_locked != best_locked)
    { vertex = locked_list[--number_locked];
      sour = part[vertex];
      dest = locked[vertex]-1;
      move_vertex (hg,vertex,sour,dest,part,cut,gain,heap);
      part_weight[sour] -= (hg->vwgt?hg->vwgt[vertex]:1.0);
      part_weight[dest] += (hg->vwgt?hg->vwgt[vertex]:1.0);
      heap_input(&(heap[dest]),vertex,gain[vertex]);
      locked[vertex] = 0;
    }

    while (number_locked)
    { vertex = locked_list[--number_locked];
      locked[vertex] = 0;
      heap_input(&(heap[part[vertex]]),vertex,gain[vertex]);
    }
   
    for (i=0; i<p; i++)
      heap_make(&(heap[i]));
  } while (best_cutsize < cutsize);

/*
  gain_check (hg, gain, part, cut);
*/

  ZOLTAN_FREE ((void **) &(cut[0]));
  ZOLTAN_FREE ((void **) &(locked));
  ZOLTAN_FREE ((void **) &(locked_list));
  ZOLTAN_FREE ((void **) &(gain));
  heap_free(&(heap[0]));
  heap_free(&(heap[1]));

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
