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
static ZOLTAN_HG_LOCAL_REF_FN local_hc;

/****************************************************************************/

ZOLTAN_HG_LOCAL_REF_FN *Zoltan_HG_Set_Local_Ref_Fn(char *str)
{
  if      (!strcasecmp(str, "hc")) return local_hc;
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

static int local_hc (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part,
  float bal_tol
)
{ int    i, j, side, vertex, edge, *cut[2], *boundary[2], boundary_n[2],
         best_vertex, *locked, round=0;
  float  total_weight, max_weight, improvement, best_improvement,
         part_weight[2], cutsize, best_cutsize;
  char   *yo="local_hc";

  if (p != 2)
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "p!=2 not yet implemented for local_hc.");
    return ZOLTAN_FATAL;
  }

  if (hg->nEdge == 0)
    return ZOLTAN_OK;
/*
  printf("START cuts:%.2f %.2f\n",hcut_size_total(hg,part),hcut_size_links (zz,hg,2,part));
*/

  best_cutsize = cutsize = hcut_size_total(hg,part);

  if (hg->vwgt)
  { total_weight = 0.0;
    for (i=0; i<hg->nVtx; i++)
      total_weight += hg->vwgt[i];
  }
  else
    total_weight = (float)(hg->nVtx);
  max_weight = (total_weight/(float)p)*bal_tol;
/*
  printf("total_weight:%f max_weight:%f \n",total_weight, max_weight);
*/

  part_weight[0] = part_weight[1] = 0.0;
  for (i=0; i<hg->nVtx; i++)
    part_weight[part[i]] += (hg->vwgt?hg->vwgt[i]:1.0);
/*
  printf("weights: %f %f\n",part_weight[0],part_weight[1]);
*/

  if (!(boundary[0] = (int *) ZOLTAN_MALLOC(2*hg->nVtx* sizeof(int))) || 
      !(cut[0]      = (int *) ZOLTAN_CALLOC(2*hg->nEdge,sizeof(int))) ||
      !(locked      = (int *) ZOLTAN_CALLOC(hg->nVtx,sizeof(int)))     )
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  boundary[1] = &(boundary[0][hg->nVtx]);
  cut[1] = &(cut[0][hg->nEdge]);

  for (i=0; i<hg->nEdge; i++)
  { cut[0][i] = cut[1][i] = 0;
    for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      (cut[part[hg->hvertex[j]]][i])++;
  }

  do
  { boundary_n[0] = boundary_n[1] = 0;
    for (i=0; i<hg->nVtx; i++)
    { for (j=hg->vindex[i]; j<hg->vindex[i+1]; j++)
        if (cut[0][hg->vedge[j]]==1 || cut[1][hg->vedge[j]]==1)
        { boundary[part[i]][boundary_n[part[i]]++] = i;
          break;
        }
    }
/*
    printf("boundaries:%d %d\n",boundary_n[0],boundary_n[1]);
*/
    best_vertex = -1;
    best_improvement = 0.0;

    for (side=0; side<2; side++)
    { for (i=0; i<boundary_n[side]; i++)
      { vertex = boundary[side][i];
        if (part_weight[1-part[vertex]]+(hg->vwgt?hg->vwgt[vertex]:1.0) < max_weight)
        { improvement = 0.0;
          for (j=hg->vindex[vertex]; j<hg->vindex[vertex+1]; j++)
          { edge = hg->vedge[j];
            if (cut[part[vertex]][edge] == 1)
              improvement += (hg->ewgt?(hg->ewgt[edge]):1.0);
            else if (cut[1-part[vertex]][edge] == 0)
              improvement -= (hg->ewgt?(hg->ewgt[edge]):1.0); 
          }

          if (improvement > best_improvement)
          { best_vertex = vertex;
            best_improvement = improvement;
          }
        }
      }
    }
/*
    printf("best: %d %f %d\n",best_vertex, best_improvement,part[best_vertex]);
*/
    if (best_improvement > 0.0)
    { if (best_vertex < 0)
        return ZOLTAN_FATAL;

      for (i=hg->vindex[best_vertex]; i<hg->vindex[best_vertex+1]; i++)
      { cut[part[best_vertex]][hg->vedge[i]] --;
        cut[1-part[best_vertex]][hg->vedge[i]] ++;
      }
      part_weight[part[best_vertex]] -= (hg->vwgt?hg->vwgt[best_vertex]:1.0);
      part[best_vertex] = 1-part[best_vertex];
      part_weight[part[best_vertex]] += (hg->vwgt?hg->vwgt[best_vertex]:1.0);
    }
  } while (best_improvement > 0.0);

  ZOLTAN_FREE ((void **) &(boundary[0]));
  ZOLTAN_FREE ((void **) &(cut[0]));
  ZOLTAN_FREE ((void **) &(locked));

/*
  printf("END   cuts:%.2f %.2f\n",hcut_size_total(hg,part),hcut_size_links (zz,hg,2,part));
*/
  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif


/*
#define gain(V) \
{ int _j; for(_j=0;_j<p;_j++) edges[_j]=0; \
  orient[(V)]=(part[(V)]+1)%p; \
  for(_j=ep[(V)];_j<ep[(V)+1];_j++) \
    edges[part[edge[_j]]] += (ew?ew[_j]:1); \
  in[(V)] = edges[part[(V)]]; \
  for(_j=0;_j<p;_j++) \
    if (_j!=part[(V)] && edges[_j]>edges[orient[(V)]]) \
      orient[(V)] = _j; \
  ex[(V)]=edges[orient[(V)]]; \
}

int local_kl ()
{ int		vertex, neig, *in, *ex, *orient, *locked_list, *edges;
  int		max_vw=1, max_w, best_max_w=0;
  BUCKETS	**B;
  BUCKET_ELE	*bucket_ele;

  CALLOC("LOCAL_KL",in,int,n);
  CALLOC("LOCAL_KL",ex,int,n);
  CALLOC("LOCAL_KL",orient,int,n);
  CALLOC("LOCAL_KL",locked_list,int,n);
  CALLOC("LOCAL_KL",edges,int,p);
  CALLOC("LOCAL_KL",B,BUCKETS*,p);
  CALLOC("LOCAL_KL",bucket_ele,BUCKET_ELE,n);

  for (i=0; i<n; i++)
    bucket_ele[i].element = (void*)i;
  for (i=0; i<p; i++)
  { if (bucket_cal(&(B[i]),-maxdegree, maxdegree))
      ERROR ("LOCAL_KL", "bucket_cal");
    best_max_w = MAX(best_max_w,part_w[i]);
  }
  for (i=0; i<n; i++)
  { gain(i);
    if (ex[i]>0)
      bucket_in(B[part[i]],&(bucket_ele[i]),ex[i]-in[i]);
  }

  do 
  { int	step=0, akt_cutsize=best_cutsize, number_locked=0, 
	best_locked=0, sour, dest, vertex_j;
    int	no_better_steps=0, vweight, max_sour, max_j;
    
    round++;
    cutsize = best_cutsize;
    if (Out > 3)
      printf("ROUND %d:\nSTEP VERTEX  PARTS MAX_WGT CHANGE CUTSIZE\n",round);

    while (step<total_w && no_better_steps<total_w/4)
    { step++;
      for(sour=0; sour<p && bucket_empty(B[sour]); sour++);
      if (sour == p)
	break;
      for (j=sour+1; j<p; j++)
        if (bucket_not_empty(B[j]))
	{ if (part_w[j] >= max_part_w)
	  { if (part_w[sour] >= max_part_w)
	    { vertex = (int)(bucket_max_element(B[sour]));
	      vweight = vw[vertex]);
	      max_sour = MAX(part_w[sour]-vweight,part_w[orient[vertex]]+vweight);
	      vertex_j = (int)(bucket_max_element(B[j]));
	      vweight = vw[vertex_j]);
	      max_j = MAX(part_w[j]-vweight,part_w[orient[vertex_j]]+vweight);
              if (max_j<max_sour || (max_j==max_sour && ex[vertex_j]-in[vertex_j]>ex[vertex]-in[vertex]))
	        sour = j;
            }
	    else
	      sour = j;
	  }
	  else
	  { if (part_w[sour] < max_part_w)
	    { vertex = (int)(bucket_max_element(B[sour]));
	      vweight = (vw?vw[vertex]:1);
	      max_sour = MAX(part_w[sour]-vweight,part_w[orient[vertex]]+vweight);
	      vertex_j = (int)(bucket_max_element(B[j]));
	      vweight = vw[vertex_j]);
	      max_j = MAX(part_w[j]-vweight,part_w[orient[vertex_j]]+vweight);
	      if (max_j >= max_part_w)
	      { if (max_sour>=max_part_w && (max_j<max_sour || (max_j==max_sour && ex[vertex_j]-in[vertex_j]>ex[vertex]-in[vertex])))
		  sour = j;
	      }
	      else if (max_sour>=max_part_w || (ex[vertex_j]-in[vertex_j]>ex[vertex]-in[vertex] || (ex[vertex_j]-in[vertex_j]==ex[vertex]-in[vertex] && max_j<max_sour)))
		sour = j;
	    }
	  }
	}

      vertex = (int)(bucket_max_element(B[sour]));
      bucket_del (&(bucket_ele[vertex]));
      locked[vertex] = part[vertex]+1;
      locked_list[number_locked++] = vertex;
      dest = part[vertex] = orient[vertex];
      akt_cutsize -= (ex[vertex]-in[vertex]);
      gain(vertex);
      part_w[sour] -= vw[vertex];
      part_w[dest] += vw[vertex];
      no_better_steps += vw[vertex];

      for (j=ep[vertex]; j<ep[vertex+1]; j++)
	if (!locked[neig=edge[j]])
        { gain(neig);
	  if (buckets(bucket_ele[neig]))
	  { if (ex[neig]==0)
	      bucket_del (&(bucket_ele[neig]));
            else
	      bucket_new_key(&(bucket_ele[neig]), ex[neig]-in[neig]);
          }
	  else if (ex[neig]>0)
	    bucket_in(B[part[neig]],&(bucket_ele[neig]),ex[neig]-in[neig]);
        }

      max_w = 0;
      for (i=0; i<p; i++)
	max_w = MAX(max_w,part_w[i]);
      if ((best_max_w>max_part_w && max_w<best_max_w) ||
	  (max_w<=max_part_w && akt_cutsize<best_cutsize))
      { best_locked = number_locked;
        best_cutsize = akt_cutsize;
	best_max_w = max_w;
        if (Out > 3)
	  printf ("KL New best cutsize : %d\n",best_cutsize);
        no_better_steps = 0;
      }
    }

    while (number_locked != best_locked)
    { vertex = locked_list[--number_locked];
      sour = part[vertex];
      dest = locked[vertex]-1;
      part[vertex] = dest;
      gain(vertex);
      locked[vertex] = 0;
      part_w[sour] -= (vw?vw[vertex]:1);
      part_w[dest] += (vw?vw[vertex]:1);
      if (ex[vertex]>0)
	bucket_in(B[part[vertex]],&(bucket_ele[vertex]),ex[vertex]-in[vertex]);
      for (j=ep[vertex]; j<ep[vertex+1]; j++)
        if (!locked[neig=edge[j]])
	{ gain(neig);
	  if (buckets(bucket_ele[neig]))
	  { if (ex[neig]==0)
	      bucket_del (&(bucket_ele[neig]));
	    else
	      bucket_new_key(&(bucket_ele[neig]),ex[neig]-in[neig]);
	  }
	  else if (ex[neig]>0)
	    bucket_in(B[part[neig]],&(bucket_ele[neig]),ex[neig]-in[neig]);
        }
    }
    while (number_locked)
    { vertex = locked_list[--number_locked];
      locked[vertex] = 0;
      gain(vertex);
      if (ex[vertex]>0)
	bucket_in(B[part[vertex]],&(bucket_ele[vertex]),ex[vertex]-in[vertex]);
    }

  } while (best_cutsize < cutsize);

  FREE(in,int,n);
  FREE(ex,int,n);
  FREE(orient,int,n);
  FREE(locked,int,n);
  FREE(locked_list,int,n);
  FREE(edges,int,p);
  FREE(part_w,int,p);
  for (i=0; i<p; i++)
    bucket_free(B[i]);
  FREE(B,BUCKETS*,p);
  FREE(bucket_ele,BUCKET_ELE,n);

  return 0;
}
*/
