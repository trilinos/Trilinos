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

static float hcut_size_total (HGraph *, Partition);
static float hcut_size_links (ZZ *, HGraph *, int, Partition);

/****************************************************************************/
/* Routine to set function pointers corresponding to input-string options. */
int Zoltan_HG_Set_Part_Options(ZZ *zz, HGPartParams *hgp)
{
  char *yo = "Zoltan_HG_Set_Part_Options";

  if (hgp->bal_tol < 0.0)
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid HG_BALANCE_TOLERANCE.");
    return ZOLTAN_FATAL;
  }

  /* Set reduction method. */
  hgp->matching = hgp->matching_opt = NULL ;
  hgp->packing  = hgp->packing_opt  = NULL ;
  hgp->grouping = hgp->grouping_opt = NULL ;

  if (!(Zoltan_HG_Set_Matching_Fn (hgp))
   && !(Zoltan_HG_Set_Packing_Fn  (hgp))
   && !(Zoltan_HG_Set_Grouping_Fn (hgp)))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid HG_REDUCTION_METHOD.");
    return ZOLTAN_FATAL;
  }

  /* Set global partitioning method */
  if (!(hgp->global_part = Zoltan_HG_Set_Global_Part_Fn(hgp->global_str)))
  { ZOLTAN_PRINT_ERROR(zz->Proc,yo,"Invalid HG_GLOBAL_PARTITIONING.");
    return ZOLTAN_FATAL;
  }

  /* Set local refinement method. */
  if (!(hgp->local_ref = Zoltan_HG_Set_Local_Ref_Fn(hgp->local_str)))
  { ZOLTAN_PRINT_ERROR(zz->Proc,yo,"Invalid HG_LOCAL_REFINEMENT.");
    return ZOLTAN_FATAL;
  }

  return ZOLTAN_OK;
}

/****************************************************************************/
/*  Main partitioning function for hypergraph partitioning. */
int Zoltan_HG_HPart_Lib (
  ZZ *zz,              /* Zoltan data structure */
  HGraph *hg,          /* Input hypergraph to be partitioned */
  int p,               /* Input:  number partitions to be generated */
  Partition part,      /* Output:  partition #s; aligned with vertex arrays. */
  HGPartParams *hgp        /* Input:  parameters for hgraph partitioning. */
)
{ int  i;
  int  ierr = ZOLTAN_OK;
  char *yo = "Zoltan_HG_HPart_Lib" ;

  if (!part)
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Output partition array is NULL.");
    return ZOLTAN_MEMERR;
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST)
  {  printf("START %3d |V|=%6d |E|=%6d %d/%s/%s-%s p=%d...\n", hg->info,hg->nVtx,
      hg->nEdge, hgp->redl, hgp->redm_str, hgp->global_str, hgp->local_str, p);
    if (zz->Debug_Level > ZOLTAN_DEBUG_LIST)
    { ierr = Zoltan_HG_Info(zz, hg);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
        return ierr;
    }
  }

  if (hgp->redl < p)
    hgp->redl = p;

  if (p <= 0)
  { if (zz->Debug_Level > ZOLTAN_DEBUG_LIST)
    { char msg[128];
      sprintf(msg, "PART ERROR...p=%d is not a positive number!\n", p);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    }
    return ZOLTAN_FATAL;
  }

  /* take care of all special cases first */
  if (p == 1)             /* everything goes in the one partition */
  { for (i=0; i<hg->nVtx; i++)
      part[i] = 0 ;
  }
  else if (p >= hg->nVtx) /* more partitions than vertices, trivial answer */
  { for (i=0; i<hg->nVtx; i++)
      part[i] = i;
  }
  else if (hg->nVtx <= hgp->redl) /* fewer vertices than desired */
  { ierr = Zoltan_HG_Global (zz, hg, p, part, hgp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      return ierr;
  }
  else /* normal multilevel situation */
  { int    *pack=NULL, *LevelMap=NULL, *c_part=NULL, limit;
    HGraph c_hg;

    /* Allocate Packing Array (used for matching, packing & grouping) */
    if (!(pack = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
    { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }

    /* Calculate Packing, Grouping or Matching */
    limit = hg->nVtx-hgp->redl;
    if (hgp->packing)
      ierr = Zoltan_HG_Packing (zz,hg,pack,hgp,&limit);
    else if (hgp->grouping)
      ierr = Zoltan_HG_Grouping (zz,hg,pack,hgp,&limit) ;
    else
    { Graph g;
      Zoltan_HG_Graph_Init(&g);
      ierr = Zoltan_HG_HGraph_to_Graph (zz, hg, &g) ;
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      { ZOLTAN_FREE((void **) &pack);
        Zoltan_HG_Graph_Free(&g);
        return ierr ;
      }
      ierr = Zoltan_HG_Matching(zz,hg,&g,pack,hgp,&limit) ;
      Zoltan_HG_Graph_Free(&g);
    }
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    { ZOLTAN_FREE((void **) &pack);
      return ierr ;
    }

    /* Allocate and initialize LevelMap */
    if (!(LevelMap = (int *) ZOLTAN_CALLOC (hg->nVtx,sizeof(int))))
    { ZOLTAN_FREE ((void **) &pack);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }

    /* Construct coarse hypergraph and LevelMap */
    ierr = Zoltan_HG_Coarsening(zz, hg,pack,&c_hg,LevelMap);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    { ZOLTAN_FREE ((void **) &pack);
      ZOLTAN_FREE ((void **) &LevelMap);
      return ierr ;
    }

    /* Check the consistency of the coarsening */
    if (limit != c_hg.nVtx-hgp->redl)
    { printf("limit %d is not %d-%d!\n",limit,c_hg.nVtx,hgp->redl);
      return ZOLTAN_FATAL;
    } 
    else if (c_hg.nVtx < hgp->redl)
    { printf("wanted coarsen to %d vertices, but reached %d vertices.\n",hgp->redl,c_hg.nVtx);
      return ZOLTAN_FATAL;
    }
    if (c_hg.nVtx > .9*hg->nVtx)
      hgp->redl = c_hg.nVtx;
    
    /* Free Packing */
    ZOLTAN_FREE ((void **) &pack);

    /* Allocate Partition of coarse graph */
    if (!(c_part = (int *) ZOLTAN_CALLOC (c_hg.nVtx,sizeof(int))))
    { ZOLTAN_FREE ((void **) &LevelMap);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }
  
    /* Recursively partition coarse hypergraph */
    ierr = Zoltan_HG_HPart_Lib (zz,&c_hg,p,c_part,hgp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    { ZOLTAN_FREE ((void **) &c_part);
      ZOLTAN_FREE ((void **) &LevelMap);
      return ierr ;
    }

    /* Project coarse partition to fine partition */
    for (i=0; i<hg->nVtx; i++)
      part[i] = c_part[LevelMap[i]];
  
    /* Free coarse graph, coarse partition and LevelMap */
    Zoltan_HG_HGraph_Free (&c_hg);
    ZOLTAN_FREE ((void **) &c_part);
    ZOLTAN_FREE ((void **) &LevelMap);
  }

  /* Locally refine partition */
  ierr = Zoltan_HG_Local (zz, hg, p, part, hgp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    return ierr ;

  /* print useful information (conditionally) */
  if (zz->Debug_Level > ZOLTAN_DEBUG_LIST)
  { ierr = Zoltan_HG_HPart_Info (zz,hg,p,part);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      return ierr ;
  }
  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST)
    printf("FINAL %3d |V|=%6d |E|=%6d %d/%s/%s-%s p=%d cut=%.2f\n",
     hg->info,hg->nVtx,hg->nEdge,hgp->redl,hgp->redm_str,hgp->global_str,
     hgp->local_str, p,hcut_size_links(zz,hg,p,part));

  return ierr ;
}

/****************************************************************************/

static float hcut_size_total (HGraph *hg, Partition part)
{ int   i, j, hpart;
  float cut=0.0;

  for (i=0; i<hg->nEdge; i++)
  { hpart = part[hg->hvertex[hg->hindex[i]]];
    for (j=hg->hindex[i]+1; j<hg->hindex[i+1] && part[hg->hvertex[j]]==hpart; j++);
    if (j != hg->hindex[i+1])
      cut += (hg->ewgt?hg->ewgt[i]:1.0);
  }
  return cut;
}

/****************************************************************************/

static float hcut_size_links (ZZ *zz, HGraph *hg, int p, Partition part)
{ int   i, j, *parts, nparts;
  float cut=0.0;
  char *yo = "hcut_size_links" ;

  if (!(parts = (int *) ZOLTAN_CALLOC (p,sizeof(int))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<hg->nEdge; i++)
  { nparts = 0;
    for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
    { if (parts[part[hg->hvertex[j]]] < i+1)
        nparts++;
      parts[part[hg->hvertex[j]]] = i+1;
    }
    cut += (nparts-1)*(hg->ewgt?hg->ewgt[i]:1.0);
  }
  ZOLTAN_FREE ((void **) &parts);
  return cut;
}

/****************************************************************************/

static int hmin_max (ZZ *zz, int P, int *q)
{ int   i, values[3];

  if (P > 0)
  { values[0] = INT_MAX;
    values[1] = values[2] = 0;
    for (i=0; i<P; i++)
    { values[2] += q[i];
      values[0] = MIN(values[0],q[i]);
      values[1] = MAX(values[1],q[i]);
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST)
      printf("%9d    %12.2f %9d    %9d\n",values[0], (float)(values[2])/P,
       values[1],values[2]);
    if (zz->Debug_Level > ZOLTAN_DEBUG_LIST)
    { for (i=0; i<P; i++)
        printf ("%d ",q[i]);
      printf("\n");
  } }
  return values[1];
}

/****************************************************************************/

static float hmin_max_float (ZZ *zz, int P, float *q)
{ int   i;
  float	values[3];

  if (P > 0)
  { values[0] = INT_MAX;
    values[1] = values[2] = 0;
    for (i=0; i<P; i++)
    { values[2] += q[i];
      values[0] = MIN(values[0],q[i]);
      values[1] = MAX(values[1],q[i]);
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST)
      printf("%12.2f %12.2f %12.2f %12.2f\n",values[0],values[2]/P,values[1],values[2]);
    if (zz->Debug_Level > ZOLTAN_DEBUG_LIST)
    { for (i=0; i<P; i++)
        printf ("%.2f ",q[i]);
      printf("\n");
  } }
  return values[1];
}

/****************************************************************************/

int Zoltan_HG_HPart_Info (ZZ *zz, HGraph *hg, int p, Partition part)
{ int	i, *size, max_size;
  char *yo = "hpart_info" ;

  if (zz->Debug_Level < ZOLTAN_DEBUG_LIST)
    return ZOLTAN_OK;

  puts("---------- Partition Information (min/ave/max/tot) ----------------");

  printf ("VERTEX-based:\n");
  if (!(size = (int *) ZOLTAN_CALLOC (p,sizeof(int))))
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<hg->nVtx; i++)
  { if (part[i]<0 || part[i]>=p)
    { fprintf(stderr, "PART ERROR...vertex %d has wrong part number %d\n",i,
       part[i]);
      return ZOLTAN_WARN;
    }
    size[part[i]]++;
  }
  printf (" Size               :"); max_size = hmin_max (zz, p, size);
  printf (" Balance            : %.3f\n", (float)max_size*p/hg->nVtx);
  ZOLTAN_FREE ((void **) &size);

  if (hg->vwgt)
  { float *size_w, max_size_w, tot_w=0.0;
    if (!(size_w = (float *) ZOLTAN_CALLOC (p,sizeof(float))))
    { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }
    for (i=0; i<hg->nVtx; i++)
    { size_w[part[i]] += (hg->vwgt[i]);
      tot_w += hg->vwgt[i];
    }
    printf (" Size w.            :");
    max_size_w = hmin_max_float (zz,p,size_w);
    printf (" Balance w.         : %.3f\n", max_size_w*p/tot_w);
    ZOLTAN_FREE ((void **) &size_w);
  }

  printf ("EDGE-based:\n");
  printf (" Cuts(total/links)  : %.3f %.3f\n",hcut_size_total(hg,part),hcut_size_links(zz,hg,p,part));
  puts("-------------------------------------------------------------------");
  return ZOLTAN_OK;
}

/****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
