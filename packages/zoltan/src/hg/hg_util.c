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


/*****************************************************************************/

void Zoltan_HG_HGraph_Init(
  HGraph *hgraph
)
{
  hgraph->info  = 0;
  hgraph->nVtx  = 0;
  hgraph->nEdge = 0;
  hgraph->nPin  = 0;
  hgraph->nDim  = 0;

  hgraph->EdgeWeightDim   = 0 ;
  hgraph->VertexWeightDim = 0 ;

  hgraph->hindex  = NULL;
  hgraph->hvertex = NULL;
  hgraph->vindex  = NULL;
  hgraph->vedge   = NULL;
  hgraph->vwgt    = NULL;
  hgraph->ewgt    = NULL;
  hgraph->coor    = NULL;
}

/*****************************************************************************/

void Zoltan_HG_Graph_Init(
  Graph *graph
)
{
  graph->info  = 0;
  graph->nVtx  = 0;
  graph->nEdge = 0;
  graph->nDim  = 0;

  graph->EdgeWeightDim   = 0 ;
  graph->VertexWeightDim = 0 ;

  graph->nindex  = NULL;
  graph->vtxdist = NULL;
  graph->neigh   = NULL;
  graph->vwgt    = NULL;
  graph->ewgt    = NULL;
  graph->coor    = NULL;
}


/****************************************************************************/
int Zoltan_HG_HGraph_Free(
  HGraph *hg
)
{
/* Frees all memory associated with a hypergraph; 
 * does not free the hypergraph itself. 
 */
  if (hg != NULL) {
    ZOLTAN_FREE ((void **) &hg->hindex);
    ZOLTAN_FREE ((void **) &hg->hvertex);
    ZOLTAN_FREE ((void **) &hg->vindex);
    ZOLTAN_FREE ((void **) &hg->vedge);
    ZOLTAN_FREE ((void **) &hg->vwgt);
    ZOLTAN_FREE ((void **) &hg->ewgt);
    ZOLTAN_FREE ((void **) &hg->coor);
  }
  return ZOLTAN_OK;
}

/****************************************************************************/
int Zoltan_HG_Graph_Free(
  Graph *g
)
{
/* Frees all memory associated with a graph; does not free the graph itself. */
  if (g != NULL) {
    ZOLTAN_FREE ((void **) &g->vtxdist);
    ZOLTAN_FREE ((void **) &g->nindex);
    ZOLTAN_FREE ((void **) &g->neigh);
    ZOLTAN_FREE ((void **) &g->vwgt);
    ZOLTAN_FREE ((void **) &g->ewgt);
    ZOLTAN_FREE ((void **) &g->coor);
  }
  return ZOLTAN_OK;
}


/****************************************************************************/
int Zoltan_HG_Info (
  ZZ *zz, 
  HGraph *hg
)
{ 
  int  	i, size, size_min=INT_MAX, size_max=0;

  if (zz->Debug_Level < ZOLTAN_DEBUG_ALL)
    return ZOLTAN_WARN;

  puts("---------- HGraph Information (min/ave/max/tot) --------------------");
  printf("info:%d |V|=%d |E|=%d |P|=%d \n",hg->info,hg->nVtx,hg->nEdge,hg->nPin);

/* print weights */
  if (hg->vwgt)
  { float vwgt_min=FLT_MAX, vwgt_max=0, vwgt_tot=0;
    for (i=0; i<hg->nVtx; i++)
    { vwgt_tot += hg->vwgt[i];
      vwgt_min = MIN(vwgt_min,hg->vwgt[i]);
      vwgt_max = MAX(vwgt_max,hg->vwgt[i]);
    }
    printf("Vertex weights      :    %9.2f %9.2f %9.2f %12.2f\n",vwgt_min,
     vwgt_tot/hg->nVtx,vwgt_max,vwgt_tot);
  }
  if (hg->ewgt)
  { float ewgt_min=FLT_MAX, ewgt_max=0.0, ewgt_tot=0.0;
    for (i=0; i<hg->nEdge; i++)
    { ewgt_tot += hg->ewgt[i];
      ewgt_min = MIN(ewgt_min,hg->ewgt[i]);
      ewgt_max = MAX(ewgt_max,hg->ewgt[i]);
    }
    printf("HEdge weights       :    %9.2f %9.2f %9.2f %12.2f\n",ewgt_min,
     ewgt_tot/hg->nEdge,ewgt_max,ewgt_tot);
  }

/* print sizes */
  if (hg->hindex)
  { size_min = INT_MAX;
    size_max = 0;
    for (i=0; i<hg->nEdge; i++)
    { size = hg->hindex[i+1] - hg->hindex[i];
      size_min = MIN(size_min,size);
      size_max = MAX(size_max,size);
    }
    printf("Edge sizes          :    %6d    %9.2f %6d    %9d\n",size_min,
     (float)(hg->nPin)/hg->nEdge,size_max,hg->nPin);
  }
  if (hg->vindex)
  { size_min = INT_MAX;
    size_max = 0;
    for (i=0; i<hg->nVtx; i++)
    { size = hg->vindex[i+1] - hg->vindex[i];
      size_min = MIN(size_min,size);
      size_max = MAX(size_max,size);
    }
    printf("Vertex sizes        :    %6d    %9.2f %6d    %9d\n",size_min,
     (float)(hg->nPin)/hg->nEdge,size_max,hg->nPin);
  }

  puts("--------------------------------------------------------------------");
  return ZOLTAN_OK;
}




/****************************************************************************/
/**  Given arrays to lookup a vertex for a given hyperedge (hindex, hvertex)
 **  or arrays to lookup a hyperedge given a vertex (vindex, vedge) create the
 **  other set of arrays (mirror).
 **  The mirror (to be filled in) should have NULL pointers in hg on entry!
 **  The calling program is responsible for freeing the mirror's memory  ***/

int Zoltan_HG_Create_Mirror (
  ZZ *zz, 
  HGraph *hg
)
{
   int i, j ;                  /* loop counters */
   int inlength, outlength ;   /* input/output array lengths */
   int *index, *data ;         /* pointers to input information */
   int *outindex, *outdata ;
   char *yo = "Zoltan_HG_Create_Mirror" ;

   ZOLTAN_TRACE_ENTER(zz, yo);

   /* determine which data to "mirror" and set corresponding data pointers. */
   if (hg != NULL &&  hg->hindex != NULL && hg->hvertex != NULL
    && hg->vindex == NULL && hg->vedge == NULL)
      {
      ZOLTAN_TRACE_DETAIL(zz, yo, "Have hindex; building vindex.");
      inlength  = hg->nEdge ;
      outlength = hg->nVtx ;
      index     = hg->hindex ;
      data      = hg->hvertex ;
      outindex  = hg->vindex = (int*) ZOLTAN_MALLOC (sizeof(int) * (hg->nVtx+1));
      outdata   = hg->vedge  = (int*) ZOLTAN_MALLOC (sizeof(int) *  hg->nPin) ;
      if (outindex == NULL || outdata == NULL)
         {
         ZOLTAN_FREE ((void **) &(hg->vindex)) ;
         ZOLTAN_FREE ((void **) &(hg->vedge)) ;
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         ZOLTAN_TRACE_EXIT(zz, yo);
         return ZOLTAN_MEMERR;
         }

      }
   else if (hg != NULL && hg->vindex != NULL && hg->vedge != NULL
    && hg->hindex == NULL && hg->hvertex == NULL)
      {
      ZOLTAN_TRACE_DETAIL(zz, yo, "Have vindex; building hindex.");
      inlength  = hg->nVtx ;
      outlength = hg->nEdge ;
      index     = hg->vindex ;
      data      = hg->vedge ;
      outindex  = hg->hindex  = (int*)ZOLTAN_MALLOC (sizeof(int) *(hg->nEdge+1));
      outdata   = hg->hvertex = (int*)ZOLTAN_MALLOC (sizeof(int) * hg->nPin) ;
      if (outindex == NULL || outdata == NULL)
         {
         ZOLTAN_FREE ((void **) &(hg->hindex)) ;
         ZOLTAN_FREE ((void **) &(hg->hvertex)) ;
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         ZOLTAN_TRACE_EXIT(zz, yo);
         return ZOLTAN_MEMERR;
         }
      }
   else {
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_WARN ;  /* unable to proceed */
   }

   /* count number of data objects in out index space */
   for (i = 0 ; i < outlength+1 ; i++)
      outindex[i] = 0 ;
   for (i = 0 ; i < inlength ; i++)
      for (j = index[i] ;  j < index[i+1] ; j++)
          (outindex[data[j]+1])++ ;

   /* compute partial sum */
   for (i = 1 ; i < outlength ; i++)
      outindex [i] += outindex[i-1] ;

   /* store data in out index location and increment index */
   for (i = 0 ; i < inlength ; i++)
      for (j = index[i] ; j < index[i+1] ; j++)
         outdata[ (outindex[data[j]])++ ] = i ;

   /* out index now points past end of sublist. Shift array one position up */
   for (i = outlength ; i > 0 ; i--)
      outindex[i] = outindex[i-1] ;
   outindex[0] = 0 ;

   ZOLTAN_TRACE_EXIT(zz, yo);
   return ZOLTAN_OK ;
}


/****************************************************************************/
/* check that (hindex, hvertex) and (vindex, vedge) are consistant mappings */
/* returns ZOLTAN_WARN if not consistant, returns ZOLTAN_OK otherwise */
int Zoltan_HG_Check (
  ZZ *zz,
  HGraph *hg
)
{
   int iedge ;               /* iedge denotes an hyperedge */
   int j ;                   /* j is the index of a vertex */
   int k ;                   /* k is an index of hyperedge */

   if (hg->hindex == NULL || hg->hvertex == NULL
    || hg->vindex == NULL || hg->vedge == NULL)
       return ZOLTAN_WARN ;

   /* starting from (hindex,hvertex), for each edge determine each associated */
   /* vertex.  Then for each vertex lookup associated edges using (vindex, vedge) */
   for (iedge = 0 ; iedge < hg->nEdge ; iedge++)                  /* for each hyperedge */
      for (j = hg->hindex[iedge] ; j < hg->hindex[iedge+1] ; j++)  /* get index to vertices */
         {
         for (k = hg->vindex[hg->hvertex[j]] ;   /* for each vertex of current hyperedge */
              k < hg->vindex[hg->hvertex[j]+1] ; k++)  /* get index to hyperedges */
                 if (hg->vedge[k] == iedge)           /* use index value to find an edge */
                    break ;                          /* still looking for agreement */
         if (k == hg->vindex[hg->hvertex[j]+1])        /* did we succeed for this hyperedge? */
            return ZOLTAN_WARN ;                     /* if yes, carry on. if no, failure exit */
         }
   return ZOLTAN_OK ;
}



/****************************************************************************/
int Zoltan_HG_HGraph_to_Graph(
  ZZ *zz,
  HGraph *hg, 
  Graph *g
)
{ 
  int 	i, j, k, e, roughly_e, pins, *_neigh, *degrees, *next, empty=0;
  float	*w, *_ewgt, weight;
  char *yo = "Zoltan_HG_HGraph_to_Graph" ;

  Zoltan_HG_Graph_Init(g);
  g->info = hg->info;
  g->nVtx = hg->nVtx;
  g->vtxdist = NULL;

/* copy vertex weights */
  if (hg->vwgt)
  { g->vwgt = (float *) ZOLTAN_MALLOC (sizeof (float) * g->nVtx);
    if (g->vwgt == NULL)
    { ZOLTAN_FREE ((void **) &g->vwgt) ;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }
    for (i=0; i<g->nVtx; i++)
      g->vwgt[i] = hg->vwgt[i];
  }
  else
    g->vwgt = NULL;

/* calculate the roughly |E| and each degree */
  roughly_e = 0;
  degrees = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx);
  for (i = 0 ; i < hg->nVtx ; i++)
     degrees[i] = 0 ;
  for (i=0; i<hg->nEdge; i++)
  /* if (hg->hindex[i+1]-hg->hindex[i] <= 10) */
  { pins = hg->hindex[i+1]-hg->hindex[i];
    roughly_e += pins*(pins-1);
    pins--;
    for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      degrees[hg->hvertex[j]] += pins;
  }

/* Calculate the initial nindex */
  e = 0;
  g->nindex = (int *) ZOLTAN_MALLOC (sizeof (int) * (hg->nVtx+1));
  next      = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx);
  for (i=0; i<hg->nVtx; i++)
  { next[i] = g->nindex[i] = e;
    e += degrees[i];
  }
  g->nindex[hg->nVtx] = e;
  ZOLTAN_FREE ((void **) &degrees);

/* Calculate the initial neigh and ewgt */
  _neigh = (int *)   ZOLTAN_MALLOC (sizeof (int) * roughly_e);
  _ewgt  = (float *) ZOLTAN_MALLOC (sizeof (float) * roughly_e);
  for (i=0; i<hg->nEdge; i++)
  /* if (hg->hindex[i+1]-hg->hindex[i] <= 10) */
  { pins = hg->hindex[i+1]-hg->hindex[i];
    /* weight = 1.0/(pins-1); */
    weight = 2.0/((pins-1)*pins);
    if (hg->ewgt)
      weight *= hg->ewgt[i];
    for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      for (k=j+1; k<hg->hindex[i+1]; k++)
      { _neigh[next[hg->hvertex[j]]] = hg->hvertex[k];
        _neigh[next[hg->hvertex[k]]] = hg->hvertex[j];
        _ewgt[next[hg->hvertex[j]]++] = _ewgt[next[hg->hvertex[k]]++] = weight;
      }
  }
  ZOLTAN_FREE ((void **) &next);

/* Compact identical icident edges and their weight */
  w = (float *) ZOLTAN_MALLOC (sizeof (float) * hg->nVtx);
  for (i=0; i<hg->nVtx; i++)
    w[i] = 0.0;
  for (i=0; i<hg->nVtx; i++)
  { int start=g->nindex[i], end=g->nindex[i+1];
    g->nindex[i] -= empty;
    for (j=start; j<end; j++)
      if (_ewgt[j] > 0.0)
        w[_neigh[j]] += _ewgt[j];
    for (j=start; j<end; j++)
      if (_ewgt[j] > 0.0  &&  w[_neigh[j]] > 0.0)
      { _neigh[j-empty] = _neigh[j];
        _ewgt[j-empty] = w[_neigh[j]];
        w[_neigh[j]] = 0.0;
      }
      else
        empty++;
  }
  g->nindex[hg->nVtx] -= empty;
  ZOLTAN_FREE ((void **) &w);

/* Reallocate to the exact size */
  g->nEdge = g->nindex[g->nVtx];

  g->neigh = (int *) ZOLTAN_MALLOC (sizeof (int) * g->nindex[hg->nVtx]);
  memcpy(g->neigh,_neigh,g->nindex[hg->nVtx]*sizeof(int));
  ZOLTAN_FREE ((void **) &_neigh);

  g->ewgt = (float *) ZOLTAN_MALLOC (sizeof (float) * g->nindex[hg->nVtx]);
  memcpy(g->ewgt,_ewgt,g->nindex[hg->nVtx]*sizeof(float));
  ZOLTAN_FREE ((void **) &_ewgt);

  return ZOLTAN_OK;
}

/****************************************************************************/
int Zoltan_HG_Graph_to_HGraph(
  ZZ *zz,   
  Graph *g,        /* Input graph */
  HGraph *hg       /* Ouput hypergraph */
)
{
/* 
 *  Converts a graph g into a hypergraph g.
 *  One hyperedge is created for each vertex of g.
 *  Hyperedge i consists of vertex i + all edge neighbors of i in graph g.
 */
char *yo = "Zoltan_HG_Graph_to_HGraph";
int i, j, k;
int *hindex = NULL, *hvertex = NULL;  /* temporary array pointers */
float *ewgts = NULL;                  /* temporary array pointers */
int ewgt_dim = zz->Edge_Weight_Dim;
int vwgt_dim = zz->Obj_Weight_Dim;
int cnt;
int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  Zoltan_HG_HGraph_Init(hg);
  hg->info  = g->info;
  hg->nVtx  = g->nVtx;
  hg->nEdge = g->nVtx;
  hg->nPin  = g->nVtx + g->nEdge;
  hg->nDim  = g->nDim;

  /* Copy coordinates from graph to hypergraph */
  if (g->coor != NULL) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Copying coordinates");
    cnt = hg->nVtx * g->nDim;        
    hg->coor = (double *) ZOLTAN_MALLOC(cnt * sizeof(double));
    if (hg->coor == NULL) {
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    memcpy(hg->coor, g->coor, cnt * sizeof(double));
  }

  /* Copy vertex weights from graph to hypergraph */
  if (g->vwgt != NULL) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Copying vertex weights");
    cnt = hg->nVtx * vwgt_dim;
    hg->vwgt = (float *) ZOLTAN_MALLOC(cnt * sizeof(float));
    if (hg->vwgt == NULL) {
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    memcpy(hg->vwgt, g->vwgt, cnt * sizeof(float));
  }

  /* Allocate memory for hypergraph edges */
  ZOLTAN_TRACE_DETAIL(zz, yo, "Allocating edge arrays");
  hindex = hg->hindex = (int *) ZOLTAN_MALLOC((hg->nEdge+1) * sizeof(int));
  hvertex = hg->hvertex = (int *) ZOLTAN_MALLOC(hg->nPin * sizeof(int));
  if (g->ewgt != NULL) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Allocating edge weights");
    cnt = hg->nEdge * ewgt_dim;
    ewgts = hg->ewgt = (float *) ZOLTAN_MALLOC(cnt * sizeof(float));
  }
  if ((hg->nEdge > 0 && hindex == NULL) || 
      (hg->nPin > 0 && hvertex == NULL) || 
      (g->ewgt != NULL && ewgts == NULL)) {
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* KDD  The following section will have to be reworked for parallel HGs. */
  ZOLTAN_TRACE_DETAIL(zz, yo, "Filling edge arrays");
  for (cnt = 0, i = 0; i < hg->nEdge; i++) {
    hindex[i] = cnt;
    hvertex[cnt++] = i;
    for (j = g->nindex[i]; j < g->nindex[i+1]; j++) {
      hvertex[cnt++] = g->neigh[j];
      for (k = 0; k < ewgt_dim; k++)
        ewgts[i * ewgt_dim + k] += g->ewgt[j * ewgt_dim + k];
    }
  } 
  hindex[hg->nEdge] = cnt;

  /* Sanity check */
  if (hg->nPin != cnt) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Sanity check failed: nPin != cnt.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Creating mirror");
  ierr = Zoltan_HG_Create_Mirror(zz, hg);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in building mirror.");
    goto End;
  }


End:
  if (ierr == ZOLTAN_MEMERR)
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory.");
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_FREE ((void **) &(hg->coor));
    ZOLTAN_FREE ((void **) &(hg->vwgt));
    ZOLTAN_FREE ((void **) &(hg->ewgt));
    ZOLTAN_FREE ((void **) &(hg->hindex));
    ZOLTAN_FREE ((void **) &(hg->hvertex));
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
