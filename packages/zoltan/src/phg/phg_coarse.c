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

#include "phg.h"


/* Procedure to coarsen a hypergraph based on a matching. All vertices of one
   match are clustered to a single vertex. All hyperedges are kept;
   identical hyperedges are not collapsed. 
   The array LevelMap is the mapping of
   the old vertices to the new vertices. It will be used to pass a partition
   of the coarse graph back to the original graph.                         */
   
int Zoltan_PHG_Coarsening
( ZZ     *zz,         /* the Zoltan data structure */
  HGraph *hg,         /* information about hypergraph, weights, etc. */
  int    *match,      /* Matching, Packing or Grouping array */
  HGraph *c_hg,       /* points to a working copy of hg structure */
  int    *LevelMap)   /* information to reverse coarsenings later */
{
  int i, j, vertex, edge, *ip, me, size, count;
  int *cmatch=NULL, *used_edges=NULL, *c_vindex=NULL, *c_vedge=NULL;
  int *listgno=NULL, *listlno=NULL, *displs=NULL, *each_size=NULL;
  float *c_ewgt=NULL, *pwgt;
  char *buffer=NULL, *rbuffer=NULL;
  PHGComm *hgc = hg->comm;
  char *yo = "Zoltan_PHG_Coarsening";
 
  ZOLTAN_TRACE_ENTER (zz, yo);  
  Zoltan_HG_HGraph_Init (c_hg);   /* inits working copy of hypergraph info */
  
  /* (over) estimate number of external matches that we to send data to */
  count = 1;
  for (i = 0; i < hg->nVtx; i++)
    if (match[i] < 0)
      ++count;
                
  if (!(cmatch    = (int*) ZOLTAN_MALLOC (hg->nVtx     * sizeof(int)))
   || !(listgno   = (int*) ZOLTAN_MALLOC (count        * sizeof(int)))
   || !(listlno   = (int*) ZOLTAN_MALLOC (count        * sizeof(int)))
   || !(displs    = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))
   || !(each_size = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int))))   {   
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     ZOLTAN_TRACE_EXIT (zz, yo);
     return ZOLTAN_MEMERR;
  }    
  
  if (hg->vwgt == NULL && hg->VtxWeightDim > 0)
     {
     hg->vwgt = (float*) ZOLTAN_MALLOC (hg->nVtx * hg->VtxWeightDim 
      * sizeof(float));
     for (i = 0; i < hg->nVtx * hg->VtxWeightDim; i++)
        hg->vwgt[i] = 1;
     }
  if (hg->VtxWeightDim > 0)
    c_hg->vwgt = (float*) ZOLTAN_MALLOC (hg->nVtx * hg->VtxWeightDim
     * sizeof(float));   
 
  for (i = 0; i < hg->nVtx; i++)
     cmatch[i] = match[i];         /* working copy of match array */
       
  /* Assume all rows in a column have the entire (column's) matching info */
  /* Calculate the number of resulting coarse vertices. */
  c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
  me = hgc->myProc_x;             /* short name, convenience variable */
  size  = 0;                      /* size (in ints) to communicate */
  count = 0;                      /* number of vertices to communicate */
  for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertice */
    if (match[i] < 0)  {               /* external processor match */
      int proc, gx = -match[i]-1;
      
      /* rule to determine "ownership" of coarsened vertices between procs */
      proc = ((gx + VTX_LNO_TO_GNO (hg,i)) & 1) ? MIN(gx, me) : MAX(gx, me);
      
      /* prepare to send data to owner or to receive data I will own */
      if (proc != me)   {             /* another processor owns this vertex */
        size += hg->vindex[i+1] - hg->vindex[i];  /* send buffer sizing */ 
        listgno[count]   = gx;                    /* listgno of vtx's to send */
        listlno[count++] = i;                     /* gno for destination */
        }                                         /* lno of my match to gno */
      else 
        c_hg->nVtx++;         /* myProc owns the matching across processors */ 
    }
    else if (cmatch[i] >= 0)  { /* allow (local only) packing and groupings */    
      c_hg->nVtx++;
      vertex = i;
      while (cmatch[vertex] >= 0)  {
        cmatch[vertex] = -cmatch[vertex] - 1;  /* flag this as done already */
        vertex         = -cmatch[vertex] - 1;  
      }
    }
  }

  /* size and allocate the send buffer */
  size += ((2 + hg->VtxWeightDim) * count);
  if (!(buffer = (char*) ZOLTAN_MALLOC (1 + size * sizeof(int))))   {  
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  }
        
  /* Message is list of <gno, gno's edge count, list of edge lno's> */
  ip = (int*) buffer;
  for (i = 0; i < count; i++)  {
    *ip++ = listgno[i];                              /* destination vertex gno */
    for (j = 0; j < hg->VtxWeightDim; j++)
       {
       pwgt = (float*) ip++;                                  /* vertex weight */
       *pwgt = hg->vwgt[listlno[i]*hg->VtxWeightDim+j] ;
       }
       
    *ip++ = hg->vindex[listlno[i]+1] - hg->vindex[listlno[i]];        /* count */
    for (j = hg->vindex[listlno[i]]; j < hg->vindex[listlno[i]+1]; j++)
      *ip++ = hg->vedge[j];    
  }    
  MPI_Allgather (&size, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);
  
  size = 0;
  for (i = 0; i < hgc->nProc_x; i++)
    size += each_size[i];
  if (!(rbuffer = (char*) ZOLTAN_MALLOC (1 + size * sizeof(int))))   {
    ZOLTAN_TRACE_EXIT (zz, yo);
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  displs[0] = 0;
  for (i = 1; i < hgc->nProc_x; i++)
     displs[i] = displs[i-1] + each_size[i-1];
                
  MPI_Allgatherv (buffer, each_size[me], MPI_INT, rbuffer, each_size, displs,
   MPI_INT, hgc->row_comm);            

  /* index all received data for rapid lookup */
  ip = (int*) rbuffer;
  for (i = 0; i < size; i++)  {
    if (VTX_TO_PROC_X (hg, ip[i]) == me)
      cmatch [VTX_GNO_TO_LNO (hg,ip[i])] = i;   
    i++;                    /* destination gno */
    i += hg->VtxWeightDim;  /* skip vertex weights */
    i += ip[i];             /* skip hyperedges */
  }
       
  if (!(used_edges = (int*) ZOLTAN_CALLOC (hg->nEdge,      sizeof(int)))
   || !(c_vindex   = (int*) ZOLTAN_MALLOC ((hg->nVtx+1)  * sizeof(int)))
   || !(c_vedge    = (int*) ZOLTAN_MALLOC (2 * hg->nPins * sizeof(int)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }
  
  if (hg->EdgeWeightDim > 0 && hg->ewgt == NULL)
     {
     hg->ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * hg->EdgeWeightDim 
      * sizeof(float));
     for (j = 0; j < hg->nEdge * hg->EdgeWeightDim; j++)
        hg->ewgt[i] = 1;
     }
  if (hg->EdgeWeightDim > 0)
     {
     c_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * hg->EdgeWeightDim 
      * sizeof(float));
     for (j = 0; j < hg->nEdge * hg->EdgeWeightDim; j++)
        c_hg->ewgt[i] = hg->ewgt[i];
     }

  /* Construct the LevelMap; match[vertex] is changed back to original value */
  /* Coarsen vertices (create vindex, vedge), sum up coarsened vertex weights */
  c_hg->nPins = 0;                      /* count of coarsened pins */
  c_hg->nVtx  = 0;                      /* count of coarsened vertices */
  for (i = 0; i < hg->nVtx; i++)  {
    if (match[i] < 0 && cmatch[i] < 0)  /* match to external vertex, don't own */
       LevelMap [i] = match[i];         /* negative value => external vtx */         
    else if (match[i] < 0) {            /* match from external vertex, I own */
       LevelMap[i] = c_hg->nVtx;        /* next available coarse vertex */
       c_vindex[c_hg->nVtx] = c_hg->nPins;
       ip =  ((int*) rbuffer) + cmatch[i];      /* point to received data */
       ip++;                                    /* skip over vtx gno */
       for (j = 0; j < hg->VtxWeightDim; j++)
          {
          pwgt = (float*) ip++;
          c_hg->vwgt[c_hg->nVtx*hg->VtxWeightDim+j] = *pwgt;
          }  
       count = *ip++;           /* extract edge count, advance to first edge */
       while (count--)  {
          edge = *ip++;           
          used_edges [edge]      = i+1;
          c_vedge[c_hg->nPins++] = edge;
       }
  
       for (j = 0; j < hg->VtxWeightDim; j++)
         c_hg->vwgt[c_hg->nVtx * hg->VtxWeightDim + j]
          = hg->vwgt[vertex    * hg->VtxWeightDim + j] ;
            
       for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)  {
         if (used_edges [hg->vedge[j]] <= i)   {
           used_edges [hg->vedge[j]]     = i+1;          
           c_vedge[c_hg->nPins++] = hg->vedge[j];
         }      
       }        
       c_hg->nVtx++;          
    }
    else if (match[i] >= 0 && cmatch[i] < 0) { /* match/pack/group local vtx's */
      c_vindex[c_hg->nVtx] = c_hg->nPins;
      vertex = i;
      while (cmatch[vertex] < 0)  {    
        LevelMap[vertex] = c_hg->nVtx; 
        
        for (j = 0; j < hg->VtxWeightDim; j++)
          c_hg->vwgt[c_hg->nVtx * hg->VtxWeightDim + j]
           = hg->vwgt[vertex    * hg->VtxWeightDim + j] ;        
           
        for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++)  {
          if (used_edges [hg->vedge[j]] <= i)  {
            used_edges [hg->vedge[j]] = i+1;                  
            c_vedge[c_hg->nPins++] = hg->vedge[j];
          }              
        }                       
        cmatch[vertex] = -cmatch[vertex] - 1;
        vertex         =  cmatch[vertex];
      }
      c_hg->nVtx++;
    }
  }
  c_vindex[c_hg->nVtx] = c_hg->nPins;
  ZOLTAN_FREE ((void**) &used_edges);
  
  /* Update the dist_x array after coarsening: allocate and fill */
  if (!(c_hg->dist_x = (int*)ZOLTAN_MALLOC ((hgc->nProc_x+1) * sizeof(int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  } 
    
  size = c_hg->nVtx;
  MPI_Allgather (&size, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);
  
  c_hg->dist_x[0] = 0;
  for (i = 1; i < hgc->nProc_x; i++)
    c_hg->dist_x[i] = c_hg->dist_x[i-1] + each_size[i-1];
  size = 0;
  for (i = 0; i < hgc->nProc_x; i++)
    size += each_size[i];
  c_hg->dist_x[hgc->nProc_x] = size;  

  /* Assuming that we do not collapse Edges, dist_y for the coarse hgraph
   * is the same as dist_y for the fine hgraph */
  if (!(c_hg->dist_y = (int*)ZOLTAN_MALLOC ((hgc->nProc_y+1) * sizeof(int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  }  
  for (i = 0; i < hgc->nProc_y+1; i++)
    c_hg->dist_y[i] = hg->dist_y[i];

  /* Done if there are no remaining vertices */
  if (c_hg->nVtx == 0)  {
    if (!(c_hg->vindex = (int*) ZOLTAN_CALLOC (1, sizeof(int))))  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
    };
  }
  else  {
    c_hg->dist_y  = hg->dist_y;
    c_hg->ratio   = hg->ratio;
    c_hg->vmap    = hg->vmap;
    c_hg->hindex  = NULL;
    c_hg->hvertex = NULL;
    c_hg->coor    = hg->coor;  /* ??? needs to be fixed */
    c_hg->nDim    = hg->nDim;    
    c_hg->ewgt    = c_ewgt; 
    c_hg->vindex  = c_vindex;
    c_hg->vedge   = c_vedge;
    c_hg->nEdge   = hg->nEdge;
    c_hg->comm    = hg->comm;
    c_hg->info    = hg->info + 1;     /* for debugging */
    c_hg->ratio   = hg->ratio;        /* for "global" recursive bisectioning */
    c_hg->redl    = hg->redl;         /* to stop coarsening near desired count */
    
    c_hg->VtxWeightDim  = MAX(hg->VtxWeightDim,1);
    c_hg->EdgeWeightDim = hg->EdgeWeightDim;
  }
  
  Zoltan_Multifree (__FILE__, __LINE__, 7, &buffer, &rbuffer, &listgno, &listlno,
   &cmatch, &displs, &each_size);  

  ZOLTAN_TRACE_EXIT (zz, yo);
  return Zoltan_HG_Create_Mirror(zz, c_hg);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
