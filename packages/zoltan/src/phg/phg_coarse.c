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
#include "zoltan_comm.h"

#define PIN_OVER_ALLOC 1.2  /* Overallocate pins by 20% */
#define PLAN_TAG 32010      /* tag for comm plan */


/* Procedure to coarsen a hypergraph based on a matching. All vertices of one
   match are clustered to a single vertex. Currently, we allow more
   than two vertices to be merged/grouped locally on a proc, but
   allow only pairs of two to be merged between processors.
   All hyperedges are kept; identical hyperedges are not collapsed. 
   The array LevelMap is the mapping of
   the old vertices to the new vertices. It will be used to pass a partition
   of the coarse graph back to the original graph.                         */
   
int Zoltan_PHG_Coarsening
( ZZ     *zz,         /* the Zoltan data structure */
  HGraph *hg,         /* information about hypergraph, weights, etc. */
  int    *match,      /* Matching, Packing or Grouping array */
  HGraph *c_hg,       /* output: the coarse hypergraph */
  int    *LevelMap,
  int    *LevelCnt,
  int   **LevelData)   /* information to reverse coarsenings later */
{
  int i, j, vertex, edge, me, size, count, nrecv, pincnt=hg->nPins;
  int *ip, *ip_old;
  int *cmatch=NULL, *used_edges=NULL, *c_vindex=NULL, *c_vedge=NULL;
  int *listgno=NULL, *listlno=NULL,  *listproc=NULL;
  int *msg_size=NULL;
  float *pwgt;
  char *buffer=NULL, *rbuffer=NULL;
  PHGComm *hgc = hg->comm;
  char *yo = "Zoltan_PHG_Coarsening";
  struct Zoltan_Comm_Obj  *comm_plan;
 
  ZOLTAN_TRACE_ENTER (zz, yo);  
  Zoltan_HG_HGraph_Init (c_hg);   /* inits working copy of hypergraph info */
  
  /* (over) estimate number of external matches that we need to send data to */
  count = 0;
  for (i = 0; i < hg->nVtx; i++)
    if (match[i] < 0)
      ++count;
 
  if (hg->nVtx > 0 && hgc->nProc_x > 0 && (
      !(cmatch    = (int*) ZOLTAN_MALLOC (hg->nVtx     * sizeof(int))))){
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     ZOLTAN_TRACE_EXIT (zz, yo);
     return ZOLTAN_MEMERR;
  }        
  if (count > 0 && (
      !(listgno   = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(listlno   = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(listproc  = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(msg_size  = (int*) ZOLTAN_MALLOC (count * sizeof(int)))
   || !(*LevelData= (int*) ZOLTAN_MALLOC (count * sizeof(int) * 2))))    {   
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     ZOLTAN_TRACE_EXIT (zz, yo);
     return ZOLTAN_MEMERR;
  }        
 
  for (i = 0; i < hg->nVtx; i++)
     cmatch[i] = match[i];         /* working copy of match array */
       
  /* Assume all rows in a column have the entire (column's) matching info */
  /* Calculate the number of resulting coarse vertices. */
  *LevelCnt   = 0;
  c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
  me = hgc->myProc_x;             /* short name, convenience variable */
  size  = 0;                      /* size (in ints) to communicate */
  count = 0;                      /* number of vertices to communicate */
  for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertice */
    if (match[i] < 0)  {               /* external processor match */
      int proc, gx;
      gx = -match[i]-1;
      proc = VTX_TO_PROC_X(hg,gx);
      
      /* rule to determine "ownership" of coarsened vertices between procs */
      proc = ((gx + VTX_LNO_TO_GNO (hg,i)) & 1) ? MIN(proc, me) : MAX(proc, me);
      
      /* prepare to send data to owner */
      if (proc != me)   {             /* another processor owns this vertex */
        size += hg->vindex[i+1] - hg->vindex[i];  /* send buffer sizing */ 
        listgno[count]   = gx;                    /* listgno of vtx's to send */
        listproc[count]  = proc;                  /* proc to send to */
        listlno[count++] = i;                     /* lno of my match to gno */
      }
      else   {       
        c_hg->nVtx++;         /* myProc owns the matching across processors */ 
        (*LevelData)[(*LevelCnt)++] = gx;
        (*LevelData)[(*LevelCnt)++] = i;        
      }
    }
    else if (cmatch[i] >= 0)  {     /* local matching, packing and groupings */    
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
  if (size > 0 && !(buffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))  {  
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  }

  /* Use communication plan for unstructured communication. */
        
  /* Message is list of <gno, vweights, gno's edge count, list of edge lno's> */
  /* We pack ints and floats into same buffer */
  ip = (int*) buffer;
  for (i = 0; i < count; i++)  {
    ip_old = ip;
    *ip++ = listgno[i];                            /* destination vertex gno */
    for (j = 0; j < hg->VtxWeightDim; j++) {
       /* EBEB Assume a float is no larger than an int. 
          This is usually true, but to be safe this trick should be avoided. */
       pwgt = (float*) ip++;                               /* vertex weight */
       *pwgt = hg->vwgt[listlno[i]*hg->VtxWeightDim+j] ;
    }
       
    *ip++ = hg->vindex[listlno[i]+1] - hg->vindex[listlno[i]];      /* count */
    for (j = hg->vindex[listlno[i]]; j < hg->vindex[listlno[i]+1]; j++)
      *ip++ = hg->vedge[j];    

    /* save mesg size in #ints */
    msg_size[i] = ip - ip_old;
  }    

  /* Create comm plan. */
  Zoltan_Comm_Create( &comm_plan, count, listproc, hgc->row_comm, PLAN_TAG, 
                      &nrecv);
  
  /* call Comm_Resize since we have variable-size messages */
  Zoltan_Comm_Resize( comm_plan, msg_size, PLAN_TAG+1, &size); 

  /* Allocate receive buffer. */
  /* size is the size of the received data, measured in #ints */
  if (size > 0 && !(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))   {
    ZOLTAN_TRACE_EXIT (zz, yo);
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Comm_Do sends personalized messages of variable sizes */
  Zoltan_Comm_Do(comm_plan, PLAN_TAG+2, buffer, sizeof(int), rbuffer);

  /* For now, destroy the plan. It would be better to save it and use
     it in the reverse for uncoarsening. */
  Zoltan_Comm_Destroy(&comm_plan);

  /* index all received data for rapid lookup */
  /* EBEB Not clear if this still works with comm plan! */
  ip = (int*) rbuffer;
  for (i = 0; i < size; i++)  {
    if (VTX_TO_PROC_X (hg, ip[i]) == me)
      cmatch [VTX_GNO_TO_LNO (hg, ip[i])] = i;     
    i++;                    /* destination gno */
    i += hg->VtxWeightDim;  /* skip vertex weights */
    i += ip[i];             /* skip hyperedges */
  }
  
  if (hg->nEdge>0 && !(used_edges = (int*)ZOLTAN_CALLOC(hg->nEdge,sizeof(int)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }
  if (hg->nVtx>0 && !(c_vindex = (int*)ZOLTAN_MALLOC((hg->nVtx+1)*sizeof(int)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }
  if (hg->nPins>0 && !(c_vedge = (int*)ZOLTAN_MALLOC (hg->nPins*sizeof(int)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }

  /* Allocate vertex weight array for coarse hgraph */
  if (hg->nVtx > 0 && hg->VtxWeightDim > 0) 
     c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx * hg->VtxWeightDim,
      sizeof(float));   
      
  /* Allocate edge weight array for coarse hgraph */
  if (hg->EdgeWeightDim > 0) {
    c_hg->ewgt =(float*)ZOLTAN_MALLOC(hg->nEdge*hg->EdgeWeightDim*sizeof(float));
    if (c_hg->ewgt == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
    }    
    for (i = 0; i < hg->nEdge * hg->EdgeWeightDim; i++)
      c_hg->ewgt[i] = hg->ewgt[i];
  }

  /* Construct the LevelMap; match[vertex] is changed back to original value */
  /* Coarsen vertices (create vindex, vedge), sum up coarsened vertex weights */
  /* EBEB This code constructs the vertex-based arrays. It might be
     bette to construct the edge-based arrays so we can easily remove edges. */

  c_hg->nPins = 0;                      /* count of coarsened pins */
  c_hg->nVtx  = 0;                      /* count of coarsened vertices */
  for (i = 0; i < hg->nVtx; i++)  {
    if (match[i] < 0 && cmatch[i] < 0)  /* match to external vertex, don't own */
      LevelMap [i] = match[i];         /* negative value => external vtx */         
    else if (match[i] < 0) {            /* match from external vertex, I own */
       LevelMap[i] = c_hg->nVtx;        /* next available coarse vertex */
       c_vindex[c_hg->nVtx] = c_hg->nPins;
       ip =  ((int*) rbuffer) + cmatch[i];      /* point to received data */               
       ip++;                                    /* skip over gno */
      
       for (j = 0; j < hg->VtxWeightDim; j++)  {
          pwgt = (float*) ip++;
          c_hg->vwgt[c_hg->nVtx*hg->VtxWeightDim+j] = *pwgt;
       }  
       count = *ip++;           /* extract edge count, advance to first edge */
       while (count--)  {
          edge = *ip++;           
          used_edges [edge] = i+1;
          if (c_hg->nPins >= pincnt)  {
             pincnt = 1 + PIN_OVER_ALLOC * pincnt;
             c_vedge = (int*) ZOLTAN_REALLOC (c_vedge, pincnt * sizeof(int));
          }
          c_vedge[c_hg->nPins++] = edge;
       }
  
      for (j = 0; j < hg->VtxWeightDim; j++)
        c_hg->vwgt[c_hg->nVtx * hg->VtxWeightDim + j]
         += hg->vwgt[i        * hg->VtxWeightDim + j] ;
            
      for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)  {
        if (used_edges [hg->vedge[j]] <= i)   {
          used_edges [hg->vedge[j]] = i+1;  
          if (c_hg->nPins >= pincnt)  {
             pincnt = 1 + PIN_OVER_ALLOC * pincnt;     
             c_vedge = (int*) ZOLTAN_REALLOC (c_vedge, pincnt * sizeof(int));
          }                  
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
           += hg->vwgt[vertex   * hg->VtxWeightDim + j] ;        
           
        for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++)  {
          if (used_edges [hg->vedge[j]] <= i)  {
            used_edges [hg->vedge[j]] = i+1; 
            if (c_hg->nPins >= pincnt)  {
               pincnt = 1 + PIN_OVER_ALLOC * pincnt;            
               c_vedge = (int*) ZOLTAN_REALLOC (c_vedge, pincnt * sizeof(int));
            }                             
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
    
  MPI_Allgather (&(c_hg->nVtx), 1, MPI_INT, &(c_hg->dist_x[1]), 1, MPI_INT, 
                 hgc->row_comm);
  
  /* dist_x is the cumulative sum */
  c_hg->dist_x[0] = 0;
  for (i = 0; i < hgc->nProc_x; i++)
    c_hg->dist_x[i+1] += c_hg->dist_x[i];

  /* Assuming that we do not collapse Edges, dist_y for the coarse hgraph
   * is the same as dist_y for the fine hgraph */
  /* EBEB We should remove hyperedges of size 1; then dist_y will change.  */
  if (!(c_hg->dist_y = (int*)ZOLTAN_MALLOC ((hgc->nProc_y+1) * sizeof(int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  }  
  for (i = 0; i < hgc->nProc_y+1; i++)
    c_hg->dist_y[i] = hg->dist_y[i];

  /* Done if there are no remaining vertices (on this proc) */
  if (c_hg->nVtx == 0)  {
    /* EBEB c_hg->vindex and vedge not properly initialized? */
    if (!(c_hg->vindex = (int*) ZOLTAN_CALLOC (1, sizeof(int))))  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
    };
  }
  else  {
    c_hg->vindex  = c_vindex;
    c_hg->vedge   = c_vedge;
  }

  /*
  c_hg->vmap    = (int *) ZOLTAN_MALLOC(hg->nVtx * sizeof(int));
  for (i = 0; i < hg->nVtx; i++) c_hg->vmap[i] = hg->vmap[i];
  */
  c_hg->vmap = NULL; /* UVC: we don't need vmap in the coarser graphs, it is only
                        needed in recursive bisection; and hence at level 0 */
  
  c_hg->hindex  = NULL;             /* is computed later by HG_Create_Mirror */
  c_hg->hvertex = NULL;             /* is computed later by HG_Create_Mirror */
  c_hg->coor    = NULL;             /* currently we don't use coordinates */
  c_hg->nDim    = hg->nDim;    
  c_hg->nEdge   = hg->nEdge;
  c_hg->comm    = hg->comm;
  c_hg->info    = hg->info + 1;     /* for debugging */
  c_hg->ratio   = hg->ratio;        /* for "global" recursive bisectioning */
  c_hg->redl    = hg->redl;         /* to stop coarsening near desired count */
    
  c_hg->VtxWeightDim  = hg->VtxWeightDim;
  c_hg->EdgeWeightDim = hg->EdgeWeightDim;
  
  Zoltan_Multifree (__FILE__, __LINE__, 7, &buffer, &rbuffer, &listgno, &listlno,
   &cmatch, &listproc, &msg_size);  
  ZOLTAN_TRACE_EXIT (zz, yo);
  return Zoltan_HG_Create_Mirror(zz, c_hg);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
