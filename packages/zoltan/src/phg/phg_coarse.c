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
   match are clustered to one vertex. Identical hyperedges are collapsed to a
   single hyperedge with combined weight. The array LevelMap is the mapping of
   the old vertices to the new vertices. It will be used to pass a partition
   of the coarse graph back to the original graph.                         */
   
int Zoltan_PHG_Coarsening
( ZZ       *zz,         /* the Zoltan data structure */
  PHGraph  *hg,         /* information about hypergraph, weights, etc. */
  int      *match,      /* Matching, Packing or Grouping array */
  PHGraph  *c_hg,       /* points to a working copy of hg structure */
  int      *LevelMap)   /* information to reverse coarsenings later */
{
  int i, j, vertex, edge, *ip, me, size, count;
  int *cmatch=NULL, *list=NULL, *used_edges=NULL, *c_vindex=NULL, *c_vedge=NULL;
  float *c_ewgt=NULL;
  char *buffer=NULL, *rbuffer=NULL;
  int *displs=NULL, *each_size=NULL;
  PHGComm *hgc = hg->comm;
  char *yo = "Zoltan_PHG_Coarsening";

  ZOLTAN_TRACE_ENTER (zz, yo);
  
  Zoltan_PHG_PHGraph_Init (c_hg);   /* inits working copy of hypergraph info */
  c_hg->info  = hg->info + 1;      /* for debugging */
  c_hg->ratio = hg->ratio;         /* for "global" recursive bisectioning */
  c_hg->redl  = hg->redl;          /* to stop coarsening near desired count */
    
  if (!(cmatch    = (int*) ZOLTAN_MALLOC (hg->nVtx     * sizeof(int)))
   || !(list      = (int*) ZOLTAN_MALLOC (hg->nVtx     * sizeof(int)))
   || !(displs    = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))
   || !(each_size = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int))))  {
     Zoltan_Multifree (__FILE__, __LINE__, 3, &cmatch, displs, list);     
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     ZOLTAN_TRACE_EXIT (zz, yo);
     return ZOLTAN_MEMERR;
  }     
  for (i = 0; i < hg->nVtx; i++)
     cmatch[i] = match[i];  /* working copy of match array */
    
  /* Assume all rows in a column have the entire (column's) matching info */
  /* Calculate the number of coarse vertices. */
  c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
  me = hgc->myProc_x;             /* convenience variable */
  size = 0;                       /* number of ints to communicate */
  count = 0;                      /* number of messages to communicate */
  for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertice */
    if (match[i] < 0)  {               /* external processor match */
      int gx, proc;
      gx = -match[i] -1;
      
      /* rule to determine "ownership" of coarsened vertices across procs */
      proc = ((gx + VTX_LNO_TO_GNO (hg,i)) & 1) ? MIN(gx, me) : MAX(gx, me);
      if (proc != me)   {             /* another processor owns this vertex */
        size += hg->vindex[i+1] - hg->vindex[i];   
        list[count++] = i;
        }
      else 
        c_hg->nVtx++;         /* myProc owns the matching across processors */ 
    }
      
    /* allow for possible (local only) packing and groupings */    
    if (match[i] >= 0 && cmatch[i] >= 0)  {
      c_hg->nVtx++;
      vertex = i;
      while (cmatch[vertex] >= 0)  {
        j              =  cmatch[vertex];
        cmatch[vertex] = -cmatch[vertex] - 1;  /* flag this as done already */     
        vertex         =  j;
      }
    }
  }

  /* size and allocate the send buffer, and weight array */
  size += count;
  if (!(buffer = (char*) ZOLTAN_MALLOC (size * sizeof(int)))  
   || !(c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx, sizeof(float))))  {
    Zoltan_Multifree (__FILE__, __LINE__, 2, c_hg->vwgt, buffer);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ZOLTAN_TRACE_EXIT (zz, yo);
    return ZOLTAN_MEMERR;
  }
              
  MPI_Allgather (&size, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);  

  for (size = 0, i = 0; i < hgc->nProc_x; i++)
    size += each_size[i];

  if (!(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))   {
    ZOLTAN_FREE (&rbuffer);         
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }  

  displs[0] = 0;
  for (i = 1; i < hgc->nProc_x; i++)
     displs[i] = displs[i-1] + each_size[i-1];
        
  /* Message is list of <gno, gno's edge count, list of edge gno's> */
  ip = (int*) buffer;
  for (i = 0; i < count; i++)  {
     *ip++ = VTX_LNO_TO_GNO (hg, list[i]);       /* destination vertex gno */        
     *ip++ = hg->vindex[list[i]+1] - hg->vindex[list[i]];   /* count */                                                             /* weights??? */
     for (j = hg->vindex[list[i]]; j < hg->vindex[list[i]+1]; j++)
       *ip++ = EDGE_LNO_TO_GNO (hg, hg->vedge[j]);        /* edges */                           /* edges */
  }
          
  MPI_Allgatherv (buffer, count, MPI_INT, rbuffer, each_size, displs, MPI_INT,
   hgc->row_comm);            

  /* index all received data */
  ip = (int*) rbuffer;
  for (i = 0; i < size; i++)  {
    vertex = ip[i];
    if (VTX_TO_PROC_X (hg, vertex) == me)
      cmatch[VTX_GNO_TO_LNO(hg,vertex)] = i;
    i++;
    i += ip[i];  /* count of hyperedges */
  }
      
  if (!(used_edges = (int*)   ZOLTAN_CALLOC (c_hg->nEdge,     sizeof(int)))
   || !(c_ewgt     = (float*) ZOLTAN_MALLOC (hg->nEdge      * sizeof(float)))
   || !(c_vindex   = (int*)   ZOLTAN_MALLOC ((hg->nVtx+1)   * sizeof(int)))
   || !(c_vedge    = (int*)   ZOLTAN_MALLOC (hg->nNonZero   * sizeof(int)))) {
      Zoltan_Multifree (__FILE__, __LINE__, 7, c_hg->vwgt,
       buffer, rbuffer, &used_edges, &c_ewgt, &c_vindex, &c_vedge);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }  

  /* Construct the LevelMap; match[vertex] is changed back to original value */
  /* Coarsen vertices (create vindex, vedge), sum up coarsened vertex weights */
  c_hg->nNonZero = 0;   /* count of coarsened pins */
  c_hg->nVtx     = 0;   /* count of coarsened vertices */
  for (i = 0; i < hg->nVtx; i++)  {
    if (match[i] < 0 && cmatch[i] < 0)        /* match to external vertex */                  
       LevelMap [i] = match[i];         /* negative value => external vtx */         
    else if (match[i] < 0) {                /* match from external vertex */
       LevelMap[i] = c_hg->nVtx;
      
/*      c_hg->vwgt[c_hg->nVtx] += vlist[k]->weight;     */
       c_hg->vindex[c_hg->nVtx] = c_hg->nNonZero;
       
       ip = ((int*) rbuffer) + i;
       count = *++ip;
       for (j = 0; j < count; j++)  {
          edge = EDGE_GNO_TO_LNO (hg, *ip++);
          used_edges [edge]     = i+1;
          c_hg->vedge[c_hg->nNonZero++] = edge;
       }
       
       c_hg->vwgt[c_hg->nVtx] += hg->vwgt ? hg->vwgt[vertex] : 1.0;
       for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)  {
         if (used_edges [hg->vedge[j]] <= i)  {
           used_edges [hg->vedge[j]]     = i+1;          
           c_hg->vedge[c_hg->nNonZero++] = hg->vedge[j];
         }      
       }        
       c_hg->nVtx++;          
    }
    else if (match[i] >= 0 && cmatch[i] < 0)   /* match, pack, group my vtx's */
      c_hg->vindex[c_hg->nVtx] = c_hg->nNonZero;
      vertex = i;
      while (cmatch[vertex] < 0)  {    
        LevelMap[vertex] = c_hg->nVtx;    
        c_hg->vwgt[c_hg->nVtx] += hg->vwgt ? hg->vwgt[vertex] : 1.0;  

        for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++)  {
          if (used_edges [hg->vedge[j]] <= i)  {
            used_edges [hg->vedge[j]]     = i+1;          
            c_hg->vedge[c_hg->nNonZero++] = hg->vedge[j];
          }              
        }                       
        cmatch[vertex] = -cmatch[vertex] - 1;
        vertex         =  cmatch[vertex];
      }
      c_hg->nVtx++;
    }
    
  ZOLTAN_FREE ((void**) &used_edges);

  /* Done if there are no remaining vertices */
  if (c_hg->nVtx == 0)  {
    c_hg->ewgt = NULL;
    if (!(c_hg->vindex = (int*) ZOLTAN_CALLOC (1, sizeof(int))))  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
    }
    c_hg->vedge = NULL;
  }
  else  {
    c_hg->ewgt   = c_ewgt;
    c_hg->vindex = c_vindex;
    c_hg->vedge  = c_vedge;
  }
  
  Zoltan_Multifree (__FILE__, __LINE__, 10, c_hg->vwgt, buffer, rbuffer,
   &c_ewgt, &c_vindex, &c_vedge, list, cmatch, displs, each_size);  

  ZOLTAN_TRACE_EXIT (zz, yo);
  return Zoltan_PHG_Create_Mirror(zz, c_hg);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
