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

#include <assert.h>
#include "phg.h"


static ZOLTAN_PHG_MATCHING_FN matching_ipm;  /* inner product matching */
static ZOLTAN_PHG_MATCHING_FN matching_loc;  /* local ipm (in other words HCM:heavy connectivity matching) */


/* static void check_upper_bound_of_matching_weight (Graph*, ZZ*, Matching); */
/* static int graph_connected_components (int, int*, int*, int);             */

static void print_matrix(int*, int*, int, int, int);

/*****************************************************************************/
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams *hgp)
{
    int exist=1;
    
    if (!strcasecmp(hgp->redm_str, "no"))        hgp->matching = NULL;
    else if (!strcasecmp(hgp->redm_str, "loc"))  hgp->matching = matching_loc;    
    else if (!strcasecmp(hgp->redm_str, "ipm"))  hgp->matching = matching_ipm;
    else {
        exist = 0;
        hgp->matching = NULL;
    }
    
    return exist;
}



/*****************************************************************************/
int Zoltan_PHG_Matching (
  ZZ *zz,
  PHGraph *hg,
  Matching match,
  PHGPartParams *hgp)
{
float *old_ewgt = NULL, *new_ewgt = NULL;
int   err = ZOLTAN_OK;
char  *yo = "Zoltan_PHG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hg->vwgt && hgp->ews) {
     if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto End;
     }
     Zoltan_PHG_Scale_HGraph_Weight (zz, hg, new_ewgt, hgp->ews);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
  }

  /* Do the matching */
  if (hgp->matching)
     err = hgp->matching (zz, hg, match);

End:

  /* Restore the old edge weights */
  if (hg->vwgt && hgp->ews)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT (zz, yo);
  return err;
}



/* UVC:
   a simple HCM/IPM variant: root of each column procs find HCM using only its local data.
   This matching is implemented just for testing purposes.
 */
static int matching_loc(ZZ *zz, PHGraph *hg, Matching match)
{
    int i, j, *eweight, *adj, *visit, degzero=0, matchcnt=0;
    char *yo = "matching_loc";
    PHGComm *hgc=hg->comm;
    struct {
        int nNonZero;
        int rank;
    } rootin, root;

    /* find the index of the proc in column group with the most #nonzeros; it will be our root
       proc for computing moves since it has better knowedge about global hypergraph */
    rootin.nNonZero = hg->nNonZero; 
    rootin.rank = hgc->myProc_y;
    MPI_Allreduce(&rootin, &root, 1, MPI_2INT, MPI_MAXLOC, hgc->col_comm);

    uprintf(hgc, "root is %d with %d nonzero\n", root.rank, root.nNonZero);

    
    if (hgc->myProc_y==root.rank) { /* only root of each column does this */
        if (!(visit = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
            || !(adj = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
            || !(eweight = (int*) ZOLTAN_CALLOC (hg->nVtx, sizeof(int))) ) {
            Zoltan_Multifree (__FILE__, __LINE__, 3, &visit, &adj, &eweight);
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
            return ZOLTAN_MEMERR;
        }
        for (i=0; i<hg->nVtx; ++i)
            visit[i] = i;
        Zoltan_PHG_Rand_Perm_Int(visit, hg->nVtx);

        for (i = 0; i<hg->nVtx; ++i) {
            int v = visit[i];
            if (match[v] == v) {
                if ( hg->vindex[v] == hg->vindex[v+1] )
                    degzero++; // isolated vertex detection
                else {
                    int maxuv=-1, maxu=-1, adjsz=0, k;

                    for (j = hg->vindex[v]; j < hg->vindex[v+1]; ++j) {
                        int edge = hg->vedge[j];
                        for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; ++k) {
                            int u = hg->hvertex[k];
                            if (u == match[u] && u != v) {
                                if (eweight[u]==0)
                                    adj[adjsz++] = u;
                                ++eweight[u];
                            }
                        }
                    }
                    for (k = 0; k < adjsz; ++k) {
                        int u = adj[k];
                        if (eweight[u] > maxuv) {
                            maxu = u;
                            maxuv = eweight[u];
                        }
                        eweight[u] = 0;
                    }
                    if (maxu!=-1) {
                        //match the maximally connected one with the current vertex
                        match[v] = maxu;
                        match[maxu] = v;
                        ++matchcnt;
                    }
                }
            }
        }


        // match isolated vertices
        if (degzero) {
            int v = -1; // mark the first isolated vertex with -1
            for (i = 0; i< hg->nVtx; ++i)
                if (hg->vindex[i]==hg->vindex[i+1]) { /* degree zero */
                    if (v == -1)
                        v = i;
                    else {
                        match[i] = v;
                        match[v] = i;
                        v = -1;
                        ++matchcnt;
                    }                
                }
        }
        uprintf(hgc, "there are %d matches\n", matchcnt);
    }
    MPI_Bcast(match, hg->nVtx, MPI_INT, root.rank, hgc->col_comm);

    Zoltan_Multifree (__FILE__, __LINE__, 3, &visit, &adj, &eweight);
    MPI_Barrier(hgc->Communicator);
    return ZOLTAN_OK;
}
    



/****************************************************************************
 * inner product matching
 * Parallelized version of the serial algorithm (see hg_match.c)
 * Based on conversations with Rob Bisseling by Aaron Becker, UIUC, summer 2004
 * completed by R. Heaphy
 */
               
/*******************************************************************************
  Bob's Notes during development
  Assumption:  hg (for example hg->hindex) contains zero-based arrays
   and information about the local portion of the hypergraph.
  Assumption: given a local, zero-based number, its corresponding gno
   can be found or computed (Karen's stuff).
  Assumption: the array "match" contains only the local (to column)
   matching information 
*/
             
static int matching_ipm (ZZ *zz, PHGraph *hg, Matching match)
{
  int i, j, lno, loop, vertex, maxpsum, *psums, *tsums;
  int count, size, *ip, bestv, bestsum, edgecount, index, *cmatch;
  int NDO, NLOOP;
  int *select, pselect;
  int *m_gno;
  int *m_vindex, *m_vedge;  /* zero-based loopup of edges by vertex */
  int *m_bestsum, *m_bestv; /* column's best results for each matched vertex */
  char *buffer, *rbuffer;    /* send and rec buffers */
  int *displs, *each_size;
  PHGComm *hgc = hg->comm;  
  char  *yo = "matching_ipm";
    
  /* compute NLOOP as 1/2 * total vertices/total columns */
  NLOOP = hg->dist_x[hgc->nProc_x+1]/(2 * hgc->nProc_x);
  NDO = 1;
  
  for (i = 0; i < hg->nVtx; i++)
     match[i] = i;
  pselect = 0;  /* marks position (count) of my processed vertices */
       
  /* local slice of global matching array.  It uses local numbering (zero-based)
     initially, match[i] = i.  After matching, match[i]=i indicates an unmatched
     vertex. A matching between vertices i & j is indicated by match[i] = j &
     match [j] = i.  NOTE: a match to an off processor vertex is indicated my a
     negative number, -(gno+1), which must use global numbers (gno's).
  */

  if (!(psums  = (int*) ZOLTAN_CALLOC (hg->nVtx,  sizeof(int)))
   || !(tsums  = (int*) ZOLTAN_CALLOC (hg->nVtx,  sizeof(int)))
   || !(select = (int*) ZOLTAN_MALLOC (NDO      * sizeof(int)))
   || !(cmatch = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
   || !(displs = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int))))
     {
     Zoltan_Multifree (__FILE__, __LINE__, 5, &psums, &tsums, &select, &cmatch,
      displs);
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
     
  if (!(m_vindex  = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))
   || !(m_vedge   = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))
   || !(m_gno     = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))
   || !(m_bestsum = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))
   || !(m_bestv   = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int))))
     {
     Zoltan_Multifree (__FILE__, __LINE__, 5, &m_vindex, &m_vedge, &m_gno,
      &m_bestsum, &m_bestv);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
        
  /* Loop processing NDO vertices per column each pass. Each loop has 4 phases:
     Phase 1: send NDO candidates for global matching - horizontal communication
     Phase 2: sum  inner products, find best         - vertical communication
     Phase 3: return best sums to owning column      - horizontal communication
     Phase 4: return actual matches                  - horizontal communication
  */
  for (loop = 0; loop < NLOOP; loop++)
     {
     /************************ PHASE 1: ***************************************/
     
     /* Select next NDO unmatched vertices to globally match.  Here many choices
        can be made for alternative algorithms: sequential, random, weight order,
        (hypergraph) vertex size, etc.  This version uses seqential: pick the
        next NDO unmatched vertices in natural lno order.
     */ 
     for (count = 0; pselect < hg->nVtx && count < NDO; pselect++)
        if (match[pselect] == pselect)    /* unmatched */
           select[count++] = pselect;     /* select it */
     if (count < NDO)                     /* what if we have a short count? */
        {
        for (i = 0; i < hg->nVtx; i++)
            if (match[i] == i)
              break;       
        while (count < NDO)             /* find an unmatched vertex */
           select[count++] = i;         /* fill the rest of the array with it */
        }                               /* assert(count == NDO); */
        
     /* determine the size of the send buffer & allocate it */   
     count = 0;
     for (i = 0; i < NDO; i++)
        count += (hg->vindex[select[i]+1] - hg->vindex[select[i]]);
     count = (2 * NDO) + count; 
          
     if (!(buffer    = (char*) ZOLTAN_MALLOC (count        * sizeof(int)))
      || !(each_size = (int *) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int))))
         {
         Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &each_size);         
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
         }  
             
     MPI_Allgather (&count, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);  

     for (size = 0, i = 0; i < hgc->nProc_x; i++)
        size += each_size[i];

     if (!(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int)))) 
         {
         ZOLTAN_FREE (&rbuffer);         
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
         }  

     displs[0] = 0;
     for (i = 1; i < hgc->nProc_x; i++)
        displs[i] = displs[i-1] + each_size[i-1];
        
     /* Message is list of <gno, gno's edge count, list of edge gno's> */
     ip = (int*) buffer;
     for (i = 0; i < NDO; i++)
        {
        *ip++ = VTX_LNO_TO_GNO (hg, select[i]);               /* vertex gno */        
        *ip++ = hg->vindex[select[i]+1] - hg->vindex[select[i]];   /* count */
        for (j = hg->vindex[select[i]]; j < hg->vindex[select[i]+1]; j++)
           *ip++ = EDGE_LNO_TO_GNO (hg, hg->vedge[j]);                                   /* edges */
        }
          
     /* send NDO vertices/edges to all row neighbors */
     /* rec all globally transmitted vertices/edges  */    
     MPI_Allgatherv (buffer, count, MPI_INT, rbuffer, each_size, displs, 
      MPI_INT, hgc->row_comm);
           
     /************************ PHASE 2: ***************************************/    

     if (!(m_gno     = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))
      || !(m_vindex  = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))
      || !(m_bestsum = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))      
      || !(m_bestv   = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))      
      || !(m_vedge   = (int*) ZOLTAN_MALLOC (size * sizeof(int))))
         {
         Zoltan_Multifree (__FILE__, __LINE__, 5, &m_gno, &m_vindex, &m_vedge,
          &m_bestsum, &m_bestv);
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
         }
     for (i = 0; i < hg->nVtx; i++)
        cmatch[i] = match[i];
        
     /* extract data from rbuffer, place in vindex & vedge arrays */ 
     ip = (int*) rbuffer;   
     for (index = i = 0; i < hgc->nProc_x * NDO; i++)
        {
        m_gno   [i] = *ip++;                      
        m_vindex[i] = index;
        edgecount = *ip++;
        while (edgecount--)
           m_vedge[index++] = EDGE_GNO_TO_LNO (hg, *ip++);
        }      
 
     Zoltan_Multifree (__FILE__, __LINE__, 4, &buffer, &rbuffer, &each_size,
      &displs);            
      
     if (!(tsums = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int)))
      || !(psums = (int*) ZOLTAN_MALLOC (NDO * hgc->nProc_x * sizeof(int))))
         {
         Zoltan_Multifree (__FILE__, __LINE__, 2, &tsums, &psums);
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
         } 
      
     /* for each match vertex, compute all local partial inner products */
     for (vertex = 0; vertex < hgc->nProc_x * NDO; vertex++) 
        {
        for (i = 0; i < hg->nVtx; i++)
           tsums[i] = psums[i] = 0;
           
        for (i = m_vindex[vertex]; i < m_vindex[vertex+1]; i++)
           for (j = hg->hindex[m_vedge[i]]; j < hg->hindex[m_vedge[i]+1]; j++)
              psums [hg->hvertex[j]]++;          /* unweighted??? */
        
        /* if local vtx, remove self inner product which is a maximum */
        if (VTX_TO_PROC_X (hg, m_gno[vertex]) == hgc->myProc_x)
           psums [VTX_GNO_TO_LNO (hg, m_gno[vertex])] = 0;   
        
        /* Want to use sparse communication with explicit summing later but
           for now, all procs in my column have same complete inner products
        */ 
        MPI_Allreduce(psums, tsums, hg->nVtx, MPI_INT, MPI_SUM, hgc->col_comm);
                         
        /* each proc computes best, all rows in a column compute same answer */ 
        maxpsum = -1;
        for (i = 0; i < hg->nVtx; i++)
           if (tsums[i] > maxpsum  &&  cmatch[i] == i)
              {
              m_bestsum [vertex]   = tsums[i];
              m_bestv   [vertex]   = i;
              }
        cmatch [m_bestv[vertex]] *= -1;      /* match[i] != -i */     
        }    
           
     /************************ PHASE 3: **************************************/

     size = 3 * NDO * hgc->nProc_x;
     if (!(buffer  = (char*) ZOLTAN_MALLOC (size * sizeof(int)))
      || !(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))
         {
         Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);         
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
         }  
                      
     /* prepare send buffer */       
     ip = (int*) buffer; 
     for (vertex = 0; vertex < NDO * hgc->nProc_x; vertex++)
        {
        *ip++ = m_gno [vertex];
        *ip++ = VTX_LNO_TO_GNO (hg, m_bestv [vertex]);
        *ip++ = m_bestsum [vertex];
        }
        
     /* send/rec to all columns in my row */     
     MPI_Allgather (buffer, 3 * NDO * hgc->nProc_x, MPI_INT, 
                   rbuffer, 3 * NDO * hgc->nProc_x, MPI_INT, hgc->row_comm);
                   
     for (i = 0; i < hg->nVtx; i++)
        m_bestsum[i] = 0;
                      
     ip = (int*) rbuffer;
     for (i = 0; i < NDO * hgc->nProc_x; i++)
        {
        vertex  = *ip++;
        bestv   = *ip++;
        bestsum = *ip++;
        lno     =  VTX_LNO_TO_GNO (hg, vertex);
        
        if (VTX_TO_PROC_X (hg, m_gno[vertex]) == hgc->myProc_x
         && bestsum > m_bestsum [lno])
           {
           m_gno     [lno] = bestv;
           m_bestsum [lno] = bestsum;
           }
        }   
 
     Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);
     size = 2 * NDO * sizeof(int);        
     if (!(buffer  = (char*) ZOLTAN_MALLOC (size))
      || !(rbuffer = (char*) ZOLTAN_MALLOC (size * hgc->nProc_x)))
         {
         Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);         
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
         }  
                             
     ip = (int*) buffer;                  
     for (i = 0; i < hg->nVtx; i++)
        if (m_bestsum [i] > 0)
           {
           if (VTX_TO_PROC_X (hg, m_gno[i]) == hgc->myProc_x)
              {
              j = VTX_GNO_TO_LNO (hg, m_gno[i]);
              match[i] = j;
              match[j] = i;
              }
           else
              match[i] = -m_gno[i] - 1;
           *ip++ = VTX_LNO_TO_GNO (hg, i);
           *ip++ = m_gno[i];
           }
            
     /************************ PHASE 4: ***************************************/
          
     MPI_Allgather (buffer, 2 * NDO, MPI_INT, 
                   rbuffer, 2 * NDO, MPI_INT, hgc->row_comm);
     
     ip = (int*) rbuffer;                   
     for (i = 0; i < NDO * hgc->nProc_x; i++)
        {
        vertex  = *ip++;
        bestv   = *ip++; 
        if (VTX_TO_PROC_X (hg, bestv)  == hgc->myProc_x
         && VTX_TO_PROC_X (hg, vertex) != hgc->myProc_x)
             match [VTX_GNO_TO_LNO (hg, bestv)] = -vertex-1;
        }
     Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);                  
     } /* end of large loop over LOOP */
     
  Zoltan_Multifree (__FILE__, __LINE__, 4, &psums, &tsums, &select, &cmatch); 
  Zoltan_Multifree (__FILE__, __LINE__, 5, &m_vindex, &m_vedge, &m_gno,
   &m_bestsum, &m_bestv);
  return ZOLTAN_OK;
}



/******************************************************************************* 
 * Draws the sparse matrix representation of a hypergraph, with vertices as
 * rows and columns as edges.  Useful for debugging on small inputs */
static void print_matrix(int *index, int *data, int x, int y, int h)
{
    char *matrix = (char*) ZOLTAN_MALLOC (x * y * sizeof(char));
    int i, j;
 
    for (i = 0; i < y; ++i)
       for (j = 0; j < x; ++j)
          matrix[i * x + j] = '.';

    if (h)
       for (i = 0; i < y; ++i)
          for (j = index[i]; j < index[i + 1]; ++j)
             matrix[i * x + data[j]] = 'x';
    else
       for (i = 0; i < x; ++i)
          for (j = index[i]; j < index[i + 1]; ++j)
             matrix[data[j] * x + i] = 'x';
 
    for (i = 0; i < y; ++i) {
       for (j = 0; j < x; ++j)
          printf("%c ", matrix[i * x + j]);
       printf("\n");
       }
    printf("\n");
   
    ZOLTAN_FREE ((void**) &matrix);
    return;
}



/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
