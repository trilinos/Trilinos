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

#include <stdlib.h>
#include "phg.h"


static ZOLTAN_PHG_MATCHING_FN matching_local;   /* function for local matching */
static ZOLTAN_PHG_MATCHING_FN matching_ipm;     /* inner product matching */
static ZOLTAN_PHG_MATCHING_FN matching_col_ipm; /* local ipm along proc columns*/


/*****************************************************************************/
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams *hgp)
{
    int exist=1;
    
    if (!strcasecmp(hgp->redm_str, "no"))
        hgp->matching = NULL;
    else if (!strncasecmp(hgp->redm_str, "l-", 2))  {
        HGPartParams hp;

        strcpy(hp.redm_str, hgp->redm_str+2);
        strcpy(hp.redmo_str, hgp->redmo_str);
        if (!Zoltan_HG_Set_Matching_Fn(&hp)) {
            exist = 0;
            hgp->matching = NULL;
        } else {   
            /* just to make sure that coarsening will continue. We'll not call
             * this code for global matching. Actually, we'll pick the best 
             * local, but code structure doesn't allow us to use a function */     
            hgp->matching = matching_local; 
            hgp->locmatching = hp.matching;
            hgp->matching_opt = hp.matching_opt;
        }
    } else if (!strcasecmp(hgp->redm_str, "c-ipm"))
        hgp->matching = matching_col_ipm;   
    else if (!strcasecmp(hgp->redm_str, "ipm"))
        hgp->matching = matching_ipm;
    else {
        exist = 0;
        hgp->matching = NULL;
    }
    
    return exist;
}



/*****************************************************************************/
int Zoltan_PHG_Matching (
  ZZ *zz,
  HGraph *hg,
  Matching match,
  PHGPartParams *hgp)
{
float *old_ewgt = NULL, *new_ewgt = NULL;
int   err = ZOLTAN_OK;
char  *yo = "Zoltan_PHG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hgp->ews) {
     if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto End;
     }
     /* EBEB We need a parallel edge weight scaling routine, 
             do not use the old serial routine! */
     Zoltan_HG_Scale_HGraph_Weight (zz, hg, new_ewgt, hgp->ews);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
  }

  /* Do the matching */
  if (hgp->locmatching) {  /* run local matching */
      int limit=hg->nVtx;
      PHGComm *hgc=hg->comm;
      int root_matchcnt, root_rank;
                
      if (hgp->matching)
          err = hgp->locmatching (zz, hg, match, &limit);

      /* Optimization */
      if (hgp->matching_opt) 
          err = hgp->matching_opt (zz, hg, match, &limit);

      /* find the index of the proc in column group with the best match
         (max #matches); it will be our root proc */
      Zoltan_PHG_Find_Root(hg->nVtx-limit, hgc->myProc_y, hgc->col_comm,
                           &root_matchcnt, &root_rank);

      MPI_Bcast(match, hg->nVtx, MPI_INT, root_rank, hgc->col_comm);

  } else if (hgp->matching) /* run global or column/row matching algorithms */
     err = hgp->matching (zz, hg, match);

End: 

  /* Restore the old edge weights */
  if (hgp->ews)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT (zz, yo);
  return err;
}


static int matching_local(ZZ *zz, HGraph *hg, Matching match)
{
    uprintf(hg->comm, "Something wrong! This function should not be called!\n");
    /* UVC: NOTE:
       The reason that we're not doing local matchin in this function, we don't
       have access to parameter structure. So there is no way to figure
       out which "local" matching needs to be called. Hence we do it
       in Zoltan_PHG_Matching */
    /* EBEB TODO: We should add the parameter struct as an input argument
       to all the matching routines. */
    return ZOLTAN_OK;
}


    
/* local inner product matching among vertices in each proc column */
/* code adapted from serial matching_ipm method */
static int matching_col_ipm(ZZ *zz, HGraph *hg, Matching match)
{
    int   i, j, v1, v2, edge, maxip, maxindex;
    int   matchcount=0;
    float *lips, *gips; /* local and global inner products */
    char  *yo = "matching_col_ipm";
    PHGComm *hgc = hg->comm;  

    lips = gips = NULL;

    if (!(lips = (float*) ZOLTAN_MALLOC(hg->nVtx * sizeof(float))) 
     || !(gips = (float*) ZOLTAN_MALLOC(hg->nVtx * sizeof(float)))){
        Zoltan_Multifree(__FILE__, __LINE__, 2, &lips, &gips);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }
    
    for (i = 0; i < hg->nVtx; i++)
        lips[i] = 0;
        
    /* for every vertex */
    for (v1 = 0; v1 < hg->nVtx; v1++) {
        if (match[v1] != v1)
            continue;

        /* for every hyperedge containing the vertex */
        for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
            edge = hg->vedge[i];
                
            /* for every other vertex in the hyperedge */
            for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
                v2 = hg->hvertex[j];
                /* 
                if(match[v2] != v2) {
                     row swapping goes here
                } 
                */
                lips[v2] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
            }
        }

        /* sum up local inner products along proc column */
        /* for now, sum up i.p. value for all vertices; this is slow! */
        /* to do: 1) ignore vertices already matched 
                  2) use sparse communication for a chunk of vertices */
        MPI_Allreduce(lips, gips, hg->nVtx, MPI_FLOAT, MPI_SUM, hgc->col_comm);              
        /* now choose the vector with greatest inner product */
        /* all processors in a column should get the same answer */
        maxip = 0;
        maxindex = -1;
        for (i = 0; i < hg->nVtx; i++) {
            v2 = i;
            if (gips[v2] > maxip && v2 != v1 && match[v2] == v2) {
                maxip = gips[v2];
                maxindex = v2;
            }
            lips[v2] = gips[v2] = 0;   /* clear for next iteration */
        }
        /* match if inner product > 0 */
        if (maxindex > -1 && maxip > 0) {
            match[v1] = maxindex;
            match[maxindex] = v1;
            matchcount++;
        } 
        
    }

    /*
    printf("Matched %d vertices\n", matchcount);
    printf("Final Matching:\n");
    for(i = 0; i < hg->nVtx; i++)
        printf("%2d ",i);
    printf("\n");
    for(i = 0; i < hg->nVtx; i++)
        printf("%2d ",match[i]);
    printf("\n");
    */

    Zoltan_Multifree(__FILE__, __LINE__, 2, &lips, &gips);
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
  
#define PHASE3 3      /* communication size is 3 ints */
#define PHASE4 2      /* communication size is 2 ints */
#define LOOP_FACTOR 8 /* forces ncandidates to have enough degrees of freedom */ 
#define PIN_FACTOR  2 /* oversizes a Malloc() to avoid future Realloc()s */
#define THRESHOLD 0   /* ignore inner products (i.p.) less than threshold */
       
static int matching_ipm (ZZ *zz, HGraph *hg, Matching match)
{
  int i, j, k, lno, loop, vertex, *order=NULL;
  int *c_order=NULL, pincnt, count, size, *ip, bestv, bestsum, edgecount, pins;
  int *cmatch=NULL, ncandidates, nrounds, *select=NULL, pselect;
  int *m_gno=NULL;
  int *m_vindex=NULL, *m_vedge=NULL;  /* zero-based loopup of edges by vertex */
  int *m_bestv=NULL;                 /* col's best results for matched vertex */
  int *b_gno=NULL;
  float *psums=NULL, *tsums=NULL, *m_bestsum=NULL, *b_bestsum=NULL;
  char *buffer=NULL, *rbuffer=NULL;    /* send and rec buffers */
  int *displs=NULL, *each_size=NULL, *each_count=NULL, total_count;
  PHGComm *hgc = hg->comm;  
  char  *yo = "matching_ipm";
       
   
  /* determine number of basic matching rounds scaled by number of procs */
  nrounds = hg->dist_x[hgc->nProc_x] / (hgc->nProc_x * LOOP_FACTOR);
  ncandidates = hg->nVtx/nrounds/2 + 1;  /* 2: each match removes 2 vertices */
        
  /* allocate storage proportional to number of local vertices */  
  if (hg->nVtx > 0 && (
      !(b_gno     = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
   || !(b_bestsum = (float*) ZOLTAN_MALLOC (hg->nVtx * sizeof(float)))  
   || !(order     = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))  
   || !(psums     = (float*) ZOLTAN_MALLOC (hg->nVtx * sizeof(float)))
   || !(tsums     = (float*) ZOLTAN_MALLOC (hg->nVtx * sizeof(float)))
   || !(cmatch    = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))))     {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
  }
    
  /* match[] is local slice of global matching array.  It uses local numbering 
   * (zero-based) initially, match[i] = i.  After matching, match[i]=i indicates
   * an unmatched vertex. A matching between vertices i & j is indicated by 
   * match[i] = j & match [j] = i.  NOTE: a match to an off processor vertex is
   * indicated my a negative number, -(gno+1), which must use global numbers
   * (gno's).        */
  for (i = 0; i < hg->nVtx; i++)
    match[i] = i;                /* initialize match array to all unmatched */
   
  /* order[] is used to implement alternative vertex visiting algorithms:
   * natural (lno), random, weight order, vertex size, etc. */
  for (i = 0; i < hg->nVtx; i++)
     order[i] = i;     /* select (visit) vertices in lno order */

  /* randomly select matching candidates from available vertices */
  Zoltan_Rand_Perm_Int (order, hg->nVtx);

  if (hgc->nProc_x > 0 && (
      !(select     = (int*) ZOLTAN_MALLOC (ncandidates  * sizeof(int)))
   || !(each_size  = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))
   || !(each_count = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))         
   || !(displs     = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))))  {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
  }
  
  /* need to determine global total of ncandidates to allocate storage */
  MPI_Allgather (&ncandidates, 1, MPI_INT, each_count, 1, MPI_INT, hgc->row_comm);
  total_count = 0;
  for (i = 0; i < hgc->nProc_x; i++)
    total_count += each_count[i];
  
  /* malloc m_vedge by PIN_FACTOR to avoid too many REALLOC()s later */
  if (hg->nPins > 0 && 
   !(m_vedge   = (int*) ZOLTAN_MALLOC  (hg->nPins * PIN_FACTOR * sizeof(int)))) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  pincnt = PIN_FACTOR * hg->nPins;
  
  /* Allocate storage proporational to global number of candidates */
  if (total_count > 0 &&  (
      !(c_order   = (int*)   ZOLTAN_MALLOC  (total_count    * sizeof(int)))
   || !(m_vindex  = (int*)   ZOLTAN_MALLOC ((total_count+1) * sizeof(int)))
   || !(m_gno     = (int*)   ZOLTAN_MALLOC  (total_count    * sizeof(int)))
   || !(m_bestsum = (float*) ZOLTAN_MALLOC  (total_count    * sizeof(float)))
   || !(m_bestv   = (int*)   ZOLTAN_MALLOC  (total_count    * sizeof(int)))))  {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
                        
  /* Loop processing ncandidates vertices per column each round.
   * Each loop has 4 phases:
   * Phase 1: send ncandidates vertices for global matching - horizontal communication
   * Phase 2: sum  inner products, find best        - vertical communication
   * Phase 3: return best sums to owning column     - horizontal communication
   * Phase 4: return actual match selections        - horizontal communication 
   *
   * No conflict resolution phase required because prelocking a potential match
   * prevents conflicts.
   */

  pselect = 0;                    /* marks position in vertices to be selected */
  for (loop = 0; loop < nrounds; loop++)  {
     
     /************************ PHASE 1: ***************************************/
     
     /* Select next ncandidates unmatched vertices to globally match. */      
     for (i = 0; i < hg->nVtx; i++)
       cmatch[i] = match[i];                                         
     for (count = 0; count < ncandidates && pselect < hg->nVtx; pselect++)
       if (cmatch[order[pselect]] == order[pselect])  {  /* unmatched */
         select[count++] = order[pselect];       /* select it */       
         cmatch[order[pselect]] = -1;            /* mark it as a pending match */
       }
     if (count < ncandidates)   {          /* what if we have a short count? */
       for (i = 0; i < hg->nVtx; i++)      /* find an unmatched vertex  */
         if (cmatch[i] == i)
           break;
       if (i < hg->nVtx)
         cmatch[i] = -1;                  
       while (count < ncandidates)       /* fill the rest of the array with it */
         select[count++] = (i == hg->nVtx) ? i-1 : i;
     }
                
     /* determine the size of the send buffer & allocate it */   
     size = 0;
     for (i = 0; i < ncandidates; i++)
       size += (hg->vindex[select[i]+1] - hg->vindex[select[i]]); 
     size += (2 * ncandidates);    /* 2: append size of vtx and counts headers */
                    
     MPI_Allgather (&size, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);  

     /* setup arrays necessary for future MPI_Allgatherv */
     displs[0] = size = 0;
     for (i = 0; i < hgc->nProc_x; i++)
       size += each_size[i];                 /* compute total size of rbuffer */
     for (i = 1; i < hgc->nProc_x; i++)
       displs[i] = displs[i-1] + each_size[i-1];    /* message displacements */
        
     /* allocate send buffer */  
     if (each_size[hgc->myProc_x] > 0
      && !(buffer =(char*)ZOLTAN_MALLOC(each_size[hgc->myProc_x]*sizeof(int)))){
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
     }          
     
     /* allocate receive buffers */
     if (size > 0 && !(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))  {
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
     }          
             
     /* fill messages as list of <gno, gno's edge count, list of edge gno's> */
     ip = (int*) buffer;
     for (i = 0; i < ncandidates; i++)   {
       *ip++ = VTX_LNO_TO_GNO (hg, select[i]);                 /* vertex gno */
       *ip++ = hg->vindex[select[i]+1] - hg->vindex[select[i]];     /* count */
       for (j = hg->vindex[select[i]]; j < hg->vindex[select[i]+1]; j++)  
         *ip++ = hg->vedge[j];                                     /* edges */
     }        
          
     /* send ncandidates vertices & their edges to all row neighbors */
     MPI_Allgatherv (buffer, each_size[hgc->myProc_x], MPI_INT, rbuffer,
      each_size, displs, MPI_INT, hgc->row_comm);
               
     /************************ PHASE 2: ***************************************/    
         
     /* extract data from rbuffer, place in vindex & vedge arrays */ 
     ip = (int*) rbuffer;
     pins = 0;
     for (i = 0; i < total_count; i++)  {
       m_vindex[i] = pins;                
       m_gno   [i] = *ip++;
       edgecount   = *ip++;
       
       if (pincnt < pins + edgecount)  {
         /* need to realloc() m_vedge array */
         m_vedge = (int*)ZOLTAN_REALLOC(m_vedge, (pins+edgecount) * sizeof(int));
         if (!m_vedge)  {
           ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
           return ZOLTAN_MEMERR;
         }     
         pincnt += edgecount;  /* pincnt is "high water" size of m_vedge array */
       }
        
       while (edgecount-- > 0)
         m_vedge[pins++] = *ip++;    
     } 
     m_vindex[i] = pins;
          
     Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);            
               
     /* initialize array to process recv'd candidates in random order */
     for (i = 0; i < total_count; i++)
       c_order[i] = i;            /* at this point, visit in recv'd order */
     Zoltan_Rand_Perm_Int (c_order, total_count);
     
     /* for each candidate vertex, compute all local partial inner products */
     for (k = 0; k < total_count; k++)  {
       vertex = c_order[k];
       for (i = 0; i < hg->nVtx; i++)
         tsums[i] = psums[i] = 0.0;

       if (hg->ewgt == NULL)
         for (i = m_vindex[vertex]; i < m_vindex[vertex+1]; i++)
           for (j = hg->hindex [m_vedge[i]]; j < hg->hindex [m_vedge[i]+1]; j++)
             ++psums [hg->hvertex[j]];
       else
         for (i = m_vindex[vertex]; i < m_vindex[vertex+1]; i++)
           for (j = hg->hindex [m_vedge[i]]; j < hg->hindex [m_vedge[i]+1]; j++)
             psums [hg->hvertex[j]] += hg->ewgt [hg->hvertex[j]];
                           
       /* if also a local vtx, remove self inner product (false maximum) */
       if (VTX_TO_PROC_X (hg, m_gno[vertex]) == hgc->myProc_x)
         psums [VTX_GNO_TO_LNO (hg, m_gno[vertex])] = 0;
                  
       /* Want to use sparse communication with explicit summing later but
          for now, all procs in my column have same complete inner products */      
       MPI_Allreduce(psums, tsums, hg->nVtx, MPI_FLOAT, MPI_SUM, hgc->col_comm);
          
       /* each proc computes best, all rows in a column compute same answer */        
       m_bestsum [vertex] = -1;
       m_bestv   [vertex] =  0;
       for (i = 0; i < hg->nVtx; i++)
         if (tsums[i] > m_bestsum[vertex]  &&  cmatch[i] == i)  {
           m_bestsum[vertex] = tsums[i];
           m_bestv  [vertex] = i;
         }    
       cmatch [m_bestv[vertex]] = -1; /* mark pending match to avoid conflicts */
     }
        
     /************************ PHASE 3: **************************************/

     /* allocate buffers to send back best results information */
     size = PHASE3 * total_count;     
     if (size > 0 && (
         !(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int) * hgc->nProc_x))
      || !(buffer  = (char*) ZOLTAN_MALLOC (size * sizeof(int)))))    {
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
     }  
                      
     /* fill send buffer with candidate's id, best match gno, and best i.p. */       
     ip = (int*) buffer; 
     for (vertex = 0; vertex < total_count; vertex++)  {
       *ip++ = m_gno [vertex];                         /* candidate's gno */
       *ip++ = VTX_LNO_TO_GNO (hg, m_bestv [vertex]);  /* best match gno */   
       *ip++ = m_bestsum [vertex];                     /* best i.p. value */
     }
      
     /* send/rec best results information to all columns in my row */     
     MPI_Allgather(buffer, size, MPI_INT, rbuffer, size, MPI_INT, hgc->row_comm);
             
     /* initialize the array to hold the best reported i.p. from each proc */
     for (i = 0; i < hg->nVtx; i++)
       b_bestsum[i] = -1.0;                /* any negative value will work */
                 
     /* read each message (candidate id, best match id, and best i.p.) */  
     ip = (int*) rbuffer;
     for (i = 0; i < total_count * hgc->nProc_x; i++)   {
       vertex  = *ip++;
       bestv   = *ip++;
       bestsum = *ip++;
       
       /* ignore messages whose candidate I don't own */
       if (VTX_TO_PROC_X (hg, vertex) != hgc->myProc_x)
          continue;
           
       /* if I own the candidate, replace current i.p. with better */
       /* Note: ties are broken to favor local vs. global matches */   
       lno =  VTX_GNO_TO_LNO (hg, vertex);  
       if ((bestsum > b_bestsum [lno]) || (bestsum == b_bestsum[lno]
        && VTX_TO_PROC_X (hg, b_gno[lno]) != hgc->myProc_x
        && VTX_TO_PROC_X (hg, bestv)      == hgc->myProc_x))    {        
           b_gno     [lno] = bestv;
           b_bestsum [lno] = bestsum;
       }                   
     }   
                               
     Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);
     
     /************************ PHASE 4: ***************************************/
     
     /* need to tell everyone, which matches I accepted/rejected */
     /* first prepare to send matches back using MPI_ALLGATHERV */
     displs[0] = 0;
     for (i = 0; i < hgc->nProc_x; i++)
        each_size[i] = PHASE4 * each_count[i];
     for (i = 1; i < hgc->nProc_x; i++)
        displs[i] = displs[i-1] + each_size[i-1];
          
     /* allocate message send/rec buffers */   
     if (!(buffer  = (char*) ZOLTAN_MALLOC(PHASE4 * ncandidates * sizeof(int)))
      || !(rbuffer = (char*) ZOLTAN_MALLOC(PHASE4 * total_count * sizeof(int)))){
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
     }  
        
     /* fill send buffer with messages. A message is two matched gno's */
     /* Note: match to self if inner product is below threshold */                     
     ip = (int*) buffer; 
     for (i = 0; i < ncandidates; i++)   {
       *ip++ = VTX_LNO_TO_GNO (hg, select[i]);
       *ip++ = (b_bestsum [select[i]] > THRESHOLD) 
         ? b_gno[select[i]] : VTX_LNO_TO_GNO (hg, select[i]); 
     }
           
     MPI_Allgatherv(buffer, PHASE4 * ncandidates, MPI_INT, rbuffer, each_size,
      displs, MPI_INT, hgc->row_comm);
     
     /* extract messages. If I own either gno, process it. */
     /* Note: -gno-1 designates an external match */ 
     ip = (int*) rbuffer;                   
     for (i = 0; i < total_count; i++)  {
       bestv  = *ip++;
       vertex = *ip++; 
                
       if (VTX_TO_PROC_X (hg, bestv)  == hgc->myProc_x
        && VTX_TO_PROC_X (hg, vertex) != hgc->myProc_x)                       
           match [VTX_GNO_TO_LNO (hg, bestv)] = -vertex-1;
             
       if (VTX_TO_PROC_X (hg, vertex) == hgc->myProc_x
        && VTX_TO_PROC_X (hg, bestv)  != hgc->myProc_x)                
           match [VTX_GNO_TO_LNO (hg, vertex)] = -bestv-1;
             
       if (VTX_TO_PROC_X (hg, bestv)  == hgc->myProc_x
        && VTX_TO_PROC_X (hg, vertex) == hgc->myProc_x)   {
           int v1 = VTX_GNO_TO_LNO (hg, bestv);             
           int v2 = VTX_GNO_TO_LNO (hg, vertex);                
           match [v1] = v2;
           match [v2] = v1;
       }                   
     }       
     Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);                       
  } /* DONE: end of large loop over rounds */
        
  Zoltan_Multifree (__FILE__, __LINE__, 7, &psums, &tsums, &select, &cmatch,
   &each_size, &each_count, &displs); 
  Zoltan_Multifree (__FILE__, __LINE__, 9, &m_vindex, &m_vedge, &m_gno, &b_gno,
   &m_bestsum, &m_bestv, &b_bestsum, &order, &c_order);     
  return ZOLTAN_OK;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
