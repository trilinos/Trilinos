/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
<<<<<<< phg_match.c
 *    $Author$
 *    $Date$
 *    $Revision$
=======
 *    $Author$
 *    $Date$
 *    $Revision$
>>>>>>> 1.39
 ****************************************************************************/

 
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdlib.h>
#include "phg.h"


static ZOLTAN_PHG_MATCHING_FN matching_local;    /* dummy function for local matching */
static ZOLTAN_PHG_MATCHING_FN matching_ipm;      /* inner product matching */
static ZOLTAN_PHG_MATCHING_FN matching_col_ipm;  /* local ipm along proc columns*/


/*****************************************************************************/
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams *hgp)
{
    int exist=1;
    
    if (!strcasecmp(hgp->redm_str, "no"))        hgp->matching = NULL;
    else if (!strncasecmp(hgp->redm_str, "l-", 2))  {
        HGPartParams hp;

        strcpy(hp.redm_str, hgp->redm_str+2);
        strcpy(hp.redmo_str, hgp->redmo_str);
        if (!Zoltan_HG_Set_Matching_Fn(&hp)) {
            exist = 0;
            hgp->matching = NULL;
        } else {            
            hgp->matching = matching_local; /* just to make sure that coarsening will
                                            continue.
                                            We'll not call this code for global matching.
                                            Actually, we'll pick the best local, but
                                            code structure doesn't allow us to use a
                                            function */
            hgp->locmatching = hp.matching;
            hgp->matching_opt = hp.matching_opt;
        }
    } else if (!strcasecmp(hgp->redm_str, "c-ipm"))  hgp->matching = matching_col_ipm;    
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
  HGraph *hg,
  Matching match,
  PHGPartParams *hgp)
{
float *old_ewgt = NULL, *new_ewgt = NULL;
int   err = ZOLTAN_OK;
char  *yo = "Zoltan_PHG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hg->ewgt && hgp->ews) {
     if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto End;
     }
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
  if (hg->ewgt && hgp->ews)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT (zz, yo);
  return err;
}


static int matching_local(ZZ *zz, HGraph *hg, Matching match)
{
    uprintf(hg->comm, " Something wrong! This function should not be called at all!\n");
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
    int   i, j, n, v1, v2, edge, maxip, maxindex;
    int   matchcount=0;
    int   *ips, *gips, *adj;
    char  *yo = "matching_col_ipm";
    PHGComm *hgc = hg->comm;  

    if (!(ips = (int*) ZOLTAN_MALLOC(hg->nVtx * sizeof(int))) 
     || !(gips = (int*) ZOLTAN_MALLOC(hg->nVtx * sizeof(int)))
     || !(adj = (int*) ZOLTAN_MALLOC(hg->nVtx * sizeof(int)))) {
        Zoltan_Multifree(__FILE__, __LINE__, 3, &ips, &gips, &adj);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }
    
    for (i = 0; i < hg->nVtx; i++)
        ips[i] = 0;
        
    /* for every vertex */
    for (v1 = 0; v1 < hg->nVtx; v1++) {
        if (match[v1] != v1)
            continue;
        
        n = 0; /* number of neighbors */
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
                if (!ips[v2]++) /* TODO: use edge weight */
                    adj[n++] = v2;
            }
        }

        /* sum up local inner products along proc column */
        /* for now, sum up i.p. value for all vertices; this is slow! */
        /* to do: 1) ignore vertices already matched 
                  2) use "superrows" with only nonzero values */
        MPI_Allreduce(ips, gips, hg->nVtx, MPI_INT, MPI_SUM, hgc->col_comm);              
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
            ips[v2] = 0;
        }
        if (maxindex != -1) {
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

    Zoltan_Multifree(__FILE__, __LINE__, 3, &ips, &gips, &adj);
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
             
static int matching_ipm (ZZ *zz, HGraph *hg, Matching match)
{
  int i, j, lno, loop, vertex, *psums, *tsums, *order;
  int count, size, *ip, bestv, bestsum, edgecount, pins, *cmatch;
  int ncandidates, nloop;
  int *select, pselect;
  int *m_gno;
  int *m_vindex, *m_vedge;  /* zero-based loopup of edges by vertex */
  int *m_bestsum, *m_bestv; /* column's best results for each matched vertex */
  int *b_gno, *b_bestsum;
  char *buffer, *rbuffer;    /* send and rec buffers */
  int *displs, *each_size, *each_count, total_count;
  PHGComm *hgc = hg->comm;  
  char  *yo = "matching_ipm";
    
  nloop = hg->dist_x[hgc->nProc_x] > 100 ? 15 : 2;
  ncandidates = MAX (1, hg->nVtx / (2 * nloop)) ; /* match impacts 2 vertices */
       
  psums = tsums = select = cmatch = each_size = displs = NULL;
  if (hg->nVtx > 0 && (
      !(order     = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))  
   || !(psums     = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
   || !(tsums     = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
   || !(cmatch    = (int*) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))))     {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
  }
  
  /* match[] is local slice of global matching array.  It uses local numbering 
   * (zero-based) initially, match[i] = i.  After matching, match[i]=i indicates
   * an unmatched vertex. A matching between vertices i & j is indicated by 
   * match[i] = j & match [j] = i.  NOTE: a match to an off processor vertex is
   * indicated my a negative number, -(gno+1), which must use global numbers
   * (gno's).        */
  /* order[] is used to impliemnet alternative vertex selection algorithms:
   * natural (lno), random, weight order, vertex size, etc. */
  for (i = 0; i < hg->nVtx; i++)
     order[i] = match[i] = i;
     


  if (hgc->nProc_x > 0 && (
      !(select    = (int*) ZOLTAN_MALLOC (ncandidates  * sizeof(int)))
   || !(each_size = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))
   || !(each_count= (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))         
   || !(displs    = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))))  {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
  }
  
  /* need to determine global total of ncandidates to allocate storage */
  MPI_Allgather (&ncandidates, 1, MPI_INT, each_count, 1, MPI_INT, hgc->row_comm);
     
  total_count= 0;
  for (i = 0; i < hgc->nProc_x; i++)
    total_count += each_count[i];
  m_vedge = m_vindex = m_gno = m_bestsum = m_bestv = b_gno = b_bestsum = 0;
  if (hg->nPins > 0 && 
   !(m_vedge   = (int*) ZOLTAN_MALLOC  (hg->nPins * 3 * sizeof(int))))  {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  if (total_count > 0 &&  (
      !(m_vindex = (int*) ZOLTAN_MALLOC ((total_count+1) * sizeof(int)))
   || !(m_gno    = (int*) ZOLTAN_MALLOC  (total_count    * sizeof(int)))
   || !(m_bestsum= (int*) ZOLTAN_MALLOC  (total_count    * sizeof(int)))
   || !(m_bestv  = (int*) ZOLTAN_MALLOC  (total_count    * sizeof(int)))))  {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
     
  if (hg->nVtx > 0 && (
      !(b_gno     = (int*) ZOLTAN_MALLOC  (hg->nVtx * sizeof(int)))
   || !(b_bestsum = (int*) ZOLTAN_MALLOC  (hg->nVtx * sizeof(int)))))  {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
     
                        
  /* Loop processing ncandidates vertices per column each pass. Each loop has 4 phases:
   * Phase 1: send ncandidates vertices for global matching - horizontal communication
   * Phase 2: sum  inner products, find best        - vertical communication
   * Phase 3: return best sums to owning column     - horizontal communication
   * Phase 4: return actual match selections        - horizontal communication */
  pselect = 0;                    /* marks position in vertices to be selected */
  for (loop = 0; loop < nloop; loop++)  {
     
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
     size += (2 * ncandidates);       /* append size of vtx and counts headers */
                    
     MPI_Allgather (&size, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);  

     displs[0] = size = 0;
     for (i = 0; i < hgc->nProc_x; i++)
       size += each_size[i];                 /* compute total size of rbuffer */
     for (i = 1; i < hgc->nProc_x; i++)
       displs[i] = displs[i-1] + each_size[i-1];    /* message displacements */
        
     if (each_size[hgc->myProc_x] > 0
      && !(buffer =(char*)ZOLTAN_MALLOC(each_size[hgc->myProc_x]*sizeof(int)))){
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
     }          
     if (size > 0 && !(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))  {
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
     }          
             
     /* Message is list of <gno, gno's edge count, list of edge gno's> */
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
       while (edgecount-- > 0)
         m_vedge[pins++] = *ip++;           
     } 
     m_vindex[i] = pins;
          
     Zoltan_Multifree (__FILE__, __LINE__, 2, &buffer, &rbuffer);            
               
     /* for each match vertex, compute all local partial inner products */
     for (vertex = 0; vertex < total_count; vertex++)  {
       for (i = 0; i < hg->nVtx; i++)
         tsums[i] = psums[i] = 0;

       for (i = m_vindex[vertex]; i < m_vindex[vertex+1]; i++)
         for (j = hg->hindex [m_vedge[i]]; j < hg->hindex [m_vedge[i]+1]; j++)
           ++psums [hg->hvertex[j]];                       /* unweighted??? */
              
       /* if local vtx, remove self inner product which is a false maximum */
       if (VTX_TO_PROC_X (hg, m_gno[vertex]) == hgc->myProc_x)
         psums [VTX_GNO_TO_LNO (hg, m_gno[vertex])] = 0;
       
       /* Want to use sparse communication with explicit summing later but
          for now, all procs in my column have same complete inner products */      
       MPI_Allreduce(psums, tsums, hg->nVtx, MPI_INT, MPI_SUM, hgc->col_comm);             
                            
       /* each proc computes best, all rows in a column compute same answer */        
       m_bestsum [vertex] = -1;
       m_bestv   [vertex] =  0;
       for (i = 0; i < hg->nVtx; i++)
         if (tsums[i] > m_bestsum[vertex]  &&  cmatch[i] == i)  {
           m_bestsum [vertex] = tsums[i];
           m_bestv   [vertex] = i;
         }    
       cmatch [m_bestv[vertex]] = -1;      /* pending match */       
     }
        
     /************************ PHASE 3: **************************************/

     size = 3 * total_count;     
     if (size > 0 && (
         !(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int) * hgc->nProc_x))
      || !(buffer  = (char*) ZOLTAN_MALLOC (size * sizeof(int)))))    {
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
     }  
                      
     /* prepare send buffer */       
     ip = (int*) buffer; 
     for (vertex = 0; vertex < total_count; vertex++)  {
       *ip++ = m_gno [vertex];
       *ip++ = VTX_LNO_TO_GNO (hg, m_bestv [vertex]);
       *ip++ = m_bestsum [vertex];
     }
      
     /* send/rec to all columns in my row */     
     MPI_Allgather(buffer, size, MPI_INT, rbuffer, size, MPI_INT, hgc->row_comm);
             
     for (i = 0; i < hg->nVtx; i++)
       b_bestsum[i] = -2;   
                 
     ip = (int*) rbuffer;
     for (i = 0; i < total_count * hgc->nProc_x; i++)   {
       vertex  = *ip++;
       bestv   = *ip++;
       bestsum = *ip++;
       
       if (VTX_TO_PROC_X (hg, vertex) != hgc->myProc_x)
          continue;    
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
     displs[0] = 0;
     for (i = 0; i < hgc->nProc_x; i++)
        each_size[i] =  2 * each_count[i];
     for (i = 1; i < hgc->nProc_x; i++)
        displs[i] = displs[i-1] + each_size[i-1];
          
     if (!(buffer  = (char*) ZOLTAN_MALLOC (2 * ncandidates * sizeof(int)))
      || !(rbuffer = (char*) ZOLTAN_MALLOC (2 * total_count * sizeof(int))))  {
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
     }  
                             
     ip = (int*) buffer; 
     for (i = 0; i < ncandidates; i++)   {
       *ip++ = VTX_LNO_TO_GNO (hg, select[i]);
       *ip++ = (b_bestsum [select[i]] > 0) 
         ? b_gno[select[i]] : VTX_LNO_TO_GNO (hg, select[i]);
     }
           
     MPI_Allgatherv(buffer, 2 * ncandidates, MPI_INT, rbuffer, each_size, displs,
      MPI_INT, hgc->row_comm);
     
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
  } /* DONE: end of large loop over LOOP */
        
  Zoltan_Multifree (__FILE__, __LINE__, 7, &psums, &tsums, &select, &cmatch,
   &each_size, &each_count, &displs); 
  Zoltan_Multifree (__FILE__, __LINE__, 7, &m_vindex, &m_vedge, &m_gno, &b_gno,
   &m_bestsum, &m_bestv, &b_bestsum);     
  return ZOLTAN_OK;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
