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
  if (hgp->edge_scaling) {
     if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto End;
     }
 
     Zoltan_PHG_Scale_Edges (zz, hg, new_ewgt, hgp);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
  }

  /* Create/update scale vector for vertices for inner product */
  if (hgp->vtx_scaling) {
     if (hgp->vtx_scal==NULL){  /* first level in V-cycle */
        if (!(hgp->vtx_scal = (float*) ZOLTAN_MALLOC (hg->nVtx * 
                               sizeof(float)))) {
           ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
           err = ZOLTAN_MEMERR;
           goto End;
        }
     }
 
     Zoltan_PHG_Scale_Vtx (zz, hg, hgp);
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
     err = hgp->matching (zz, hg, match, hgp);

End: 

  /* Restore the old edge weights if scaling was used. */
  if (hgp->edge_scaling)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT (zz, yo);
  return err;
}


static int matching_local(ZZ *zz, HGraph *hg, Matching match, PHGPartParams *hgp)
{
    uprintf(hg->comm, "Something wrong! This function should not be called!\n");
    /* UVC: NOTE:
       The reason that we're not doing local matchin in this function, we don't
       have access to parameter structure. So there is no way to figure
       out which "local" matching needs to be called. Hence we do it
       in Zoltan_PHG_Matching */
    /* EBEB: Added the parameter struct as an input argument, so now we 
       could change the structure.  */
    return ZOLTAN_OK;
}


    
/* local inner product matching among vertices in each proc column */
/* code adapted from serial matching_ipm method */
#define MAX_NNZ 50  /* Max number of nonzeros to store for each inner product */
static int matching_col_ipm(ZZ *zz, HGraph *hg, Matching match, PHGPartParams *hgp)
{
    int   i, j, k, v1, v2, edge, best_vertex;
    int   nadj, dense_comm;
    float maxip= 0.;
    int   *adj=NULL;
    int   *order=NULL;
    float *lips, *gips; /* local and global inner products */
    float *ptr;
    char  *sendbuf, *recvbuf; /* comm buffers */
    char  *yo = "matching_col_ipm";
    PHGComm *hgc = hg->comm;  
    float lquality[3] = {0,0,0}; /* local  matchcount, matchweight */
    float gquality[3] = {0,0,0}; /* global matchcount, matchweight */

    lips = gips = NULL;
    sendbuf = recvbuf = NULL;

    if (!(lips = (float*) ZOLTAN_MALLOC(hg->nVtx * sizeof(float))) 
     || !(gips = (float*) ZOLTAN_MALLOC(hg->nVtx * sizeof(float)))
     || !(adj =  (int*)  ZOLTAN_MALLOC(hg->nVtx * sizeof(int))) 
     || !(order = (int*)  ZOLTAN_MALLOC(hg->nVtx * sizeof(int)))){
        Zoltan_Multifree(__FILE__, __LINE__, 4, &lips, &gips, &adj, &order);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }
    
    for (i = 0; i < hg->nVtx; i++){
        lips[i] = gips[i] = .0;
    }

    /* Do dense communication for small problems. */
    dense_comm = (hg->nVtx < 2*MAX_NNZ);
    /* dense_comm = 1;  EBEB */

    if (!dense_comm){

      /* allocate buffers */
      if (!(sendbuf = (char*) ZOLTAN_MALLOC(2*MAX_NNZ * sizeof(float)))) {
        ZOLTAN_FREE(&sendbuf);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
      }
      if (hgc->myProc_y==0 &&
        !(recvbuf = (char*) ZOLTAN_MALLOC(2*MAX_NNZ*hgc->nProc_y* sizeof(float)))) {
        ZOLTAN_FREE(&sendbuf);
        ZOLTAN_FREE(&recvbuf);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
      }
    }

    /* Compute vertex visit order. */
    Zoltan_PHG_Vertex_Visit_Order(zz, hg, hgp, order);

    /* for every vertex */
    for (k = 0; k < hg->nVtx; k++) {
        v1 = order[k];

        if (match[v1] != v1)
            continue;

        nadj = 0;  /* number of neighbors */

        /* for every hyperedge containing the vertex */
        for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
            edge = hg->vedge[i];
                
            /* for every other vertex in the hyperedge */
            for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
                v2 = hg->hvertex[j];
                if (match[v2] == v2) {
                    /* v2 is not matched yet */
                    if (lips[v2]==0.0)   /* v2 is a new neighbor */
                        adj[nadj++] = v2;
                    lips[v2] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
                }
            }
        }

        /* sum up local inner products along proc column */
        /* currently we communicate for just one vertex at a time */
        /* future: use sparse communication for a chunk of vertices */

#ifdef DEBUG_EB
        printf("[%1d] Debug: %d out of %d entries are nonzero\n", 
          zz->Proc, nadj, hg->nVtx);
#endif

        if (dense_comm){
          /* dense communication */
          MPI_Reduce(lips, gips, hg->nVtx, MPI_FLOAT, MPI_SUM, 0,
            hgc->col_comm);              
          /* clear lips array for next iter  */
          for (i = 0; i < hg->nVtx; i++) {
            lips[i] = 0.0;
          }
        }
        else {
          /* sparse communication */
          /* we send partial inner products to root row */

          /* pack data into send buffer */
          ptr = (float *) sendbuf;
#ifdef DEBUG_EB
          printf("Debug: Send data; myProc_y=%d, v1=%d\n", hgc->myProc_y, v1);
#endif
          if (nadj <= MAX_NNZ){
            for (i=0; i<MIN(nadj,MAX_NNZ); i++){
              *ptr++ = (float) adj[i];
              *ptr++ = lips[adj[i]];
#ifdef DEBUG_EB
              printf("%d %4.2f, ", adj[i], *(ptr-1));
#endif
            }
            if (i<MAX_NNZ-1){
              /* marker to say there is no more data */
              *ptr++ = -1.;
              *ptr++ = -1.;
#ifdef DEBUG_EB
              printf("%d %d, ", -1, -1);
#endif
            }
#ifdef DEBUG_EB
            printf("\n");
#endif
          }
          else {
            /* pick highest values if too many nonzeros */
            /* naive algorithm to find top MAX_NNZ values */ 
#ifdef DEBUG_EB
            printf("Debug: nadj= %d > MAX_NNZ= %d, inexact inner product!\n",
                nadj, MAX_NNZ);
#endif
            for (i=0; i<MAX_NNZ; i++){
              maxip = 0.0;
              for (j=0; j<nadj; j++){
                if (lips[adj[j]]>maxip){
                  best_vertex = adj[j];
                  maxip = lips[best_vertex];
                  lips[best_vertex] = 0.0;
                }
              }
              *ptr++ = (float) best_vertex;
              *ptr++ = maxip;
            }
          }

          /* send partial inner product values to root row */
          /* use fixed size, probably faster than variable sized Gatherv */
          MPI_Gather (sendbuf, 2*MAX_NNZ, MPI_FLOAT, recvbuf,
           2*MAX_NNZ, MPI_FLOAT, 0, hgc->col_comm);

          /* root unpacks data into gips array */
          if (hgc->myProc_y==0){
            nadj = 0;
            for (i=0; i<hgc->nProc_y; i++){
#ifdef DEBUG_EB
              printf("Debug: Received data, v1=%d, i=%d\n", v1, i);
#endif
              ptr = (float *) recvbuf;
              ptr += i*2*MAX_NNZ;
              for (j=0; j<MAX_NNZ; j++){
                v2 = *ptr++;
                if (v2<0) break; /* skip to data from next proc */
#ifdef DEBUG_EB
                printf(" %d %4.2f, ", v2, *ptr);
#endif
                if (gips[v2]==0.0)
                  adj[nadj++] = v2;
                gips[v2] += *ptr++;
              }
#ifdef DEBUG_EB
              printf("\n");
#endif
            }
          }
        }

        /* now choose the vector with greatest inner product */
        /* do only on root proc, then broadcast to column */
        if (hgc->myProc_y==0){
          maxip = 0.;
          best_vertex = -1;
          if (dense_comm){
#ifdef DEBUG_EB
            printf("Debug: v1=%d, gips=\n", v1);
#endif
            for (i = 0; i < hg->nVtx; i++) {
              v2 = order[i];
#ifdef DEBUG_EB
              if (gips[v2]>0) printf(" %d=%4.2f, ", v2, gips[v2]);
#endif
              if (gips[v2] > maxip && v2 != v1 && match[v2] == v2) {
                  maxip = gips[v2];
                  best_vertex = v2;
              }
            }
#ifdef DEBUG_EB
            printf("\n");
#endif
          }
          else { /* sparse */
#ifdef DEBUG_EB
            printf("Debug: v1=%d, gips=\n", v1);
#endif
            for (i = 0; i < nadj; i++) {
              v2 = adj[i];
#ifdef DEBUG_EB
              if (gips[v2]>0) printf(" %d=%4.2f, ", v2, gips[v2]);
#endif
              /* Tie-breaking strategy makes a difference!
               * Random or using the order array seem to work well. */
              if (((gips[v2] > maxip) || (gips[v2]==maxip && 
                order[v2] < order[best_vertex]))
                && v2 != v1 && match[v2] == v2) {
                  maxip = gips[v2];
                  best_vertex = v2;
              }
              lips[v2] = gips[v2] = .0; /* clear arrays for next iter */
            }
#ifdef DEBUG_EB
            printf("\n");
#endif
          }
        } 

        /* broadcast the winner, best_vertex */
        MPI_Bcast(&best_vertex, 1, MPI_INT, 0, hgc->col_comm); 

        /* match if inner product > 0 */
        if (best_vertex > -1) {
            match[v1] = best_vertex;
            match[best_vertex] = v1;
            if (hgc->myProc_y==0){
              /* match quality computation */
              lquality[0] += maxip;
              lquality[1] += 1.0;
            }
        } 

        
    }

    if (hgc->myProc_y==0){
      lquality[2] = hg->nVtx; /* to find global number of vertices */
      MPI_Allreduce(lquality, gquality, 3, MPI_FLOAT, MPI_SUM, hgc->row_comm); 
       
      uprintf (hgc, "LOCAL (GLOBAL) i.p. sum %.2f (%.2f), matched pairs %d (%d), "
       "total vertices %d\n", lquality[0], gquality[0], (int)lquality[1],
       (int)gquality[1], (int)gquality[2]);  
    }

    /*
    printf("DEBUG: Final Matching:\n");
    for(i = 0; i < hg->nVtx; i++)
        printf("%d, %d\n", i, match[i]);
    printf("\n");
    */

    if (!dense_comm){
      ZOLTAN_FREE(&sendbuf);
      ZOLTAN_FREE(&recvbuf);
    }
    Zoltan_Multifree(__FILE__, __LINE__, 4, &lips, &gips, &adj, &order);
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
#define THRESHOLD 0.0 /* ignore inner products (i.p.) less than threshold */
       
static int matching_ipm (ZZ *zz, HGraph *hg, Matching match, PHGPartParams *hgp)
{
  int i, j, k, lno, loop, vertex, *order=NULL;
  int *c_order=NULL, pincnt, count, size, *ip, bestv, edgecount, pins;
  int *cmatch=NULL, ncandidates, nrounds, *select=NULL, pselect;
  int *m_gno=NULL;
  int *m_vindex=NULL, *m_vedge=NULL;  /* zero-based loopup of edges by vertex */
  int *m_bestv=NULL;                 /* col's best results for matched vertex */
  int *b_gno=NULL;
  float *psums=NULL, *tsums=NULL, *m_bestsum=NULL, *b_bestsum=NULL, bestsum;
  char *buffer=NULL, *rbuffer=NULL;    /* send and rec buffers */
  int *displs=NULL, *each_size=NULL, *each_count=NULL, total_count;
  float local_quality[3] = {0,0,0}, global_quality[3] = {0,0,0},  gno;  
  PHGComm *hgc = hg->comm;  
  char  *yo = "matching_ipm";
int bobcount=0;


  /* determine number of basic matching rounds scaled by number of procs */
  nrounds = hgc->nProc_x * LOOP_FACTOR;
  ncandidates = hg->nVtx/(2 * nrounds) + 1;  /* 2: each match pairs 2 vertices */
  
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
   * (gno's) and not local numbers, lno.        */
  for (i = 0; i < hg->nVtx; i++)
    match[i] = i;                /* initialize match array to all unmatched */
   
  /* order[] is used to implement alternative vertex visiting algorithms:
   * natural (lno), random, weight order, vertex size, etc. */

  /* Compute vertex visit order. Random is default. */
  Zoltan_PHG_Vertex_Visit_Order(zz, hg, hgp, order);
  
  if (hgc->nProc_x > 0 && (
      !(select     = (int*) ZOLTAN_MALLOC (ncandidates  * sizeof(int)))
   || !(each_size  = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))
   || !(each_count = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))         
   || !(displs     = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))))  {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
  }
  
  /* need to determine global total of ncandidates to allocate storage */
  MPI_Allreduce (&ncandidates, &total_count, 1, MPI_INT, MPI_SUM, hgc->row_comm);
  
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
     
    /* initialize temporary copy of the match array, cmatch  */
    for (i = 0; i < hg->nVtx; i++)
       cmatch[i] = match[i];   
            
     /* Select next ncandidates unmatched vertices to globally match. */                                            
     for (count = 0; count < ncandidates && pselect < hg->nVtx; pselect++)
       if (cmatch[order[pselect]] == order[pselect])  {  /* unmatched */
         select[count++] = order[pselect];       /* select it */       
         cmatch[order[pselect]] = -1;            /* mark it as a pending match */
       }  
       
     MPI_Allgather (&count, 1, MPI_INT, each_count, 1, MPI_INT, hgc->row_comm);
     
     total_count = 0;
     for (i = 0; i < hgc->nProc_x; i++)
        total_count += each_count[i];
     if (total_count == 0)
        break;
                
     /* determine the size of the send buffer & allocate it */   
     size = 2 * count;     /* per vertex: vertex id and count of edges */
     for (i = 0; i < count; i++)
       size += (hg->vindex[select[i]+1] - hg->vindex[select[i]]); 
                    
     MPI_Allgather (&size, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);  

     /* setup arrays necessary for future MPI_Allgatherv */
     displs[0] = size = 0;
     for (i = 0; i < hgc->nProc_x; i++)
       size += each_size[i];                 /* compute total size of rbuffer */
     for (i = 1; i < hgc->nProc_x; i++)
       displs[i] = displs[i-1] + each_size[i-1];    /* message displacements */
        
     /* allocate send buffer */  
     if (each_size[hgc->myProc_x] > 0
      && !(buffer =(char*)ZOLTAN_MALLOC(1+each_size[hgc->myProc_x]*sizeof(int)))){
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
     }          
     
     /* allocate receive buffers */
     if (size > 0 && !(rbuffer = (char*) ZOLTAN_MALLOC (1+size * sizeof(int))))  {
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
     }          
             
     /* fill messages as list of <gno, gno's edge count, list of edge gno's> */
     ip = (int*) buffer;
     for (i = 0; i < count; i++)   {
       *ip++ = VTX_LNO_TO_GNO (hg, select[i]);                 /* vertex gno */
       *ip++ = hg->vindex[select[i]+1] - hg->vindex[select[i]];     /* count */
       for (j = hg->vindex[select[i]]; j < hg->vindex[select[i]+1]; j++)  
         *ip++ = hg->vedge[j];                                     /* edges */
     }        
          
     /* send (by size) count vertices & their edges to all row neighbors */
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
             psums [hg->hvertex[j]] += 1.0;
       else
         for (i = m_vindex[vertex]; i < m_vindex[vertex+1]; i++)
           for (j = hg->hindex [m_vedge[i]]; j < hg->hindex [m_vedge[i]+1]; j++)
             psums [hg->hvertex[j]] += hg->ewgt [m_vedge[i]];
                           
       /* if local vtx, remove self inner product (false maximum) and matched */
       if (VTX_TO_PROC_X (hg, m_gno[vertex]) == hgc->myProc_x)
         psums [VTX_GNO_TO_LNO (hg, m_gno[vertex])] = 0;
       for (i = 0; i < hg->nVtx; i++)
         if (match[i] != i)
           psums[i] = 0.0;                  
     
       /* Want to use sparse communication with explicit summing later but
          for now, all procs in my column have same complete inner products */      
       MPI_Allreduce(psums, tsums, hg->nVtx, MPI_FLOAT, MPI_SUM, hgc->col_comm);
          
       /* each proc computes best, all rows in a column compute same answer */        
       m_bestsum [vertex] = -1.0;             /* any negative float will do */
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
     if (
         !(rbuffer = (char*) ZOLTAN_MALLOC (1+size * sizeof(int) * hgc->nProc_x))
      || !(buffer  = (char*) ZOLTAN_MALLOC (1+size * sizeof(int))))    {
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
     }  
                      
     /* fill send buffer with candidate's id, best match gno, and best i.p. */       
     ip = (int*) buffer; 
     for (vertex = 0; vertex < total_count; vertex++)  {
       *ip++ = m_gno [vertex];                         /* candidate's gno */
       *ip++ = VTX_LNO_TO_GNO (hg, m_bestv [vertex]);  /* best match gno */   
       * (float*) ip++ = m_bestsum [vertex];           /* best i.p. value */
     }
      
     /* send/rec best results information to all columns in my row */     
     MPI_Allgather(buffer, size, MPI_INT, rbuffer, size, MPI_INT, hgc->row_comm);
             
     /* initialize the array to hold the best reported i.p. from each proc */
     for (i = 0; i < hg->nVtx; i++)
       b_bestsum[i] = -2.0;            /* any negative value < -1.0 will work */
                 
     /* read each message (candidate id, best match id, and best i.p.) */  
     for (ip = (int*) rbuffer; ip < ((int*) rbuffer) + size * hgc->nProc_x; ) {
       vertex  = *ip++;
       bestv   = *ip++;
       bestsum = * (float*)ip++;
               
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
     
     /* need to tell everyone, which matches I accepted (rejected implied) */
     /* first prepare to send matches back using MPI_ALLGATHERV */
     displs[0] = 0;
     for (i = 0; i < hgc->nProc_x; i++)
        each_size[i] = PHASE4 * each_count[i];
     for (i = 1; i < hgc->nProc_x; i++)
        displs[i] = displs[i-1] + each_size[i-1];
          
     /* allocate message send/rec buffers */   
     if (!(buffer  = (char*) ZOLTAN_MALLOC(1+PHASE4 * count * sizeof(int)))
      || !(rbuffer = (char*) ZOLTAN_MALLOC(1+PHASE4 * total_count * sizeof(int)))){
         ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
     }  
        
     /* fill send buffer with messages. A message is two matched gno's */
     /* Note: match to self if inner product is below threshold */                     
     ip = (int*) buffer; 
     for (i = 0; i < count; i++)   {
       *ip++ = gno = VTX_LNO_TO_GNO (hg, select[i]);
       *ip++ = (b_bestsum [select[i]] > THRESHOLD) 
         ? b_gno[select[i]] : VTX_LNO_TO_GNO (hg, select[i]); 
         
if (b_bestsum[select[i]] > THRESHOLD)  {
  local_quality[0] += b_bestsum[select[i]];  
  local_quality[1] += 1.0;
}       
     }
           
     MPI_Allgatherv(buffer, PHASE4 * count, MPI_INT, rbuffer, each_size,
      displs, MPI_INT, hgc->row_comm);
     
     /* extract messages. If I own either gno, process it. */
     /* Note: -gno-1 designates an external match as a negative number */ 
     for (ip = (int*) rbuffer; ip < (int*) rbuffer + (PHASE4 * total_count); ) {
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

/*   
for (i = 0; i < hg->nVtx; i++)
  uprintf (hgc, "BOBBOB match[%d] = %d\n", VTX_LNO_TO_GNO (hg, i),
   (match[i] < 0) ? -match[i]-1 : VTX_LNO_TO_GNO (hg, match[i]));
*/
    
bobcount = 0;  
for (i = pselect; i < hg->nVtx; i++)
  if (match[order[pselect]] == order[pselect])
    bobcount++;
uprintf (hgc, "RTHRTH %d of %d, missed %d vertices with ncandidates %d\n",
 pselect, hg->nVtx, bobcount, ncandidates);

  
  local_quality[2] = hg->nVtx;    /* to compute the global number of vertices */
  
  MPI_Allreduce(local_quality,global_quality,3, MPI_FLOAT, MPI_SUM, hgc->row_comm); 
       
  uprintf (hgc, "LOCAL (GLOBAL) i.p. sum %.2f (%.2f), matched pairs %d (%d), "
   "total vertices %d\n", local_quality[0], global_quality[0], (int)local_quality[1],
   (int)global_quality[1], (int)global_quality[2]);  
  
  Zoltan_Multifree (__FILE__, __LINE__, 7, &psums, &tsums, &select, &cmatch,
   &each_size, &each_count, &displs); 
  Zoltan_Multifree (__FILE__, __LINE__, 9, &m_vindex, &m_vedge, &m_gno, &b_gno,
   &m_bestsum, &m_bestv, &b_bestsum, &order, &c_order);     
  return ZOLTAN_OK;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
