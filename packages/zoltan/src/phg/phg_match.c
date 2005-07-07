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


static ZOLTAN_PHG_MATCHING_FN pmatching_local; /* function for local matching */
static ZOLTAN_PHG_MATCHING_FN pmatching_ipm;   /* inner product matching */
static ZOLTAN_PHG_MATCHING_FN pmatching_col_ipm;   /* OLD column ipm, will be phased out */

#ifdef ALT_IPM
static ZOLTAN_PHG_MATCHING_FN pmatching_alt_ipm;   /* alternating ipm */
#endif

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
            hgp->matching = pmatching_local; 
            hgp->locmatching = hp.matching;
            hgp->matching_opt = hp.matching_opt;
        }
    } else if (!strcasecmp(hgp->redm_str, "c-ipm"))
        hgp->matching = pmatching_ipm;   
    else if (!strcasecmp(hgp->redm_str, "col-ipm")) /* old c-ipm */
        hgp->matching = pmatching_col_ipm;          /* will be removed later */
    else if (!strcasecmp(hgp->redm_str, "ipm"))
        hgp->matching = pmatching_ipm;
#ifdef ALT_IPM
    else if (!strcasecmp(hgp->redm_str, "alt-ipm"))
        hgp->matching = pmatching_alt_ipm;
#endif
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
     if (hg->nEdge && !(new_ewgt = (float*) 
                      ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
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
     Zoltan_PHG_Scale_Vtx (zz, hg, hgp);
  }

  /* Do the matching */
  if (hgp->matching) 
     err = hgp->matching (zz, hg, match, hgp);

End: 

  /* Restore the old edge weights if scaling was used. */
  if (hgp->edge_scaling)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT (zz, yo);
  return err;
}


static int pmatching_local(
  ZZ *zz,
  HGraph *hg,
  Matching match,
  PHGPartParams *hgp
)
{
    int limit=hg->nVtx, err=ZOLTAN_OK;
    PHGComm *hgc=hg->comm;
    int root_matchcnt, root_rank;
    
    err = hgp->locmatching (zz, hg, match, &limit);
    
    /* Optimization */
    if (hgp->matching_opt) 
        err = hgp->matching_opt (zz, hg, match, &limit);
    
    /* find the index of the proc in column group with the best match
       (max #matches); it will be our root proc */
    Zoltan_PHG_Find_Root(hg->nVtx-limit, hgc->myProc_y, hgc->col_comm,
                         &root_matchcnt, &root_rank);
    
    MPI_Bcast(match, hg->nVtx, MPI_INT, root_rank, hgc->col_comm);
    
    return err;
}


    
/* local inner product matching among vertices in each proc column */
/* this is (usually) faster than the full ipm but quality may be worse. */

#define MAX_NNZ 50  /* Max number of nonzeros to store for each inner product */
                    /* Reduce this value to save memory and comm volume. */
static int pmatching_col_ipm(
  ZZ *zz,
  HGraph *hg,
  Matching match,
  PHGPartParams *hgp
)
{
    int   i, j, k, v1, v2, edge, best_vertex;
    int   nadj, dense_comm;
    float maxip= 0.;
    int   *adj=NULL;
    int   *order=NULL;
    float *lips, *gips; /* local and global inner products */
    float *ptr;
    char  *sendbuf, *recvbuf; /* comm buffers */
    char   msg[160];          /* error messages */
    char  *yo = "pmatching_col_ipm";
    PHGComm *hgc = hg->comm;  
    float lquality[3] = {0,0,0}; /* local  matchcount, matchweight */

    lips = gips = NULL;
    sendbuf = recvbuf = NULL;
    adj = order = NULL;

    if (hg->nVtx)  
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
                if (v2 > hg->nVtx){
                  sprintf(msg, "vertex %d > %d is out of range!\n", v2, hg->nVtx);
                  ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
                }
                if (match[v2] == v2) {
                    /* v2 is not matched yet */
                    if (lips[v2]==0.0)   /* v2 is a new neighbor */
                        adj[nadj++] = v2;
                    /* Check for vertex scaling. For efficiency we may 
                       have to move the 'if' test out of the inner loop. 
                       Only scale for v2 not v1 to save flops.              */ 
                    if (hgp->vtx_scal)
                      lips[v2] += hgp->vtx_scal[v2]*
                                  (hg->ewgt ? hg->ewgt[edge] : 1.0);
                    else
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
          /* dense communication: treat lips & gips as dense vectors */
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

          /* pack data into send buffer as pairs (index, value) */
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
            if (nadj < MAX_NNZ){
              /* insert marker to say there is no more data */
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
            /* ptr = sendbuf */
#ifdef DEBUG_EB
            printf("Debug: nadj= %d > MAX_NNZ= %d, inexact inner product!\n",
                nadj, MAX_NNZ);
#endif
            /* Pick random selection of vertices if more than MAX_NNZ. */
            Zoltan_Rand_Perm_Int(adj, nadj, &(hgc->RNGState_col));

            /* Make sure highest value is among the selected vtx. */
            maxip = 0.0; 
            best_vertex = 0;
            for (j=0; j<nadj; j++){
              if (lips[adj[j]]>maxip){
                best_vertex = j;
                maxip = lips[adj[j]];
              }
            }
            if (best_vertex >= MAX_NNZ){
              /* put best in front */
              adj[0] = adj[best_vertex];
            }

            /* Pack the MAX_NNZ first vertices into send buffer */
            for (i=0; i<MAX_NNZ; i++){
              *ptr++ = (float) adj[i];
              *ptr++ = lips[adj[i]];
            }
          }

#ifdef DEBUG_EB
          /* EBEB Sanity check for debugging */
          ptr = (float *) sendbuf;
          for (i=0; i<MAX_NNZ; i++, ptr+=2){
            if (*ptr < 0) break;
            if (*ptr > hg->nVtx){
              sprintf(msg, "vertex %f > %d is out of range!\n", *ptr, hg->nVtx);
              ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
            }
          }
#endif

          /* send partial inner product values to root row */
          /* use fixed size, probably faster than variable sized Gatherv */
          MPI_Gather (sendbuf, 2*MAX_NNZ, MPI_FLOAT, recvbuf,
           2*MAX_NNZ, MPI_FLOAT, 0, hgc->col_comm);

          /* root unpacks data into gips array */
          if (hgc->myProc_y==0){
            nadj = 0;
            for (i=0; i<hgc->nProc_y; i++){
#ifdef DEBUG_EB
              printf("Debug: Received data, v1=%d, col=%d\n", v1, i);
#endif
              ptr = (float *) recvbuf;
              ptr += i*2*MAX_NNZ;
              for (j=0; j<MAX_NNZ; j++){
                v2 = *ptr++;
                if (v2<0) break; /* skip to data from next proc */
                /* Sanity check for debugging */
                if (v2 > hg->nVtx){
                  sprintf(msg, "vertex %d > %d is out of range!\n", v2, hg->nVtx);
                  ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
                }
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

    if (hgp->output_level >= PHG_DEBUG_LIST && hgc->myProc_y==0){
      float gquality[3] = {0,0,0}; /* global matchcount, matchweight */
      lquality[2] = hg->nVtx; /* to find global number of vertices */
      MPI_Allreduce(lquality, gquality, 3, MPI_FLOAT, MPI_SUM, hgc->row_comm); 
      uprintf (hgc, 
        "LOCAL (GLOBAL) i.p. sum %.2f (%.2f), matched pairs %d (%d), "
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

#ifdef ALT_IPM
/**************************************************************************
  Alternating ipm method. Alternate between full ipm and c-ipm. 
 *************************************************************************/
#define DO_FULL_IPM 2         /* Do full ipm every 2 levels */
                              /* This could be a parameter. */
static int pmatching_alt_ipm(
  ZZ *zz,
  HGraph* hg,
  Matching match,
  PHGPartParams *hgp
)
{
  int i;
  char *str;
  static int level=0;

  ++level;  /* we don't have access to level data, so keep track this way */

  str = hgp->redm_str; /* save original parameter string */

  if (level%DO_FULL_IPM == 0) 
    hgp->redm_str = "ipm";  // Need strcpy!
  else
    hgp->redm_str = "c-ipm";

  i = pmatching_ipm(zz, hg, match, hgp);

  /* set redm parameter back to original (alt-ipm) */
  hgp->redm_str = str;
  
  return i;
}
#endif

/****************************************************************************
 * inner product matching
 * Parallelized version of the serial algorithm (see hg_match.c)
 * Based on conversations with Rob Bisseling by Aaron Becker, UIUC, summer 2004
 * completed by R. Heaphy
 */
               
#define ROUNDS_CONSTANT 8     /* forces candidates to have enough freedom */ 
#define IPM_TAG        28731  /* MPI message tag, arbitrary value */
#define HEADER_COUNT    4     /* number of integers in header for messages
                                 in fixed sized in Phase 2 send buffer  */
                                 
/* these thresholds need to become parameters in phg - later ??? */
#define PSUM_THRESHOLD 0.0    /* ignore inner products (i.p.) < threshold */
#define TSUM_THRESHOLD 0.0    /* ignore inner products (i.p.) < threshold */
                                 
/* This routine encapsulates the common calls to use the Zoltan unstructured
** communications library for the matching code */
static int communication_by_plan (ZZ* zz, int sendcnt, int* dest, int* size, 
 int scale, int* send, int* reccnt, int* recsize, int* nRec, int** rec,
 MPI_Comm comm, int tag);
 
/****************************************************************************/
/* Because this calculation is done in two locations it has been converted to
** a subroutine to assure it is always consistant. Inline is not yet legal! */
static int calc_nCandidates (int num_vtx, int procs)
{
/* 2 below because each match pairs 2 vertices */
return num_vtx ? 1 + num_vtx/(2 * procs * ROUNDS_CONSTANT) : 0;
}
 

/****************************************************************************/
static int pmatching_ipm(
  ZZ *zz,
  HGraph* hg,
  Matching match,
  PHGPartParams *hgp
)
{
  int i, j, k, n, m, round, pvisit, kstart, *r, *s;   /* loop counters */
  int lno, gno, count, old_kstart;                    /* temp variables */
  int sendcnt, sendsize, reccnt, recsize, msgsize;        /* temp variables */
  int nRounds;                   /* # of matching rounds to be performed; 
                                    identical on all procs in 
                                    hgc->Communicator.*/
  int nCandidates;               /* # of candidates on this proc; identical
                                    on all procs in hgc->col_comm. */
  int nTotal;                    /* Total # of candidates; sum of nCandidates
                                    across proc rows. */
  int *send, *dest, *size, *rec, *index, *aux; /* working buffers */
  int nSend, nDest, nSize, nRec, nIndex, nAux; /* currently allocated size of
                                                  the above working buffers */
  float *sums, *f, bestsum;                           /* working buffer, float*/
  int *visit, *cmatch, *select, *permute, *edgebuf;   /* fixed usage arrays */
  int nEdgebuf, nPermute;                             /* array size in ints */
  PHGComm *hgc = hg->comm;
  int err = ZOLTAN_OK, old_row, row, col;
  int gmax_nPins, gmax_nVtx;     /* Global max # pins/proc and vtx/proc */
  int **rows = NULL;             /* used only in merging process */
  int bestlno, vertex, nselect;
  char *yo = "pmatching_ipm";
  int *master_data = NULL, *master_procs = NULL, *mp = NULL, nmaster = 0;
  int cFLAG, edge;                 /* column match only if user requested */
  static int development_timers[2] = {-1, -1};

  
  /* this restriction will be removed later, but for now NOTE this test */
  if (sizeof(int) < sizeof (float))  {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Code must be modified before using");
    err = ZOLTAN_FATAL;
    goto fini;
  }

  /* test if user only wants a column matching */
  cFLAG = strcasecmp(hgp->redm_str, "c-ipm") ? 0 : 1;

  /* determine basic working parameters */
  /* ERIK: need to calculate nCandidates based on # of unmatched vertices */
  nRounds     = hgc->nProc_x * ROUNDS_CONSTANT;
  nCandidates = calc_nCandidates (hg->nVtx, hgc->nProc_x);
  
  /* determine maximum global number of Vtx and Pins for storage allocation */
  /* determine initial sum of all candidates, nTotal for storage allocation */
  MPI_Allreduce(&hg->nPins, &gmax_nPins, 1, MPI_INT, MPI_MAX, hgc->Communicator);
  gmax_nVtx = nTotal = 0;
  for (i = 0; i < hgc->nProc_x; i++)  {
    count = hg->dist_x[i+1]-hg->dist_x[i];    /* number of vertices on proc i */
    if (count > gmax_nVtx)
       gmax_nVtx = count;
    nTotal += calc_nCandidates (count, hgc->nProc_x);
  }
                 
  /* allocate working and fixed sized array storage */
  sums = NULL;
  send = dest = size = rec = index = aux = visit = cmatch = select = NULL;
  permute = edgebuf = NULL;
  nPermute = nIndex = nAux = 1 + MAX(nTotal, gmax_nVtx);
  nDest = nSize = 1 + MAX (hgc->nProc_x, MAX(nTotal, gmax_nVtx));
  nSend = nRec = nEdgebuf = MAX (1000, MAX(gmax_nPins, gmax_nVtx+2));
  
  if (hg->nVtx)  
    if (!(cmatch = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof (int)))
     || !(visit  = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof (int)))
     || !(sums   = (float*) ZOLTAN_CALLOC (hg->nVtx,  sizeof (float))))  {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
       err = ZOLTAN_MEMERR;
       goto fini;
    }
  
  if (!(permute =    (int*)  ZOLTAN_MALLOC (nPermute     * sizeof (int)))
   || !(edgebuf =    (int*)  ZOLTAN_MALLOC (nEdgebuf     * sizeof (int)))
   || !(select  =    (int*)  ZOLTAN_MALLOC ((1+nCandidates)  * sizeof (int)))        
   || !(send    =    (int*)  ZOLTAN_MALLOC (nSend        * sizeof (int)))
   || !(dest    =    (int*)  ZOLTAN_MALLOC (nDest        * sizeof (int)))
   || !(size    =    (int*)  ZOLTAN_MALLOC (nSize        * sizeof (int)))
   || !(rec     =    (int*)  ZOLTAN_MALLOC (nRec         * sizeof (int)))
   || !(index   =    (int*)  ZOLTAN_MALLOC (nIndex       * sizeof (int)))
   || !(aux     =    (int*)  ZOLTAN_MALLOC (nAux         * sizeof (int)))
   || !(rows    =    (int**) ZOLTAN_MALLOC ((hgc->nProc_y + 1) * sizeof (int*)))
   || !(master_data =(int*)  ZOLTAN_MALLOC (3 * nTotal   * sizeof (int)))
   || !(master_procs=(int*)  ZOLTAN_MALLOC (nTotal       * sizeof (int))))  {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
     err = ZOLTAN_MEMERR;
     goto fini;
  }
    
  /* match[] is local slice of global matching array.  It uses local numbering 
   * (zero-based). Initially, match[i] = i. After matching, match[i]=i indicates
   * an unmatched vertex. A matching between vertices i & j is indicated by 
   * match[i] = j & match [j] = i.  NOTE: a match to an off processor vertex is
   * indicated by a negative number, -(gno+1), and must use global numbers
   * (gno's) and not local numbers, lno's which are zero based.        */

  /* Compute candidates' vertex visit order (selection). Random is default. */
  Zoltan_PHG_Vertex_Visit_Order(zz, hg, hgp, visit);
  
  /* Loop processing ncandidates vertices per column each round.
   * Each loop has 4 phases:
   * Phase 1: send ncandidates vertices for global matching - horizontal comm
   * Phase 2: sum  inner products, find best in column - vertical communication
   * Phase 3: return best sums to owner's column    - horizontal communication
   * Phase 4: return actual match selections        - horizontal communication
   *
   * No conflict resolution required because temp locking prevents conflicts. */
  
  pvisit = 0;                                    /* marks position in visit[] */
  for (round = 0; round < nRounds; round++)  {
    if (cFLAG)  {
      for (nTotal = i = 0; i < hg->nVtx; i++)
        if (match[i] == i)
          permute[nTotal++] = i;
      goto skip_phase1;
    } 

    
    /************************ PHASE 1: ***************************************/
    
    
    mp = master_data;
    nmaster = 0;                   /* count of data accumulted in master row */    
    memcpy (cmatch, match, hg->nVtx * sizeof(int));  /* for temporary locking */
                
    /* Select upto nCandidates unmatched vertices to globally match. */
    for (sendcnt = 0; sendcnt < nCandidates && pvisit < hg->nVtx; pvisit++) {
      if (cmatch[visit[pvisit]] == visit[pvisit])  {         /* unmatched */
        select[sendcnt++] = visit[pvisit];      /* select it as a candidate */
        cmatch[visit[pvisit]] = -1;             /* mark it as a pending match */
      }  
    }
    nselect = sendcnt;                          /* save for later use */
                        
    /* fill send buff as list of <gno, gno's edge count, list of edge lno's> */
    /* NOTE: can't overfill send buffer by definition of initial sizing */
    s = send;
    for (i = 0; i < sendcnt; i++)   {
      lno = select[i];
      *s++ = VTX_LNO_TO_GNO (hg, lno);                                 /* gno */
      *s++ = hg->vindex[lno+1] - hg->vindex[lno];                    /* count */
      for (j = hg->vindex[lno]; j < hg->vindex[lno+1]; j++)  
        *s++ = hg->vedge[j];                                   /* lno of edge */
    }        
    sendsize = s - send;           
    
    /* determine actual global number of candidates this round */
    MPI_Allreduce (&sendcnt, &nTotal, 1, MPI_INT, MPI_SUM, hgc->row_comm);     
    if (nTotal == 0)   
      break;                            /* globally all work is done, so quit */

    /* communication to determine global size and displacements of rec buffer */
    MPI_Allgather (&sendsize, 1, MPI_INT, size, 1, MPI_INT, hgc->row_comm); 
     
    /* determine the size of the rec buffer & reallocate bigger iff necessary */
    recsize = 0;
    for (i = 0; i < hgc->nProc_x; i++)
      recsize += size[i];            /* compute total size of edgebuf in ints */
    if (recsize > nEdgebuf)  {
      nEdgebuf = recsize;
      if (!(edgebuf = (int*) ZOLTAN_REALLOC(edgebuf, nEdgebuf * sizeof(int)))) {
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto fini;
      }
    }                          
    
    /* setup displacement array necessary for MPI_Allgatherv */
    dest[0] = 0;
    for (i = 1; i < hgc->nProc_x; i++)
      dest[i] = dest[i-1] + size[i-1];

    /* send vertices & their edges to all row neighbors */
    MPI_Allgatherv (send, sendsize, MPI_INT, edgebuf, size, dest, MPI_INT,
     hgc->row_comm);

    /************************ PHASE 2: ***************************************/
         
    /* create random permutation of index into the edge buffer */
    i = 0;
    for (j = 0 ; j < nTotal  &&  i < recsize; j++)   {
      permute[j] = i++;             /* save index of gno in permute[] */
      count      = edgebuf[i++];    /* count of edges */
      i += count;                   /* skip over count edges */
    }

skip_phase1:
    /* Communication grouped candidates by processor, scramble them! */
    /* Otherwise all candidates from proc column 0 will be matched first. */
    if (hgc->nProc_x > 1 || cFLAG)  {
      /* Future: Instead of Zoltan_Rand_Perm_Int, we could use 
         Zoltan_PHG_Vertex_Visit_Order() to reorder the candidates
         but that routine uses a local hg so won't work on the candidates. */
      Zoltan_Srand_Sync(Zoltan_Rand(NULL), &(hgc->RNGState_col), hgc->col_comm);
      Zoltan_Rand_Perm_Int (permute, nTotal, &(hgc->RNGState_col));
    }

    /* for each candidate vertex, compute all local partial inner products */
    kstart = old_kstart = 0;         /* next candidate (of nTotal) to process */
    while (kstart < nTotal)  {
      sendsize = 0;                      /* position in send buffer */
      sendcnt = 0;                       /* count of messages in send buffer */
      s = send;                          /* start at send buffer origin */
      for (k = kstart; k < nTotal; k++)   {   
        r     = &edgebuf[permute[k]];     
        gno   = *r++;                        /* gno of candidate vertex */
        count = *r++;                        /* count of following hyperedges */
        
        if (cFLAG)
          gno = permute[k];                  /* need to use next local vertex */
 
        if (hgp->use_timers > 3)  {
          if (development_timers[0] < 0)
            development_timers[0] = Zoltan_Timer_Init(zz->ZTime, 0, 
                                                      "inner products");
          ZOLTAN_TIMER_START(zz->ZTime, development_timers[0], 
                             hg->comm->Communicator);
        }
                  
        /* now compute the row's nVtx inner products for kth candidate */
        m = 0;
        if (!cFLAG && (hg->ewgt == NULL) && (hgp->vtx_scal == NULL))
          for (i = 0; i < count; r++, i++)
            for (j = hg->hindex [*r]; j < hg->hindex [*r + 1]; j++) {
              if (sums [hg->hvertex[j]] == 0.0) 
                index [m++] = hg->hvertex[j];
              sums [hg->hvertex[j]] += 1.0;
            }
        else if (!cFLAG && (hg->ewgt == NULL) && (hgp->vtx_scal != NULL))
          for (i = 0; i < count; r++, i++)
            for (j = hg->hindex [*r]; j < hg->hindex [*r + 1]; j++) {
              if (sums [hg->hvertex[j]] == 0.0)
                index [m++] = hg->hvertex[j];
              sums [hg->hvertex[j]] += hgp->vtx_scal[hg->hvertex[j]];
            }
        else if (!cFLAG && (hg->ewgt != NULL) && (hgp->vtx_scal == NULL))
          for (i = 0; i < count; r++, i++)
            for (j = hg->hindex [*r]; j < hg->hindex [*r + 1]; j++)  {
              if (sums [hg->hvertex[j]] == 0.0)
                index [m++] = hg->hvertex[j];
              sums [hg->hvertex[j]] += hg->ewgt [*r];
            }
        else if (!cFLAG && (hg->ewgt != NULL) && (hgp->vtx_scal != NULL))
          for (i = 0; i < count; r++, i++)
            for (j = hg->hindex [*r]; j < hg->hindex [*r + 1]; j++)  {
              if (sums [hg->hvertex[j]] == 0.0)
                index [m++] = hg->hvertex[j];
              sums [hg->hvertex[j]] +=hgp->vtx_scal[hg->hvertex[j]]*hg->ewgt[*r];
            }

        else if (cFLAG && (hg->ewgt == NULL) && (hgp->vtx_scal == NULL))
          for (i = hg->vindex[gno]; i < hg->vindex[gno+1]; i++)  {
            edge = hg->vedge[i];
            for (j = hg->hindex [edge]; j < hg->hindex [edge+1]; j++) {
              if (sums [hg->hvertex[j]] == 0.0)
                index [m++] = hg->hvertex[j];
              sums [hg->hvertex[j]] += 1.0;
            }
          }  
        else if (cFLAG && (hg->ewgt == NULL) && (hgp->vtx_scal != NULL))
          for (i = hg->vindex[gno]; i < hg->vindex[gno+1]; i++)  {
            edge = hg->vedge[i];
            for (j = hg->hindex [edge]; j < hg->hindex [edge+1]; j++) {
              if (sums [hg->hvertex[j]] == 0.0)
                index [m++] = hg->hvertex[j];
              sums [hg->hvertex[j]] += hgp->vtx_scal[hg->hvertex[j]];
            }
          }
        else if (cFLAG && (hg->ewgt != NULL) && (hgp->vtx_scal == NULL))
          for (i = hg->vindex[gno]; i < hg->vindex[gno+1]; i++)  {
            edge = hg->vedge[i];
            for (j = hg->hindex [edge]; j < hg->hindex [edge+1]; j++) {
              if (sums [hg->hvertex[j]] == 0.0)
                index [m++] = hg->hvertex[j];
              sums [hg->hvertex[j]] += hg->ewgt [edge];
            }
          }
        else if (cFLAG && (hg->ewgt != NULL) && (hgp->vtx_scal != NULL))
          for (i = hg->vindex[gno]; i < hg->vindex[gno+1]; i++)  {
            edge = hg->vedge[i];
            for (j = hg->hindex [edge]; j < hg->hindex [edge+1]; j++) {
              if (sums [hg->hvertex[j]] == 0.0)
                index [m++] = hg->hvertex[j];
              sums [hg->hvertex[j]] +=hgp->vtx_scal[hg->hvertex[j]]
               *hg->ewgt[edge];
            }
          }

        if (hgp->use_timers > 3)
          ZOLTAN_TIMER_STOP(zz->ZTime, development_timers[0]);
          
        /* if local vtx, remove self inner product (useless maximum) */
        if (cFLAG)
          sums [gno] = 0.0;             /* here gno is really a local id */
        else if (VTX_TO_PROC_X (hg, gno) == hgc->myProc_x)
          sums [VTX_GNO_TO_LNO (hg, gno)] = 0.0;
         
        /* count unmatched partial sums exceeding PSUM_THRESHOLD */   
        count = 0;
        for (i = 0; i < m; i++)  {
          lno = index[i];
          if (match[lno] == lno  &&  sums[lno] > PSUM_THRESHOLD)
            aux[count++] = lno;      /* save lno for significant partial sum */
          else
            sums[lno] = 0.0;         /* clear unwanted entries */  
        }     
        if (count == 0)
          continue;

        /* HEADER_COUNT (row, col, gno, count of <lno, psum> pairs) */                    
        msgsize = HEADER_COUNT + 2 * count;
        
        /* iff necessary, resize send buffer to fit at least first message */
        if (sendcnt == 0 && (msgsize > nSend))  {
          nSend += msgsize;         
          send = (int*) ZOLTAN_REALLOC(send, nSend * sizeof(int));
          if (send == NULL) {
            ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
            err = ZOLTAN_MEMERR;
            goto fini;
          }
          s = send;    
        }

        if (sendsize + msgsize <= nSend)  {
          /* current partial sums fit, so put them into the send buffer */
          dest[sendcnt]   = gno % hgc->nProc_y;  /* proc to compute tsum */
          size[sendcnt++] = msgsize;             /* size of message */
          sendsize     += msgsize;             /* cummulative size of message */
          
          *s++ = hgc->myProc_y;      /* save my row (for merging) */
          *s++ = hgc->myProc_x;      /* save my col (for debugging) */
          *s++ = gno;          
          *s++ = count;
          for (i = 0; i < count; i++)  {          
            *s++ = aux[i];                          /* lno of partial sum */
             f = (float*) s++;
            *f = sums[aux[i]];                      /* partial sum */           
            sums[aux[i]] = 0.0;
          }          
        }
        else  {           /* psum message doesn't fit into buffer */
          for (i = 0; i < count; i++)              
            sums[aux[i]] = 0.0;        
          break;
        }  
      }                  /* DONE: loop over k */                    

      if (hgp->use_timers > 3)  {
        if (development_timers[1] < 0)
          development_timers[1] = Zoltan_Timer_Init(zz->ZTime, 0, 
                                                   "build totals");
        ZOLTAN_TIMER_START(zz->ZTime, development_timers[1], 
                           hg->comm->Communicator);
      }      
      
      /* synchronize all rows in this column to next kstart value */
      old_kstart = kstart;      
      MPI_Allreduce (&k, &kstart, 1, MPI_INT, MPI_MIN, hgc->col_comm);
            
      /* Send inner product data in send buffer to appropriate rows */
      err = communication_by_plan (zz, sendcnt, dest, size, 1, send, &reccnt, 
       &recsize, &nRec, &rec, hgc->col_comm, IPM_TAG);
      if (err != ZOLTAN_OK)
        goto fini;
      
      /* build index into receive buffer pointer for each new row of data */
      old_row = -1;
      k = 0;
      for (r = rec; r < rec + recsize  &&  k < hgc->nProc_y; )  {     
        row = *r++;        
        col = *r++;               /* column only for debugging, may go away */
        if (row != old_row)  {
          index[k++] = r - rec;   /* points at gno, not row or col */
          old_row = row;
        }
        gno   = *r++;
        count = *r++;
        r += (count * 2);
      }
     
      /* save current positions into source rows within rec buffer */
      for (i = 0; i < k; i++)
        rows[i] = &rec[index[i]];
      for (i = k; i < hgc->nProc_y; i++)
        rows[i] = &rec[recsize];       /* in case no data came from a row */
      rows[i] = &rec[recsize];         /* sentinel */
      
      /* merge partial i.p. sum data to compute total inner products */ 
      s = send; 
      for (n = old_kstart; n < kstart; n++)  {
        m = 0;        
        gno = (cFLAG) ? permute[n] : edgebuf [permute[n]];
        
#ifdef RTHRTH        
        /* Not sure if this test makes any speedup ???, works without! */
        if (gno % hgc->nProc_y != hgc->myProc_y)
          continue;                           /* this gno is not on this proc */
#endif

        /* merge step: look for target gno from each row's data */
        for (i = 0; i < hgc->nProc_y; i++)  {
          if (rows[i] < &rec[recsize] && *rows[i] == gno)  {       
            count = *(++rows[i]);
            for (j = 0; j < count; j++)  {
              lno = *(++rows[i]);         
              if (sums[lno] == 0.0)       /* is this first time for this lno? */
                aux[m++] = lno;           /* then save the lno */          
              sums[lno] += *(float*) (++rows[i]);    /* sum the psums */
            }
            rows[i] += 3;                 /* skip past current lno, row, col */              
          }
        }
          
        /* determine how many total inner products exceed threshold */  
        count = 0;
        for (i = 0; i < m; i++)
          if (sums[aux[i]] > TSUM_THRESHOLD)
            count++;   

        /* create <gno, count, <lno, tsum>> in send array */           
        if (count > 0)  {
          if ( (s - send) + (2 + 2 * count) > nSend ) {
            sendsize = s - send;
            nSend += (2 + 2 * count);
            send = (int*) ZOLTAN_REALLOC (send, nSend * sizeof(int));
            if (send == NULL)  {
              ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory");
              return ZOLTAN_MEMERR;
            }
           s = send + sendsize;   /* since realloc buffer could move */ 
          }      
          *s++ = gno;
          *s++ = count;
        }  
        for (i = 0; i < m; i++)   {
          lno = aux[i];             
          if (sums[lno] > TSUM_THRESHOLD)  {
            *s++ = lno;
             f = (float*) s++;
            *f = sums[lno];
          }  
          sums[lno] = 0.0;  
        }     
      }
      sendsize = s - send;   /* size (in ints) of send buffer */
      
      /* Communicate total inner product results to MASTER ROW */
      MPI_Allgather (&sendsize, 1, MPI_INT, size, 1, MPI_INT, hgc->col_comm);
      recsize = 0;
      for (i = 0; i < hgc->nProc_y; i++)
        recsize += size[i];        
          
      dest[0] = 0;
      for (i = 1; i < hgc->nProc_y; i++)
        dest[i] = dest[i-1] + size[i-1];
        
      if (recsize > nRec) {
        nRec = recsize;
        ZOLTAN_FREE (&rec);
        rec = (int*) ZOLTAN_MALLOC (nRec * sizeof(int));
        if (rec == NULL)  {
          ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory");
          return ZOLTAN_MEMERR;
        }
      }    
      if (recsize)  
        MPI_Gatherv (send, sendsize, MPI_INT, rec, size, dest, MPI_INT, 0,
         hgc->col_comm);
       
      /* Determine best vertex and best sum for each candidate */
      if (hgc->myProc_y == 0)
        for (r = rec; r < rec + recsize;)  {
          gno   = *r++;                    /* candidate's GNO */
          count = *r++;                    /* count of nonzero pairs */
          bestsum = -1.0;                  /* any negative value will do */
          bestlno = -1;                    /* any negative value will do */
          for (i = 0; i < count; i++)  {
            lno =          *r++;
            f   =  (float*) r++;       
            if (cFLAG  && *f > bestsum  &&  match[lno] == lno)  {
              bestsum = *f;
              bestlno = lno;
            }
                                         
            if (!cFLAG && *f > bestsum  &&  cmatch[lno] == lno)  {
              bestsum = *f;
              bestlno = lno;
            }
          }

          if (cFLAG  && match[gno] == gno && (bestsum > TSUM_THRESHOLD))
            {
            match[bestlno] = gno;
            match[gno]     = bestlno;
            }
                        
          if (!cFLAG && bestsum > TSUM_THRESHOLD)  {
            cmatch[bestlno] = -1;  /* mark pending match to avoid conflicts */
            *mp++ = gno;
            *mp++ = VTX_LNO_TO_GNO (hg, bestlno);
             f = (float*) mp++;
            *f = bestsum;
            master_procs[nmaster++] = VTX_TO_PROC_X (hg, gno);
          }
        } 
      if (cFLAG && hgc->nProc_y > 1)  {
        /* Broadcast what we matched so far */
        MPI_Bcast (match, hg->nVtx, MPI_INT, 0, hgc->col_comm); 
      }
      
      if (hgp->use_timers > 3)
        ZOLTAN_TIMER_STOP(zz->ZTime, development_timers[1]);
          
    }                                       /* DONE: kstart < nTotal loop */
    if (cFLAG)  {
      break;      /* done, no more phases (3 or 4) or rounds */
    }    

    
         /************************ PHASE 3: **********************************/
    
    /* MASTER ROW only: send best results to candidates' owners */
    err = communication_by_plan (zz, nmaster, master_procs, NULL, 3,
     master_data, &reccnt, &recsize, &nRec, &rec, hgc->row_comm, IPM_TAG+5);
    if (err != ZOLTAN_OK)
      goto fini;  

    /* read each message (candidate id, best match id, and best i.p.) */ 
    if (hgc->myProc_y == 0) 
      for (r = rec; r < rec + 3 * reccnt; )   {
        gno    = *r++;
        vertex = *r++;
        f      = (float*) r++;
        bestsum = *f;            
                             
        /* Note: ties are broken to favor local over global matches */   
        lno =  VTX_GNO_TO_LNO (hg, gno);  
        if ((bestsum > sums [lno]) || (bestsum == sums[lno]
         && VTX_TO_PROC_X (hg, gno)    != hgc->myProc_x
         && VTX_TO_PROC_X (hg, vertex) == hgc->myProc_x))    {        
            index [lno] = vertex;
            sums  [lno] = bestsum;
        }                   
      }   
                                    
    /************************ PHASE 4: ************************************/
        
    /* fill send buffer with messages. A message is two matched gno's */
    /* Note: match to self if inner product is below threshold */
    s = send; 
    sendcnt = 0;
    if (hgc->myProc_y == 0)
      for (i = 0; i < nselect; i++)   {
        int d1, d2;
        lno = select[i];
        *s++ = gno = VTX_LNO_TO_GNO (hg, lno);
        *s++ = vertex = (sums [lno] > TSUM_THRESHOLD) ? index[lno] : gno;
            
        /* each distict owner (gno or vertex) needs its copy of the message */
        /*  KDDKDD:  Matching gno to vertex with IP=sums[lno] */
        d1 = VTX_TO_PROC_X (hg, gno);
        d2 = VTX_TO_PROC_X (hg, vertex);
        dest[sendcnt++] = d1;
        if (d1 != d2)  {
          *s++ = gno;
          *s++ = vertex;        
          dest[sendcnt++] = d2;
        }
      }
        
    /* send match results only to impacted parties */
    err = communication_by_plan (zz, sendcnt, dest, NULL, 2, send, &reccnt,
     &recsize, &nRec, &rec, hgc->row_comm, IPM_TAG+10);
    if (err != ZOLTAN_OK)
      goto fini;

    /* update match array with current selections */
    /* Note: -gno-1 designates an external match as a negative number */
    if (hgc->myProc_y == 0)
      for (r = rec; r < rec + 2 * reccnt; )  {   
        gno    = *r++;
        vertex = *r++;

        if (VTX_TO_PROC_X (hg, gno)    == hgc->myProc_x
         && VTX_TO_PROC_X (hg, vertex) == hgc->myProc_x)   {
            int v1 = VTX_GNO_TO_LNO (hg, vertex);             
            int v2 = VTX_GNO_TO_LNO (hg, gno);                
            match [v1] = v2;
            match [v2] = v1;
        }                         
        else if (VTX_TO_PROC_X (hg, gno) == hgc->myProc_x)
          match [VTX_GNO_TO_LNO (hg, gno)]    = -vertex - 1;
        else              
          match [VTX_GNO_TO_LNO (hg, vertex)] = -gno - 1;
      }      
    
    /* update match array to the entire column */   
    MPI_Bcast (match, hg->nVtx, MPI_INT, 0, hgc->col_comm);
  }                                             /* DONE: loop over rounds */
   
if (0 && hgc->myProc_x == 0 && hgc->myProc_y == 0)
{
int local = 0, global = 0, unmatched = 0;
for (i = 0; i < hg->nVtx; i++)
  {
  if      (match[i] == i)  unmatched++;
  else if (match[i] < 0)   global++;
  else                     local++;
  }
uprintf (hgc, "RTHRTH %d unmatched, %d external, %d local of %d\n",
 unmatched, global, local, hg->nVtx);
}



if (0)
{
/* The following tests that the global match array is a valid permutation */
/* NOTE:  THESE TESTS ARE NOT MANDATORY; THEY CAN BE EXCLUDED AFTER WE 
 * COMPLETE TESTING OF matching_ipm. 
 */

for (i = 0; i < hg->nVtx; i++)
  if (match[i] < 0)  cmatch[i] = -match[i] - 1;
  else               cmatch[i] = VTX_LNO_TO_GNO (hg, match[i]);

MPI_Allgather (&hg->nVtx, 1, MPI_INT, size, 1, MPI_INT, hgc->row_comm); 

recsize = 0;
for (i = 0; i < hgc->nProc_x; i++)
  recsize += size[i];
  
dest[0] = 0;
for (i = 1; i < hgc->nProc_x; i++)
  dest[i] = dest[i-1] + size[i-1];
  
if (nRec < recsize)
  rec = (int*) ZOLTAN_REALLOC (rec, recsize * sizeof(int));
MPI_Allgatherv (cmatch, hg->nVtx, MPI_INT, rec, size, dest, MPI_INT, hgc->row_comm);

if (nSend < recsize)
  send = (int*) ZOLTAN_REALLOC (send, recsize * sizeof(int));
  
for (i = 0; i < recsize; i++)
  send[i] = 0;
for (i = 0; i < recsize; i++)
  ++send[rec[i]];

count = 0;
for (i = 0; i < recsize; i++)
  if (send[i] != 1)
    count++;
if (count)    
  uprintf (hgc, "RTHRTH %d FINAL MATCH ERRORS of %d\n", count, recsize); 
}

  if (hgp->use_timers > 3) {
    Zoltan_Timer_Print(zz->ZTime, development_timers[0], zz->Proc, 
                       hg->comm->Communicator, stdout);
    Zoltan_Timer_Print(zz->ZTime, development_timers[1], zz->Proc, 
                       hg->comm->Communicator, stdout);
  }
     
fini:
  Zoltan_Multifree (__FILE__, __LINE__, 15, &cmatch, &visit, &sums, &send,
   &dest, &size, &rec, &index, &aux, &permute, &edgebuf, &select, &rows,
   &master_data, &master_procs);
  return err;
}


  
/****************************************************************************/
static int communication_by_plan (ZZ* zz, int sendcnt, int* dest, int* size, 
 int scale, int* send, int* reccnt, int *recsize, int* nRec, int** rec,
 MPI_Comm comm, int tag)
{
   ZOLTAN_COMM_OBJ *plan = NULL;
   int err;
   char *yo = "communication_by_plan";
   
   /* communicate send buffer messages to other row/columns in my comm */  
   err = Zoltan_Comm_Create (&plan, sendcnt, dest, comm, tag, reccnt);
   if (err != ZOLTAN_OK) {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "failed to create plan");
     return err;
   }
        
   /* resize plan if necessary */
   if (size != NULL) {
     err = Zoltan_Comm_Resize (plan, size, tag+1, recsize);
     if (err != ZOLTAN_OK) {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "failed to resize plan");
       return err;
     }
     scale = 1;       
   }
   
   /* realloc rec buffer if necessary */  
   if (*recsize > *nRec)  {   
     *nRec = *recsize;
     if (!(*rec = (int*) ZOLTAN_REALLOC (*rec, *nRec * sizeof(int))))  {
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory");
       return ZOLTAN_MEMERR;
     }
   }
   
   /* send messages from send buffer to destinations */      
   err = Zoltan_Comm_Do (plan, tag+2, (char*) send, scale * sizeof(int),
    (char*) *rec);
   if (err != ZOLTAN_OK)  {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "failed in Comm_Do");
     return err;
   }
   
   /* free memory associated with the plan */
   Zoltan_Comm_Destroy (&plan); 
   return ZOLTAN_OK;
}
   
   

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif


