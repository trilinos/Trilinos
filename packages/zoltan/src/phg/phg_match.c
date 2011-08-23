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
#include "g2l_hash.h"
#include "zz_util_const.h"

static ZOLTAN_PHG_MATCHING_FN pmatching_ipm;         /* inner product matching */
static ZOLTAN_PHG_MATCHING_FN pmatching_agg_ipm;     /* agglomerative IPM */
static ZOLTAN_PHG_MATCHING_FN pmatching_geom;        /* geometric matching */
static ZOLTAN_PHG_MATCHING_FN pmatching_local_ipm;   /* local ipm */
static ZOLTAN_PHG_MATCHING_FN pmatching_alt_ipm;     /* alternating ipm */
static ZOLTAN_PHG_MATCHING_FN pmatching_hybrid_ipm;  /* hybrid ipm */

static int Zoltan_PHG_match_isolated(ZZ* zz, HGraph* hg, PHGPartParams* hgp,
                                     ZOLTAN_GNO_TYPE *match, int small_degree);
static int Zoltan_PHG_compute_esize(ZZ* zz, HGraph* hg, PHGPartParams* hgp);

/* Functions to be used in geometric matching */
static ZOLTAN_NUM_OBJ_FN geometric_get_num_obj;
static ZOLTAN_OBJ_LIST_FN geometric_get_obj_list;
static ZOLTAN_NUM_GEOM_FN geometric_get_num_geom;
static ZOLTAN_GEOM_MULTI_FN geometric_get_geom_multi;
  
typedef struct triplet {
    ZOLTAN_GNO_TYPE candidate;      /* gno of candidate vertex */
    ZOLTAN_GNO_TYPE partner;        /* gno of best match found so far */
    float ip;           /* total inner product between candidate and partner */
} Triplet;

    
static HGraph *HG_Ptr;

/*
  Matching is NOT restricted (called as MATCH_OK) if
  1- we are not doing fixed vertex partitioning, or partitioning with preferred parts
     (i.e. !hgp->UsePrefPart)
  2- if one of the vertices is free (i.e. fv? < 0);
     in other words, it is OK if only one of them is fixed
  3- if both vertices are fixed in the "same" part; note that since we are achieving k-way
     partitioning via recursive bisection; "same" is defined as both vertices are on the
     same side of that level's bisection. For example for 4-way partitioning since
     bisec_split will be 2 in the first bisection; we'll allow vertices fixed in parts
     {0, 1} are matched among themselves, also the vertices fixed in {2, 3} are allowed to
     be matched among themselves. 
 */   
#define MATCH_OK(hgp, hg, fv1, fv2) \
        (!((hgp)->UsePrefPart)   ||     \
         ((fv1) < 0) ||     \
         ((fv2) < 0) ||     \
         ((((fv1) < (hg)->bisec_split) ? 0 : 1) == (((fv2) < (hg)->bisec_split) ? 0 : 1)))
    
/*
 In addition to above conditions, in Agglomerative matching, in order to reduce the extra
 communication, we will only allow matching if candidate vertex is free or its fixed part
 matches with the other vertex's fixed part. Hence the signature of the AGG_MATCH_OK is
 changed to identify which of two argument is candidate vertex's fixed part.
*/
#define AGG_MATCH_OK(hgp, hg, candf, fv2) \
    (!((hgp)->UsePrefPart) || (((candf) < 0) ||  \
     (((fv2)>=0) && ((((candf) < (hg)->bisec_split) ? 0 : 1) == (((fv2) < (hg)->bisec_split) ? 0 : 1)))))

    

/*****************************************************************************/
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams *hgp)
{
    int exist=1;

    if (!strcasecmp(hgp->redm_str, "no"))
        hgp->matching = NULL;
    else if (!strcasecmp(hgp->redm_str, "none"))
        hgp->matching = NULL;
    else if (!strcasecmp(hgp->redm_str, "ipm"))
        hgp->matching = pmatching_ipm;
    else if (!strncasecmp(hgp->redm_str, "agg", 3)) { /* == "agg-ipm" */
        hgp->matching = pmatching_agg_ipm;
        hgp->match_array_type = 1;
    }
    else if (!strncasecmp(hgp->redm_str, "rcb", 3) || 
             !strncasecmp(hgp->redm_str, "rib", 3)) {
        hgp->matching = pmatching_geom;
	hgp->match_array_type = 1;
    }
    else if (!strcasecmp(hgp->redm_str, "l-ipm"))
        hgp->matching = pmatching_local_ipm;           
    else if (!strcasecmp(hgp->redm_str, "c-ipm"))
        hgp->matching = pmatching_ipm;   
    else if (!strcasecmp(hgp->redm_str, "a-ipm"))
        hgp->matching = pmatching_alt_ipm;
    else if (!strcasecmp(hgp->redm_str, "h-ipm"))
        hgp->matching = pmatching_hybrid_ipm;
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
  ZOLTAN_GNO_TYPE *match,
  PHGPartParams *hgp)
{
float *old_ewgt = NULL, *new_ewgt = NULL;
int   ierr = ZOLTAN_OK;
char  *yo = "Zoltan_PHG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hgp->edge_scaling) {
     if (hg->nEdge && !(new_ewgt = (float*) 
                      ZOLTAN_MALLOC(hg->nEdge * sizeof(float))))
         MEMORY_ERROR;
 
     Zoltan_PHG_Scale_Edges (zz, hg, new_ewgt, hgp->edge_scaling);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
  }

  /* Create/update scale vector for vertices for inner product */
  if (hgp->vtx_scaling) 
     Zoltan_PHG_Scale_Vtx (zz, hg, hgp);
  
  /* Do the matching */
  if (hgp->matching) {
    /* first match isolated vertices */
    /* if matching is not geometric, then run isolated */
    if (hgp->matching != pmatching_geom)
      Zoltan_PHG_match_isolated(zz, hg, hgp, match, 0);

    /* now do the real matching */
    if ((ierr = Zoltan_PHG_compute_esize(zz, hg, hgp))==ZOLTAN_OK)
        ierr = hgp->matching (zz, hg, match, hgp);
    /* clean up by matching "near-isolated" vertices of degree 1 */
    /* only useful in special cases (e.g. near-diagonal matrices).
       in many cases there is a slight increase in cuts, so 
       turn it off for now.                                      */
    /* Zoltan_PHG_match_isolated(zz, hg, hgp, match, 1); */
  }

End: 

  /* Restore the old edge weights if scaling was used. */
  if (hgp->edge_scaling)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE (&new_ewgt);
  ZOLTAN_TRACE_EXIT (zz, yo);
  return ierr;
}



/***************************************************************************/
/* MPI_Op for computing the maximum inner-product for each candidate */
static void phasethreereduce (
  void *tin, 
  void *tinout, 
  int  *tnum, 
  MPI_Datatype *mytype)
{
  int num = *tnum;
  Triplet *in    = (Triplet*) tin;
  Triplet *inout = (Triplet*) tinout;
  int i;

  for (i = 0; i < num; i++) {
    if (in[i].candidate == -1 && inout[i].candidate == -1) {
      continue;                         /* No values set for this candidate */
    }

/* Assumption - candidates are the same in in and inout, unless one is -1 */
/*
    if (inout[i].candidate == -1){
      inout[i] = in[i];
    }
    else if (in[i].candidate != -1){
*/
      if (in[i].ip > inout[i].ip)
        inout[i] = in[i];                 /* in has larger inner product */
      else if (in[i].ip == inout[i].ip) {
        int in_proc    = VTX_TO_PROC_X (HG_Ptr, in[i].partner);
        int inout_proc = VTX_TO_PROC_X (HG_Ptr, inout[i].partner);
        int cand_proc  = VTX_TO_PROC_X (HG_Ptr, in[i].candidate);
  
  
        /* Give preference to partners on candidate's processor */
        if (((in_proc == cand_proc) && (inout_proc == cand_proc))
         || ((in_proc != cand_proc) && (inout_proc != cand_proc))) {
          /* Both partners are on candidate's proc OR
             neither partner is on candidate's proc.
             Break ties by larger partner gno. */
          if (in[i].partner > inout[i].partner)
            inout[i] = in[i];
        }
        else if (in_proc == cand_proc) {
          inout[i] = in[i];   /* Give preference to local partner */
        }
      } 
/*
    }
*/
  }
}


/* UVCUVC TODO CHECK: later consider about reusing edge size computed in
   recursive bisection split; but that also requires communication of
   those values in redistribution */ 
static int Zoltan_PHG_compute_esize(
  ZZ *zz,
  HGraph *hg,
  PHGPartParams *hgp
)
{
    int  i;
    MPI_Datatype zoltan_gno_mpi_type;
    ZOLTAN_GNO_TYPE *lsize = NULL;
    static char *yo = "Zoltan_PHG_compute_esize";
    int ierr = ZOLTAN_OK;

    if (hg->nEdge && !hg->esize) {
      if (!(lsize = (ZOLTAN_GNO_TYPE*)  ZOLTAN_MALLOC(hg->nEdge*sizeof(ZOLTAN_GNO_TYPE))))
        MEMORY_ERROR;
      if (!(hg->esize = (ZOLTAN_GNO_TYPE*)  ZOLTAN_MALLOC(hg->nEdge*sizeof(ZOLTAN_GNO_TYPE))))
        MEMORY_ERROR;

      for (i=0; i<hg->nEdge; ++i)
          lsize[i] = (ZOLTAN_GNO_TYPE)(hg->hindex[i+1] - hg->hindex[i]);

      zoltan_gno_mpi_type = Zoltan_mpi_gno_type();
      MPI_Allreduce(lsize, hg->esize, hg->nEdge, zoltan_gno_mpi_type, MPI_SUM, hg->comm->row_comm);
      
End:
      ZOLTAN_FREE(&lsize);
    }
    return ierr;
}


static int Zoltan_PHG_match_isolated(
  ZZ *zz,
  HGraph *hg,
  PHGPartParams *hgp,
  ZOLTAN_GNO_TYPE *match,
  int small_degree /* 0 or 1; 0 corresponds to truely isolated vertices */
)
{
    int v=-1, i, unmatched=0, *ldeg, *deg;
#ifdef _DEBUG
    int cnt=0;
#endif
    static char *yo = "Zoltan_PHG_match_isolated";
    int ierr = ZOLTAN_OK;

    if (hg->nVtx) {
      if (!(ldeg = (int*)  ZOLTAN_MALLOC(2*hg->nVtx*sizeof(int))))
        MEMORY_ERROR;
      deg = ldeg + hg->nVtx;
      /* match isolated vertices.
         UVCUVC: right now we match in the natural order,
         I don't think we need random matching but if needed
         we can match in random order. */
      for (i=0; i<hg->nVtx; ++i)
          ldeg[i] = hg->vindex[i+1] - hg->vindex[i];
      MPI_Allreduce(ldeg, deg, hg->nVtx, MPI_INT, MPI_SUM, hg->comm->col_comm);
      
      if (small_degree>0){
          /* Only match on procs with many unmatched vertices */
          unmatched= 0;
          for (i=0; i<hg->nVtx; ++i)
              if (match[i]==i) unmatched++;
      }
      if ((small_degree==0) || (unmatched > 0.8*hg->nVtx))
          for (i=0; i<hg->nVtx; ++i){
              if ((match[i]==i) && (deg[i] <= small_degree)) { 
#ifdef _DEBUG
                  ++cnt;
#endif
                  /* match with previous unmatched vertex */
                  /* EBEB For degree-1 vertices, we could be more clever
                     and match vertices that share a common neighbor */
                  if (v==-1)
                      v = i;
                  else if (MATCH_OK(hgp, hg, hg->pref_part[i], hg->pref_part[v])) {
                      match[v] = i;
                      match[i] = v;
                      v = -1;
                  }
              }
          }
#ifdef _DEBUG
      if (cnt)
          uprintf(hg->comm, "Local H(%d, %d, %d) and there were %d isolated vertices\n", hg->nVtx, hg->nEdge, hg->nPins, cnt);           
#endif
End:
      ZOLTAN_FREE(&ldeg);
    }
    return ierr;
}


/*****************************************************************************/
/* inner product matching                                                    */ 
/* based on implementations by Rob Bisseling and Umit Catalyurek             */
/* for each vertex, we match with the unmatched vertex which has the most    */
/* hyperedges in common with it (ie, the pair with greatest inner product).  */
/* by Aaron Becker, UIUC, Summer 2004                                        */
/* 8/5/04  Erik says matching_ipm is nearly equivalent to matching_rhm;
   but rhm uses a scaled inner product. */
static int matching_ipm(ZZ *zz, HGraph *hg, PHGPartParams *hgp,
                        ZOLTAN_GNO_TYPE *match, int *limit)
{
    int   i, j, k, n;
    float maxip;
    int   *order = NULL;
    float *ips = NULL; 
    ZOLTAN_GNO_TYPE v1, v2, edge, maxindex;
    ZOLTAN_GNO_TYPE *adj = NULL;
    char  *yo = "matching_ipm";

    if (hg->nVtx && 
        (!(ips = (float*) ZOLTAN_MALLOC(hg->nVtx * sizeof(float))) 
     || !(adj = (ZOLTAN_GNO_TYPE*) ZOLTAN_MALLOC(hg->nVtx * sizeof(ZOLTAN_GNO_TYPE)))
     || !(order = (int*) ZOLTAN_MALLOC(hg->nVtx * sizeof(int))))) {
        Zoltan_Multifree(__FILE__, __LINE__, 3, &ips, &adj, &order);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }
    
    for (i = 0; i < hg->nVtx; i++){
        ips[i] = .0;
        order[i] = i;
    }
 
    /* random node visit order is default */
    /* other options may be added later */
    Zoltan_Rand_Perm_Int (order, hg->nVtx, NULL);

    /* for every vertex */
    for (k = 0; k < hg->nVtx  &&  *limit > 0; k++) {
        v1 = order[k];

        if (match[v1] != v1)
            continue;
        
        n = 0; /* number of neighbors */
        /* for every hyperedge containing the vertex */
        for (i = hg->vindex[v1]; i < hg->vindex[v1+1]; i++) {
            edge = hg->vedge[i];

            if (hg->esize[edge] < hgp->MatchEdgeSizeThreshold) {
                /* for every other vertex in the hyperedge */
                for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
                    v2 = hg->hvertex[j];
                    if (match[v2] != v2) {
                        /* row swapping goes here */
                    }
                    else {
                        if (ips[v2]==0.0)
                            adj[n++] = v2;
                        ips[v2] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
                    }
                }
            }
        }

        /* now choose the vector with greatest inner product */
        maxip = 0;
        maxindex = -1;
        for (i = 0; i < n; i++) {
            v2 = adj[i];
            if (ips[v2] > maxip && v2 != v1 && match[v2] == v2
                && MATCH_OK(hgp, hg, hg->pref_part[v1], hg->pref_part[v2])) {
                maxip = ips[v2];
                maxindex = v2;
            }
            ips[v2] = 0;
        }
        if (maxindex != -1) {
            match[v1] = maxindex;
            match[maxindex] = v1;
            (*limit)--;
        } 
        
        /*
          printf("Done with %d, best match is %d with product %d\n",
                  v1, maxindex, maxip);
         */
    }

    /*
      printf("Final Matching:\n");
      for(i = 0; i < hg->nVtx; i++)
          printf("%2d ",i);
      printf("\n");
      for(i = 0; i < hg->nVtx; i++)
          printf("%2zd ",match[i]);
      printf("\n");
    */

    Zoltan_Multifree(__FILE__, __LINE__, 3, &ips, &adj, &order);
    return ZOLTAN_OK;
}


static int pmatching_local_ipm(
  ZZ *zz,
  HGraph *hg,
  ZOLTAN_GNO_TYPE *match,
  PHGPartParams *hgp
)
{
    int limit=hg->nVtx, err=ZOLTAN_OK;
    PHGComm *hgc=hg->comm;
    int root_matchcnt, root_rank;
    MPI_Datatype zoltan_gno_mpi_type;

    /* UVC: In order to simplify adding/testing old sequential matching codes easy
       I'm keeping the interface same; just move your favorite matching code to this file
       and change the function name in the next line */ 
    err = matching_ipm(zz, hg, hgp, match, &limit);
    
    /* UVC: again to simplify testing Optimization; I'm leaving the call here, just
     move your function and rename the call
     if (hgp->matching_opt) 
        err = hgp->matching_opt (zz, hg, match, &limit);
    */
    
    /* find the index of the proc in column group with the best match
       (max #matches); it will be our root proc */
    Zoltan_PHG_Find_Root(hg->nVtx-limit, hgc->myProc_y, hgc->col_comm,
                         &root_matchcnt, &root_rank);
    
    zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

    MPI_Bcast(match, hg->nVtx, zoltan_gno_mpi_type, root_rank, hgc->col_comm);
    
    return err;
}


/**************************************************************************
  Alternating ipm method. Alternate between full ipm and fast method. 
 *************************************************************************/
static int pmatching_alt_ipm(
  ZZ *zz,
  HGraph* hg,
  ZOLTAN_GNO_TYPE *match,
  PHGPartParams *hgp
)
{
  int ierr = ZOLTAN_OK;
  char redm_orig[MAX_PARAM_STRING_LEN];

  strcpy(redm_orig, hgp->redm_str); /* save original parameter string */

  /* first level is 0 */
  if ((hg->info&1) == 0)  /* alternate even-odd levels */
    strcpy(hgp->redm_str, hgp->redm_fast); /* fast method is typically l-ipm */
  else
    strcpy(hgp->redm_str, "ipm");  

  Zoltan_PHG_Set_Matching_Fn (hgp);  /* set pointer to the desired matching function */
  ierr = hgp->matching (zz, hg, match, hgp);
  hgp->matching = pmatching_alt_ipm;  /* reset function pointer */

  /* set redm parameter back to original */
  strcpy(hgp->redm_str, redm_orig);
  
  return ierr;
}

/**************************************************************************
  Hybrid ipm method. First partial c-ipm, then full ipm on unmatched vtxs.
 *************************************************************************/
static int pmatching_hybrid_ipm(
  ZZ *zz,
  HGraph* hg,
  ZOLTAN_GNO_TYPE *match,
  PHGPartParams *hgp
)
{
  int ierr = ZOLTAN_OK;
  char redm_orig[MAX_PARAM_STRING_LEN];

  strcpy(redm_orig, hgp->redm_str); /* save original parameter ("h-ipm") */

  /* Set parameter for how good ip values to keep. */
  /* This could perhaps be a user parameter. */
  hgp->hybrid_keep_factor = 1.0;

  /* First do (partial) c-ipm. */
  strcpy(hgp->redm_str, "c-ipm");
  ierr = pmatching_ipm(zz, hg, match, hgp);  

  /* Reset hybrid_keep_factor to be safe. */
  hgp->hybrid_keep_factor = 0.0;

  /* Then do full ipm on remaining unmatched vertices. */
  strcpy(hgp->redm_str, "ipm");  
  ierr = pmatching_ipm(zz, hg, match, hgp);  

  /* Reset redm parameter back to original */
  strcpy(hgp->redm_str, redm_orig);
  
  return ierr;
}


/****************************************************************************
 * inner product matching (with user selectable column variant, c-ipm)
 * Parallelized version of the serial algorithm (see hg_match.c)
 * Based on conversations with Rob Bisseling by Aaron Becker, UIUC, summer 2004
 * completed by R. Heaphy. New features by Karen Devine & Erik Boman.
 */
               
 /* match[] is local slice of global matching array. Initially, match[i] = i.
  * After matching, match[i]=i indicates an unmatched vertex. A matching between
  * vertices i & j is indicated by match[i] = j & match[j] = i.  NOTE: Positive
  * integers are local numbering (zero-based).  A match to an off processor
  * vertex is indicated by a negative number, -(gno+1), and must use global
  * numbers (gno's) and not local numbers, lno's.                             */ 
 
#define ROUNDS_CONSTANT 8     /* controls the number of candidate vertices */ 
#define IPM_TAG        28731  /* MPI message tag, arbitrary value */
#define CONFLICT_TAG   IPM_TAG+1 

/* these thresholds need to become parameters in phg - maybe ??? */
#define PSUM_THRESHOLD 0.0    /* ignore partial inner products < threshold */
#define TSUM_THRESHOLD 0.0    /* ignore total inner products < threshold */
                                 
/* Forward declaration for a routine that encapsulates the common calls to use
** the Zoltan unstructured communications library for the matching code */
static int communication_by_plan (ZZ* zz, int sendcnt, int* dest, int* size, 
 int scale, char * send, int* reccnt, int* recsize, int* nRec, char ** rec,
 MPI_Comm comm, int tag);


/* Actual inner product calculations between candidates (in rec buffer)    */
/* and local vertices.  Not inlined because ARG is changing for each edge */
/* in inner loop and also uses hgp->MatchEdgeSizeThreshold */
#define INNER_PRODUCT1(ARG)\
  for (i = 0; i < count; intptr++, i++) {\
    if (((ARG)>0.0) && (hg->esize[*intptr] < hgp->MatchEdgeSizeThreshold)) {\
      for (j = hg->hindex[*intptr]; j < hg->hindex[*intptr + 1]; j++) {\
        if (cmatch[hg->hvertex[j]] == hg->hvertex[j])  {\
          if (sums[hg->hvertex[j]] == 0.0)\
            index[m++] = hg->hvertex[j];\
          sums[hg->hvertex[j]] += (ARG);\
        }\
      }\
    }\
  }

#define AGG_INNER_PRODUCT1(ARG)\
  for (i = 0; i < count; intptr++, i++) {\
    if (((ARG)>0.0) && (hg->esize[*intptr]<hgp->MatchEdgeSizeThreshold)) {\
      for (j = hg->hindex[*intptr]; j < hg->hindex[*intptr + 1]; j++) {\
        int v=lhead[hg->hvertex[j]];\
        if (sums[v] == 0.0)\
          aux[m++] = v;\
        sums[v] += (ARG);\
      }\
    }\
  }
  

/* Mostly identical inner product calculation to above for c-ipm variant. Here */
/* candidates are a subset of local vertices and are not in a separate buffer  */
#define INNER_PRODUCT2(ARG)\
   for (i = hg->vindex[candidate_gno]; i < hg->vindex[candidate_gno+1]; i++)  {\
     edge = hg->vedge[i];\
     if (((ARG)>0.0) && (hg->esize[edge]<hgp->MatchEdgeSizeThreshold)) {\
       for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++)  {\
         if (cmatch[hg->hvertex[j]] == hg->hvertex[j])    {\
           if (sums[hg->hvertex[j]] == 0.0)\
             index[m++] = hg->hvertex[j];\
           sums[hg->hvertex[j]] += (ARG);\
         }\
       }\
     }\
   }  



/* simple macro to start timer */
#define MACRO_TIMER_START(arg, message, sync) \
  if (hgp->use_timers > 3)  {\
    if (timer->matchstage[arg] < (arg))\
      timer->matchstage[arg] = Zoltan_Timer_Init(zz->ZTime, sync, message);\
    ZOLTAN_TIMER_START(zz->ZTime, timer->matchstage[arg], hg->comm->Communicator);\
  }   

/* simple corresponding macro to stop timer */
#define MACRO_TIMER_STOP(arg) \
  if (hgp->use_timers > 3) \
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->matchstage[arg], hg->comm->Communicator);

/* convenience macro to encapsulate resizing a buffer when necessary. Note: */
/* currently ZOLTAN_REALLOC aborts on any error and doesn't return - But... */
#define MACRO_REALLOC(new_size, old_size, buffer)  {\
  old_size = (new_size); \
  if (!(buffer = (char*) ZOLTAN_REALLOC (buffer, (old_size) * sizeof(char) ))) \
    MEMORY_ERROR; \
  } 

/* instead of realloc; just free and alloc new one to avoid memcpy in realloc */
#define MACRO_RESIZE(new_size, old_size, buffer)  {\
  if ((new_size)>(old_size)) {\
    old_size = (new_size);\
    ZOLTAN_FREE(&buffer);\
    if (!(buffer = (char*) ZOLTAN_MALLOC (old_size * sizeof(char)))) \
        MEMORY_ERROR; \
    } \
} 

/****************************************************************************/
/* Because this calculation is done in two locations it has been converted to a
** subroutine to assure consistancy. Inline is not yet always available!
** ERIK: need to calculate nCandidates based on # of unmatched vertices     */
static int calc_nCandidates (int num_vtx, int procs)
{
  /* Constant 2 below because each match pairs 2 vertices */
  return num_vtx ? 1 + num_vtx/(2 * procs * ROUNDS_CONSTANT) : 0;
}
 

/****************************************************************************/
static int pmatching_ipm (ZZ *zz,
  HGraph* hg,
  ZOLTAN_GNO_TYPE *match,
  PHGPartParams *hgp)
{
  int k, kstart, old_kstart, edge;
  int i, j = 0, n, m, round, vindex;                    /* loop counters  */
  int lno, bestlno, count = 0;                        /* temp variables */
  int nRounds;                /* # of matching rounds to be performed;       */
                              /* identical on all procs in hgc->Communicator.*/

  int nCandidates;            /* # of candidates on this proc; identical     */
                                          /* on all procs in hgc->col_comm.              */
  int total_nCandidates;      /* Sum of nCandidates across row. */
  int candidate_index = 0, first_candidate_index = 0;
  int nTotal;
  MPI_Datatype zoltan_gno_mpi_type;

  ZOLTAN_GNO_TYPE candidate_gno = 0;             /* gno of current candidate */

  int *dest = NULL,    nDest,  
      *size = NULL,    nSize,
      *index = NULL,   nIndex,
      *aux = NULL;

  int *visit = NULL,      /* candidate visit order (all candidates) */
      *permute = NULL,    /* reorder of candidates after global communication */
      *select = NULL;

  intptr_t *gno_locs=NULL;

  ZOLTAN_GNO_TYPE  partner_gno;
  ZOLTAN_GNO_TYPE *cmatch = NULL;  /* working copy of match array */

  char *sendbuf=NULL, *recvbuf=NULL, *s, *r, *edgebuf = NULL;
  ZOLTAN_GNO_TYPE *gnoptr=NULL, *gno=NULL;
  int *intptr=NULL;
  float *floatptr=NULL;
  int gno_size, int_size, float_size;


  /* TODO64 These seven variables are ints, but it seems likely that if the global
   * number of vertices or edges requires a 64 bit integer, that these counts, sizes,
   * etc. may at some point also need those storage sizes.  But they are passed to MPI
   * and Zoltan_Comm_*, which both expect ints.
   */

  int nRec, nSend, reccnt=0, sendcnt, recsize, sendsize, msgsize, nEdgebuf;

  float f, bestsum;      /* holds current best inner product */
  float *sums = NULL; /* holds candidate's inner products with each local vtx */
  PHGComm *hgc = hg->comm;
  int ierr = ZOLTAN_OK;
  int max_nPins, max_nVtx;       /* Global max # pins/proc and vtx/proc */
  intptr_t *rows = NULL;              /* used only in merging process */
  Triplet *master_data = NULL, *global_best = NULL;
  int *master_procs = NULL;
  int cFLAG;                    /* if set, do only a column matching, c-ipm */
  MPI_Op phasethreeop;
  MPI_Datatype phasethreetype;
  int pref = 0, num_matches_considered = 0;
  double ipsum = 0.;
  struct phg_timer_indices *timer = Zoltan_PHG_LB_Data_timers(zz);
  char *yo = "pmatching_ipm";
  
   
  ZOLTAN_TRACE_ENTER (zz, yo);
  MACRO_TIMER_START (0, "matching setup", 0);
  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();
  Zoltan_Srand_Sync (Zoltan_Rand(NULL), &(hgc->RNGState_col), hgc->col_comm);
     
  /* set a flag if user wants a column matching or a full matching */
  cFLAG = strcasecmp (hgp->redm_str, "c-ipm") ? 0 : 1;
  if (!cFLAG) {
    MPI_Type_contiguous (sizeof(Triplet), MPI_CHAR, &phasethreetype);
    MPI_Type_commit (&phasethreetype);
    MPI_Op_create (&phasethreereduce, 1, &phasethreeop);
  }

  /* determine basic working parameters */
  nRounds     = cFLAG ? ROUNDS_CONSTANT : hgc->nProc_x * ROUNDS_CONSTANT;
  nCandidates = calc_nCandidates (hg->nVtx, cFLAG ? 1 : hgc->nProc_x); 

  /* determine maximum number of Vtx and Pins for storage allocation */
  /* determine initial sum of all candidates = total_nCandidates==>allocation */
  if (cFLAG)  {
    max_nVtx          = hg->nVtx;
    max_nPins         = hg->nPins;
    total_nCandidates = nCandidates;    
  }
  else  {
    MPI_Allreduce(&hg->nPins, &max_nPins, 1, MPI_INT,MPI_MAX,hgc->Communicator);
    max_nVtx = 0;
    total_nCandidates = 0;
    for (i = 0; i < hgc->nProc_x; i++)  {
      count = (int)(hg->dist_x[i+1] - hg->dist_x[i]); /* number of vertices on proc i */
      if (count > max_nVtx)
        max_nVtx = count;
      if (i == hgc->myProc_x)
        first_candidate_index = total_nCandidates;
      total_nCandidates += calc_nCandidates (count, hgc->nProc_x);
    }
  }

  /* allocate "complicated" fixed-sized array storage */
  msgsize = MAX (total_nCandidates, max_nVtx);
  nIndex = 1 + MAX (msgsize, MAX (hgc->nProc_x, hgc->nProc_y));
  nDest  = nIndex;
  nSize  = nIndex;

  nSend = nEdgebuf = nRec = max_nPins * sizeof(int);   /* nSend/nRec for all other paired communication */

  if (hg->nVtx)  
    if (!(cmatch = (ZOLTAN_GNO_TYPE*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(ZOLTAN_GNO_TYPE)))
     || !(visit  = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
     || !(aux    = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))     
     || !(sums   = (float*) ZOLTAN_CALLOC (hg->nVtx,  sizeof(float))))
        MEMORY_ERROR;

  if (!cFLAG && total_nCandidates && (hgc->myProc_y == 0)) {  /* Master row */
    if (!(master_data=(Triplet*)ZOLTAN_CALLOC(total_nCandidates,sizeof(Triplet)))
     || !(global_best=(Triplet*)ZOLTAN_CALLOC(total_nCandidates,sizeof(Triplet))))
        MEMORY_ERROR;
    for (i = 0; i < total_nCandidates; i++) {
      master_data[i].candidate = -1;
      master_data[i].partner   = -1;
      master_data[i].ip        = -1.0;
    }
  } 

  if (!cFLAG)
    if (!(edgebuf = (char*) ZOLTAN_MALLOC (nEdgebuf)))
        MEMORY_ERROR;
 

  if (!(dest    = (int*) ZOLTAN_MALLOC (nDest              * sizeof(int)))
   || !(size    = (int*) ZOLTAN_MALLOC (nSize              * sizeof(int)))
   || !(index   = (int*) ZOLTAN_MALLOC (nIndex             * sizeof(int)))
   || !(rows    = (intptr_t*) ZOLTAN_MALLOC ((hgc->nProc_y + 1) * sizeof(intptr_t)))
   || (nSend && !(sendbuf  = (char *)ZOLTAN_MALLOC (nSend)))
   || (nRec  && !(recvbuf  = (char *)ZOLTAN_MALLOC (nRec)))
   || (total_nCandidates && !(permute = (int*) ZOLTAN_MALLOC (total_nCandidates * sizeof(int))))
   || (total_nCandidates && !(gno_locs = (intptr_t *) ZOLTAN_MALLOC (total_nCandidates * sizeof(intptr_t))))
   || (total_nCandidates && !(select = (int*) ZOLTAN_MALLOC (total_nCandidates * sizeof(int)))))
      MEMORY_ERROR;
  
  /* Compute candidates' vertex visit order (selection). Random is default. */
  Zoltan_PHG_Vertex_Visit_Order (zz, hg, hgp, visit);

  
  /* Loop processing ncandidates vertices per column each round.
   * Each loop has 3 phases, phase 3 may be repeated as necessary
   * Phase 1: send ncandidates vertices for global matching - horizontal comm
   * Phase 2: sum  inner products, find best in column - vertical communication
   * Phase 3: return best match selections        - horizontal communication
   *
   * No conflict resolution is required since temp locking prevents conflicts. */

  MACRO_TIMER_STOP (0);
  vindex = 0;                        /* marks current position in visit array */
  for (round = 0; round < nRounds; round++) {
    MACRO_TIMER_START (1, "matching phase 1", 0);
    
    /************************ PHASE 1: ****************************************/
    
    memcpy (cmatch, match, hg->nVtx * sizeof(ZOLTAN_GNO_TYPE));  /* for temporary locking */    
    if (cFLAG)  {
      /* select upto nCandidates unmatched vertices to locally match */
      for (i = 0; i < total_nCandidates; i++)
        permute[i] = -1;                       /* to flag missing candidates  */
        
      for (nTotal = 0; nTotal < nCandidates && vindex < hg->nVtx; vindex++)
        if (cmatch [visit[vindex]] == (ZOLTAN_GNO_TYPE)visit[vindex])  {         /* unmatched */
          permute [nTotal++] = visit[vindex];    /* select it as a candidate */
          cmatch [visit[vindex]] = (ZOLTAN_GNO_TYPE)(-visit[vindex])-1;  /* mark as pending match */
        }

      total_nCandidates = nTotal;
      for (i = 0; i < total_nCandidates; i++)
        select[i] = i;      
    }
    else  {       
      /* Select upto nCandidates unmatched vertices to globally match. */
      for (sendcnt = 0; sendcnt < nCandidates && vindex < hg->nVtx; vindex++)
        if (cmatch[visit[vindex]] == (ZOLTAN_GNO_TYPE)visit[vindex])  {         /* unmatched */
          select [sendcnt++] = visit[vindex];    /* select it as a candidate */
          cmatch [visit[vindex]] = (ZOLTAN_GNO_TYPE)(-visit[vindex])-1;  /* mark as pending match */
        }
                        
      /* assure send buffer is large enough by first computing required size */

      gno_size = int_size = float_size = 0;

      for (i = 0; i < sendcnt; i++)  {
        lno = select[i];
        n = hg->vindex[lno+1] - hg->vindex[lno];
        if (n > 0){
            gno_size += sizeof(ZOLTAN_GNO_TYPE);         /* gno of candidate */
            int_size += (2 + n) * sizeof(int);   /* candidate index, edge count, edge lnos */
            if (hgp->UsePrefPart) int_size += sizeof(int);    /* preferred partition */
        }
      }

      sendsize = gno_size + int_size;     /* total size in chars */

      if (sendsize > nSend)
        MACRO_REALLOC (1.2 * sendsize, nSend, sendbuf);    /* resize send buffer */    

      /* put <candidate_gno, candidate_index, count, <edge>> into send buffer */

      s = sendbuf;

      for (i = 0; i < sendcnt; i++)   {

        lno = select[i];
        n = hg->vindex[lno+1] - hg->vindex[lno];

        if (n > 0){
          gnoptr = (ZOLTAN_GNO_TYPE *)s;
          intptr = (int *)(gnoptr + 1);

          *gnoptr++ = VTX_LNO_TO_GNO(hg, lno);                  /* gno of candidate */
          *intptr++ = i + first_candidate_index;                 /* candidate index */

          if (hgp->UsePrefPart)          
              *intptr++ = hg->pref_part[lno];                /* pref partition info */
          *intptr++ = n;                                              /* edge count */
          for (j = hg->vindex[lno]; j < hg->vindex[lno+1]; j++)  
            *intptr++ = hg->vedge[j];                                   /* edge lno */

          s = (char *)intptr;
        }
      }
         
      /* communication to determine global size of rec buffer */
      MPI_Allgather (&sendsize, 1, MPI_INT, size, 1, MPI_INT, hgc->row_comm);
     
      /* determine size of the rec buffer & reallocate bigger iff necessary */
      recsize = 0;
      for (i = 0; i < hgc->nProc_x; i++)
        recsize += size[i];          /* compute total size of edgebuf in chars */
      if (recsize > nEdgebuf)
        MACRO_REALLOC (1.2 * recsize, nEdgebuf, edgebuf);  /* enlarge edgebuf */
    
      /* setup displacement array necessary for MPI_Allgatherv */
      dest[0] = 0;
      for (i = 1; i < hgc->nProc_x; i++)
        dest[i] = dest[i-1] + size[i-1];

      /* communicate vertices & their edges to all row neighbors */
      MPI_Allgatherv(sendbuf, sendsize, MPI_CHAR, edgebuf, size, dest, MPI_CHAR, hgc->row_comm);
         
      /* Communication has grouped candidates by processor, rescramble!     */
      /* Otherwise all candidates from proc column 0 will be matched first, */

      for (i = 0; i < total_nCandidates; i++)
        select[i] = i;
      Zoltan_Rand_Perm_Int (select, total_nCandidates, &(hgc->RNGState_col));

      for (i = 0; i < total_nCandidates; i++)
        permute[i] = -1;                 /* to flag missing sparse candidates */

      r = edgebuf;
      n = 0;

      while (r < edgebuf + recsize){

        gnoptr = (ZOLTAN_GNO_TYPE *)r;
        intptr = (int *)(gnoptr + 1);

        gno_locs[n] = (char *)gnoptr - edgebuf;

        candidate_index = *intptr++;

        permute[candidate_index] = n;

        if (hgp->UsePrefPart) intptr++;       
        count           =  *intptr++;       /* count of edges */      

        r = (char *)(intptr + count);  /* skip of edge lnos */
        n++;
      }
    }            /* DONE:  if (cFLAG) else ...  */                          


    MACRO_TIMER_STOP (1);
    
    /************************ PHASE 2: ***************************************/
      
    /* for each candidate vertex, compute all local partial inner products */    


    kstart = old_kstart = 0;         /* next candidate (of nTotal) to process */
    while (kstart < total_nCandidates) {

      MACRO_TIMER_START (2, "Matching kstart A", 0);
      for (i = 0; i < hgc->nProc_y; i++)
        rows[i] = -1;                  /* to flag data not found for that row */
      sendsize = 0;                    /* position in send buffer */
      sendcnt  = 0;                    /* count of messages in send buffer */
      s = sendbuf;

      for (k = kstart; k < total_nCandidates; k++)  {

        n = permute[select[k]];

        if (n == -1) 
          continue;                /* don't have this sparse candidate locally */
        
        if (!cFLAG)  {

          gnoptr = (ZOLTAN_GNO_TYPE *)(edgebuf + gno_locs[n]);
          intptr = (int *)(gnoptr + 1);

          candidate_gno = *gnoptr;

          candidate_index = *intptr++;          /* candidate_index of vertex */
          if (hgp->UsePrefPart)          
              pref = *intptr++;                 /* pref vertex information */          
          count           = *intptr++;          /* count of following hyperedges */

        }
        else  {
          candidate_index = k;
          candidate_gno   = permute[k];  /* need to use next local vertex */

          if (hgp->UsePrefPart)          
              pref = hg->pref_part[candidate_gno];          
        }                          /* here candidate_gno is really a local id */

                  
        /* now compute the row's nVtx inner products for kth candidate */
        m = 0;
        if (!cFLAG) {
          if ((hg->ewgt != NULL) && (hgp->vtx_scal == NULL)){
            INNER_PRODUCT1(hg->ewgt[*intptr])

          }
          else if ((hg->ewgt == NULL) && (hgp->vtx_scal == NULL)){
            INNER_PRODUCT1(1.0)
          }
          else if ((hg->ewgt != NULL) && (hgp->vtx_scal != NULL)){
            INNER_PRODUCT1(hgp->vtx_scal[hg->hvertex[j]] * hg->ewgt[*intptr])
          }
          else {/* UVC: no need: if ((hg->ewgt == NULL) && (hgp->vtx_scal != NULL)) */
            INNER_PRODUCT1(hgp->vtx_scal[hg->hvertex[j]])
          }
        } else   {                                            /* cFLAG */
          if      ((hg->ewgt == NULL) && (hgp->vtx_scal == NULL))
            INNER_PRODUCT2(1.0)
          else if ((hg->ewgt == NULL) && (hgp->vtx_scal != NULL))
            INNER_PRODUCT2(hgp->vtx_scal[hg->hvertex[j]])
          else if ((hg->ewgt != NULL) && (hgp->vtx_scal == NULL))
            INNER_PRODUCT2(hg->ewgt[edge])
          else if ((hg->ewgt != NULL) && (hgp->vtx_scal != NULL))
            INNER_PRODUCT2(hgp->vtx_scal[hg->hvertex[j]] * hg->ewgt[edge])
        }
          
        /* if local vtx, remove self inner product (useless maximum) */
        if (cFLAG)
          sums[candidate_gno] = 0.0;     /* since candidate_gno is really lno */
        else if (VTX_TO_PROC_X(hg, candidate_gno) == hgc->myProc_x)
          sums[VTX_GNO_TO_LNO(hg, candidate_gno)] = 0.0;
         
        /* count partial sums exceeding PSUM_THRESHOLD */   
        count = 0;
        for (i = 0; i < m; i++)  {
          lno = index[i];
          if (sums[lno] > PSUM_THRESHOLD
              && MATCH_OK(hgp, hg, hg->pref_part[lno], pref))
            aux[count++] = lno;      /* save lno for significant partial sum */
          else
            sums[lno] = 0.0;         /* clear unwanted entries */  
        }     
        if (count == 0)
          continue;         /* no partial sums to append to message */       

        /* iff necessary, resize send buffer to fit at least first message */

        gno_size = sizeof(ZOLTAN_GNO_TYPE);
        int_size = (count + 2) * sizeof(int);
        float_size = count * sizeof(float);

        msgsize = gno_size + int_size + float_size;

        if (sendcnt == 0 && (msgsize > nSend)) {
          MACRO_REALLOC (1.2 * msgsize, nSend, sendbuf);  /* increase buffer size */
          s = sendbuf;         
        }
        
        /* message is <candidate_gno, candidate_index, count, <lno>, <psum>> */
        if (sendsize + msgsize <= nSend)  {
          /* flag first data in each destination row for merging */
          if (rows[candidate_gno % hgc->nProc_y] != 1)  {
            rows[candidate_gno % hgc->nProc_y] = 1;
            candidate_index = -candidate_index -1;
          }
          
          /* current partial sums fit, so put them into the send buffer */
          dest[sendcnt]   = candidate_gno % hgc->nProc_y;    /* destination */
          size[sendcnt++] = msgsize;          /* size of message */

          s = sendbuf + sendsize;
          gnoptr = (ZOLTAN_GNO_TYPE *)(s);
          s += gno_size;
          intptr = (int *)(s);
          s += int_size;
          floatptr = (float*)(s);
  
          *gnoptr++ = candidate_gno;
          *intptr++ = candidate_index;        
          *intptr++ = count;
          for (i = 0; i < count; i++)  {          
            *intptr++ = aux[i];                   /* lno of partial sum */
          }
          for (i = 0; i < count; i++)  {          
            *floatptr++ = sums[aux[i]];           /* partial sum */
            sums[aux[i]] = 0.0;
          }
          sendsize       += msgsize;          /* cummulative size of message */
        }
        else  {           /* psum message doesn't fit into buffer */
          for (i = 0; i < count; i++)              
            sums[aux[i]] = 0.0;
          break;   
        }  
      }                  /* DONE: loop over k */                    

      MACRO_TIMER_STOP (2);
      MACRO_TIMER_START (3, "Matching kstart B", 0);
    
      /* synchronize all rows in this column to next kstart value */
      old_kstart = kstart;
      MPI_Allreduce (&k, &kstart, 1, MPI_INT, MPI_MIN, hgc->col_comm);

      /* Send inner product data in send buffer to appropriate rows */
      ierr = communication_by_plan (zz, sendcnt, dest, size, 1, sendbuf, &reccnt, 
       &recsize, &nRec, &recvbuf, hgc->col_comm, IPM_TAG);

      if (ierr != ZOLTAN_OK)
        goto End;
    
      /* build index into receive buffer pointer for each proc's row of data */
      for (i = 0; i < hgc->nProc_y; i++)
        rows[i] = recsize;       /* sentinal in case no row's data available */

      k = 0;
      r = recvbuf;

      while(r < recvbuf+recsize){
        gnoptr = (ZOLTAN_GNO_TYPE *)r;
        intptr = (int *)(gnoptr + 1);
        count = intptr[1];
        candidate_index = intptr[0];
        if (candidate_index < 0){                    /* candidate_index */
          *intptr = -candidate_index - 1;         /* make sentinal candidate_index positive */
          rows[k++] = (char *)gnoptr - recvbuf;   /* so we can find gno later */
        }

        r += sizeof(ZOLTAN_GNO_TYPE) + ((2 + count) * sizeof(int)) + (count * sizeof(float));
      }

      /* merge partial i.p. sum data to compute total inner products */
      s = sendbuf; 
      for (n = old_kstart; n < kstart; n++) {
        m = 0;       
        for (i = 0; i < hgc->nProc_y; i++) {

          if (rows[i] < recsize){
            gnoptr = (ZOLTAN_GNO_TYPE *)(recvbuf + rows[i]);
            intptr = (int *)(gnoptr + 1);

            if (intptr[0] == select[n])  {   /* candidate index */
              candidate_gno   = gnoptr[0];
              candidate_index = intptr[0];
              count = intptr[1];
              floatptr = (float *)(intptr + 2 + count);

              for (j = 0; j < count; j++)  {
                lno = intptr[2 + j];
                if (sums[lno] == 0.0){       /* is this first time for this lno? */
                  aux[m++] = lno;           /* then save the lno */
                }
  
                sums[lno] += *floatptr++;
              }
              rows[i] = (char *)floatptr - recvbuf;
            }
          }
        }

        /* determine how many total inner products exceed threshold */
        count = 0;
        for (i = 0; i < m; i++)
          if (sums[aux[i]] > TSUM_THRESHOLD)
            count++;

        /* Put <candidate_gno, candidate_index, count, <lno>, <tsum>> into send */
        msgsize = 0;

        if (count > 0)  {

          msgsize = sizeof(ZOLTAN_GNO_TYPE) + ((2 + count) * sizeof(int)) + (count * sizeof(float));

          if (s - sendbuf + msgsize > nSend ) {
            sendsize = s - sendbuf;
            MACRO_REALLOC (1.2*(sendsize + msgsize), nSend, sendbuf);
            s = sendbuf + sendsize;         /* since realloc buffer could move */
          }      

          gnoptr = (ZOLTAN_GNO_TYPE *)s;
          intptr = (int *)(gnoptr + 1);
          floatptr = (float *)(intptr + 2 + count);

          *gnoptr = candidate_gno;
          *intptr++ = candidate_index;
          *intptr++ = count;
        }  
        for (i = 0; i < m; i++)   {
          lno = aux[i];             
          if (sums[lno] > TSUM_THRESHOLD)  {
            *intptr++ = lno;
            *floatptr++ = sums[lno];
          }  
          sums[lno] = 0.0;  
        }     
        s += msgsize;
      }
      sendsize = s - sendbuf;
    
      /* Communicate total inner product results to MASTER ROW */
      MPI_Gather(&sendsize, 1, MPI_INT, size, 1, MPI_INT, 0, hgc->col_comm);
    
      if (hgc->myProc_y == 0) {
        recsize = 0;
        for (i = 0; i < hgc->nProc_y; i++)
          recsize += size[i];        
        
        dest[0] = 0;
        for (i = 1; i < hgc->nProc_y; i++)
          dest[i] = dest[i-1] + size[i-1];
        
        if (recsize > nRec)
          MACRO_REALLOC (1.2 * recsize, nRec, recvbuf);      /* make rec buffer bigger */
      }

      MPI_Gatherv (sendbuf, sendsize, MPI_CHAR, recvbuf, size, dest, MPI_CHAR, 0, hgc->col_comm);

      /* Determine best vertex and best sum for each candidate */
      if (hgc->myProc_y == 0) {   /* do following only if I am the MASTER ROW */
        for (r = recvbuf; r < recvbuf + recsize;)  {
          gnoptr = (ZOLTAN_GNO_TYPE *)r;
          intptr = (int *)(gnoptr + 1);

          candidate_gno   = *gnoptr++;
          candidate_index = *intptr++;
          count           = *intptr++;                    /* count of nonzero pairs */
          bestsum = -1.0;                        /* any negative value will do */
          bestlno = -1;                          /* any negative value will do */

          floatptr = (float *)(intptr + count);

          for (i = 0; i < count; i++)  {
            lno =  *intptr++;
            f   =  *floatptr++;     

            if ((f > bestsum) && cmatch[lno] == lno)  { 
              bestsum = f;
              bestlno = lno;
            }      
          }

          r = (char *)(floatptr);
         
          /* For hybrid ipm, keep matches that are above average in c-ipm */
          if (bestsum>0){
             /* ipsum is cumulative sum of best inner products (bestsum) */
             ipsum += bestsum;
             num_matches_considered++;        
          }
          if (cFLAG && bestsum > MAX(TSUM_THRESHOLD, 
              hgp->hybrid_keep_factor*ipsum/num_matches_considered))  {            
            cmatch[bestlno] = -1;                   
            match[bestlno]       = candidate_gno;
            match[candidate_gno] = bestlno;
          }
                        
          if (!cFLAG && bestsum > TSUM_THRESHOLD)  {
            cmatch[bestlno] = -1;  /* mark pending match to avoid conflicts */
            master_data[candidate_index].candidate = candidate_gno;
            master_data[candidate_index].partner = VTX_LNO_TO_GNO (hg, bestlno);
            master_data[candidate_index].ip = bestsum;
          }
        }
      } 
      MACRO_TIMER_STOP (3);    
    }            /* DONE: kstart < max_nTotal loop */

    if (cFLAG)  {
      MPI_Bcast (match, hg->nVtx, zoltan_gno_mpi_type, 0, hgc->col_comm);          
      continue;      /* skip phases 3 and 4, continue rounds */ 
    }    

    /************************ NEW PHASE 3: ********************************/
    
    MACRO_TIMER_START (4, "Matching Phase 3", 1);
    
    /* Only MASTER ROW computes best global match for candidates */
    /* EBEB or perhaps we can do this fully distributed? */
    if (hgc->myProc_y == 0) {
      HG_Ptr = hg;
      MPI_Allreduce(master_data, global_best, total_nCandidates, phasethreetype,
                    phasethreeop, hgc->row_comm);

      /* Look through array of "winners" and update match array. */
      /* Local numbers are used for local matches, otherwise
         -(gno+1) is used in the match array.                    */
      for (i = 0; i < total_nCandidates; i++) {
        int cproc, vproc;
        candidate_gno = global_best[i].candidate;

        /* Reinitialize master_data for next round */
        master_data[i].candidate = -1;
        master_data[i].partner   = -1;
        master_data[i].ip        = -1.0;
        if (candidate_gno == -1)
          continue;

        partner_gno = global_best[i].partner;
        cproc = VTX_TO_PROC_X(hg, candidate_gno);
        vproc = VTX_TO_PROC_X(hg, partner_gno);
        if (cproc == hgc->myProc_x) {
          if (vproc == hgc->myProc_x)   {
            int v1 = VTX_GNO_TO_LNO(hg, partner_gno);
            int v2 = VTX_GNO_TO_LNO(hg, candidate_gno);
            match[v1] = (ZOLTAN_GNO_TYPE)v2;
            match[v2] = (ZOLTAN_GNO_TYPE)v1;
          }
          else 
            match[VTX_GNO_TO_LNO(hg, candidate_gno)] = -partner_gno - 1;
        }                         
        else if (vproc == hgc->myProc_x)
          match[VTX_GNO_TO_LNO(hg, partner_gno)] = -candidate_gno - 1;
      }
    } /* End (hgc->myProc_y == 0) */

    /* broadcast match array to the entire column */
    MPI_Bcast (match, hg->nVtx, zoltan_gno_mpi_type, 0, hgc->col_comm);
    MACRO_TIMER_STOP (4);                       /* end of phase 3 */
  }                                             /* DONE: loop over rounds */

  ZOLTAN_FREE(&aux);
  ZOLTAN_FREE(&sums);
  ZOLTAN_FREE(&gno_locs);
  ZOLTAN_FREE(&permute);
  ZOLTAN_FREE(&select);
  ZOLTAN_FREE(&edgebuf);
  ZOLTAN_FREE(&rows);
  ZOLTAN_FREE(&index);
  ZOLTAN_FREE(&visit);
  ZOLTAN_FREE(&global_best);
  ZOLTAN_FREE(&master_data);
  
  MACRO_TIMER_START (6, "Matching Cleanup", 0);

  /* optional sanity tests */
  if (zz->Debug_Level > 4 && hgc->myProc_x == 0 && hgc->myProc_y == 0)  {
    int local = 0, global = 0, unmatched = 0;
    for (i = 0; i < hg->nVtx; i++)  {
      if      (match[i] == i)  unmatched++;
      else if (match[i] < 0)   global++;
      else                     local++;
    }
    uprintf (hgc, "%d RTHRTH %d unmatched, %d external, %d local of %d\n",
     hg->info, unmatched, global, local, hg->nVtx);
  }

  if (zz->Debug_Level > 4 && hgc->myProc_x==0 && hgc->myProc_y==0)
    fprintf (stdout, "RTHRTH rounds %d\n", nRounds);

  if (zz->Debug_Level > 4)  {
    /* The following tests that the global match array is a valid permutation */
    /* NOTE:  THESE TESTS ARE NOT MANDATORY; THEY CAN BE EXCLUDED AFTER WE    */
    /* COMPLETE TESTING OF matching_ipm.                                      */

    for (i = 0; i < hg->nVtx; i++)
      if (match[i] < 0)  cmatch[i] = -match[i] - 1;
      else               cmatch[i] = VTX_LNO_TO_GNO(hg, match[i]);

    MPI_Allgather (&hg->nVtx, 1, MPI_INT, size, 1, MPI_INT, hgc->row_comm); 

    recsize = 0;
    for (i = 0; i < hgc->nProc_x; i++)
      recsize += size[i];

    msgsize = recsize * sizeof(ZOLTAN_GNO_TYPE);
  
    dest[0] = 0;
    for (i = 1; i < hgc->nProc_x; i++)
      dest[i] = dest[i-1] + size[i-1];
  
    if (nRec < msgsize){
      MACRO_REALLOC (msgsize, nRec, recvbuf);
    }

    gno = (ZOLTAN_GNO_TYPE *)recvbuf;

    MPI_Allgatherv (cmatch, hg->nVtx, zoltan_gno_mpi_type, recvbuf, size, dest, zoltan_gno_mpi_type,
     hgc->row_comm);

    if (nSend < msgsize){
      MACRO_REALLOC (msgsize, nSend, sendbuf);  /* make send buffer bigger */
    }

    gnoptr = (ZOLTAN_GNO_TYPE *)sendbuf;

    /* TODO64  Here we are indexing an array by a global number.  Is this a problem? */
    /*  From 5/4/10 meeting: return an error in phg build if it looks like a proc could have
        more than 2*10^9 vertices in the 2D distribution 
        (so globalNumVertex/sqrt(numProcs) < 2*10^9) */

    for (i = 0; i < recsize; i++)
      gnoptr[i] = 0;
    for (i = 0; i < recsize; i++)
      ++gnoptr[gno[i]];

    count = 0;
    for (i = 0; i < recsize; i++)
      if (gnoptr[i] != 1)
        count++;
    if (count)    
      uprintf (hgc, "RTHRTH %d FINAL MATCH ERRORS of %d\n", count, recsize); 
  }
  MACRO_TIMER_STOP (6);

End:
  if (!cFLAG) {
    MPI_Op_free(&phasethreeop);
    MPI_Type_free(&phasethreetype);
    ZOLTAN_FREE(&global_best);
    ZOLTAN_FREE(&master_data);
  }

  Zoltan_Multifree (__FILE__, __LINE__, 15, &cmatch, &visit, &sums, &sendbuf,
   &dest, &size, &recvbuf, &index, &aux, &permute, &edgebuf, &select, &rows,
   &gno_locs, &master_procs);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}


  
/****************************************************************************/

static int communication_by_plan (ZZ* zz, int sendcnt, int* dest, int* size, 
 int scale, char * send, 
 int* reccnt, int *recsize, int * nRec, char ** rec, MPI_Comm comm, int tag)
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
   else
     *recsize = *reccnt * scale;
   
   /* realloc rec buffer if necessary */  
   if (*recsize > *nRec)  {   
     if (!(*rec = (char*) ZOLTAN_REALLOC (*rec, *recsize * sizeof(char)))){
       ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error");
       return ZOLTAN_MEMERR;
     }
     *nRec = *recsize;
   }
   
   /* send messages from send buffer to destinations */      
   err = Zoltan_Comm_Do (plan, tag+2, (char*) send, scale * sizeof(char), (char*) *rec);
   if (err != ZOLTAN_OK)  {
     ZOLTAN_PRINT_ERROR (zz->Proc, yo, "failed in Comm_Do");
     return err;
   }
   
   /* free memory associated with the plan */
   Zoltan_Comm_Destroy (&plan); 
   return ZOLTAN_OK;
}




/****************************************************************************/
static int pmatching_agg_ipm (ZZ *zz,
                              HGraph* hg,
                              ZOLTAN_GNO_TYPE *match,
                              PHGPartParams *hgp)
{
  int ierr = ZOLTAN_OK;    
  int i, j, n, m, round, vindex;                        /* loop counters  */
  int sendcnt, sendsize, reccnt=0, recsize, msgsize;       /* temp variables */
  int nRounds;                /* # of matching rounds to be performed;       */
  int nSend,             /* working buffers and their sizes */
    *dest = NULL,    nDest,  
    *size = NULL,    nSize,
     nRec,
    *aux = NULL,   
     nEdgebuf;  /* holds candidates for processing (ipm)   */
  int *visit = NULL,       /* candidate visit order (all candidates) */
    *lhead = NULL,       /* to accumulate ipm values correctly */
    *lheadpref = NULL,
    *idxptr = NULL;    /* reorder of candidates after global communication */
  float bestsum;      /* holds current best inner product */
  float *sums = NULL, /* holds candidate's inner products with each local vtx */
    *cw = NULL,   /* current vertex weight */
    *tw=NULL, *maxw = NULL, *candw = NULL,
    *fptr;
  char *visited = NULL;
  PHGComm *hgc = hg->comm;
  int max_nPins, max_nVtx;       /* Global max # pins/proc and vtx/proc */
  intptr_t *rows = NULL;              /* used only in merging process */
  int bestlno, lno;
  Triplet *master_data = NULL, *global_best = NULL;
  int *master_procs = NULL;
  MPI_Op phasethreeop;
  MPI_Datatype phasethreetype;
  int VtxDim = (hg->VtxWeightDim>0) ? hg->VtxWeightDim : 1;
  int pref = 0;
  int replycnt, header_size;
  struct phg_timer_indices *timer = Zoltan_PHG_LB_Data_timers(zz);
  char *yo = "pmatching_agg_ipm";
  KVHash hash;

  ZOLTAN_GNO_TYPE candidate_gno = 0;         /* gno of current candidate */
                                             /* identical on all procs in hgc->Communicator.*/
  int nCandidates;            /* # of candidates on this proc; identical     */
  /* on all procs in hgc->col_comm.              */
  int total_nCandidates = 0;      /* Sum of nCandidates across row. */
  int candidate_index = 0, *candIdx;
  int *locCandidates = NULL, locCandCnt,      /* current selected candidates (this round) & number */
                  *candvisit=NULL;                      /* to randomize visit order of candidates*/

  int count=0, k, kstart, old_kstart;                      /* temp variables */
  ZOLTAN_GNO_TYPE partner_gno;

  char *sendbuf=NULL, *recvbuf=NULL, *edgebuf=NULL;
  char *r, *s;

  ZOLTAN_GNO_TYPE *gnoptr;
  int *intptr;
  float *floatptr;
  int gno_size, int_size, float_size;

  intptr_t *gno_locs=NULL;

  ZOLTAN_TRACE_ENTER (zz, yo);
  MACRO_TIMER_START (0, "matching setup", 0);
  Zoltan_Srand_Sync (Zoltan_Rand(NULL), &(hgc->RNGState_col), hgc->col_comm);
  
  /* this restriction may be removed later, but for now NOTE this test */
  if (sizeof(int) < sizeof(float))  {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Code must be modified before using");
    ierr = ZOLTAN_FATAL;
    goto End;
  }
  
  
  MPI_Type_contiguous (sizeof(Triplet), MPI_CHAR, &phasethreetype);
  MPI_Type_commit (&phasethreetype);
  MPI_Op_create (&phasethreereduce, 1, &phasethreeop);
  
  /* compute nCandidates per process;
     unless hypergraph is too small use user given parameter
     otherwise use nVtx/3 to allow some matching
  */
  nCandidates = MAX(1, MIN(hgp->nCand, hg->nVtx/3));
  if (!(candIdx = (int*) ZOLTAN_MALLOC ((1+hgc->nProc_x) * sizeof(int)))) 
    MEMORY_ERROR;

    
  /* determine maximum number of Vtx and Pins for storage allocation */
  /* determine initial sum of all candidates = total_nCandidates==>allocation */
  MPI_Allreduce(&hg->nPins, &max_nPins, 1, MPI_INT,MPI_MAX,hgc->Communicator);
  max_nVtx = 0;
  total_nCandidates = 0;
  for (i = 0; i < hgc->nProc_x; i++)  {
    count = (int)(hg->dist_x[i+1] - hg->dist_x[i]); /* number of vertices on proc i TODO: good test for phg build*/
    if (count > max_nVtx)
      max_nVtx = count;
    candIdx[i] = total_nCandidates;
    total_nCandidates += MIN(hgp->nCand, count);
  }
  candIdx[i] = total_nCandidates;
                 
  /* allocate "complicated" fixed-sized array storage */
  msgsize = MAX (total_nCandidates, max_nVtx);
  nSize = nDest = 1 + MAX (msgsize, MAX (hgc->nProc_x, hgc->nProc_y));

  header_size = 2 + (hgp->UsePrefPart ? 1 : 0);

  nSend = nRec = nEdgebuf = 
    max_nPins * sizeof(int) + 
    total_nCandidates * (sizeof(ZOLTAN_GNO_TYPE) + (header_size * sizeof(int)) + (VtxDim * sizeof(float)));

  if (hg->nVtx)  
    if (!(lhead  = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
        || (hgp->UsePrefPart && !(lheadpref = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int))))
        || !(visited= (char*)  ZOLTAN_MALLOC (hg->nVtx * sizeof(char)))
        || !(cw     = (float*) ZOLTAN_MALLOC (VtxDim * hg->nVtx * sizeof(float)))
        || !(visit  = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
        || !(aux    = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))     
        || !(sums   = (float*) ZOLTAN_CALLOC (hg->nVtx,  sizeof(float))))
      MEMORY_ERROR;

  if (VtxDim) 
    if (!(tw = (float*) ZOLTAN_MALLOC (VtxDim * sizeof(float)))
        || !(maxw = (float*) ZOLTAN_MALLOC (VtxDim * sizeof(float))))
      MEMORY_ERROR;

  if (hgc->myProc_y==0 && total_nCandidates ) {  
    if (!(master_data=(Triplet*)ZOLTAN_CALLOC(total_nCandidates,sizeof(Triplet)))
     || !(global_best=(Triplet*)ZOLTAN_CALLOC(total_nCandidates,sizeof(Triplet))))
      MEMORY_ERROR;
    for (i = 0; i < total_nCandidates; i++) {
      master_data[i].candidate = -1;
      master_data[i].partner   = -1;
      master_data[i].ip        = -1.0;
    }
  } 

  if (!(edgebuf = (char*) ZOLTAN_MALLOC (nEdgebuf))
      || !(dest    = (int*) ZOLTAN_MALLOC (nDest              * sizeof(int)))
      || !(size    = (int*) ZOLTAN_MALLOC (nSize              * sizeof(int)))
      || !(rows    = (intptr_t*) ZOLTAN_MALLOC ((hgc->nProc_y + 1) * sizeof(intptr_t)))
      || (nRec  && !(recvbuf  = (char *) ZOLTAN_MALLOC (nRec)))
      || (nSend && !(sendbuf  = (char *) ZOLTAN_MALLOC (nSend))))
    MEMORY_ERROR;
  if (total_nCandidates) {
    if (!(idxptr   = (int*)   ZOLTAN_MALLOC (total_nCandidates * sizeof(int)))
        || !(candvisit = (int*)   ZOLTAN_MALLOC (total_nCandidates * sizeof(int)))                
        || !(gno_locs = (intptr_t*)   ZOLTAN_MALLOC (total_nCandidates * sizeof(intptr_t)))
        || !(locCandidates = (int*)   ZOLTAN_MALLOC (nCandidates * sizeof(int)))
        || !(candw  = (float*) ZOLTAN_MALLOC (VtxDim * total_nCandidates * sizeof(float))))
      MEMORY_ERROR;
    if (hgc->myProc_y==0) {
        if (Zoltan_KVHash_Create(&hash, 2*(1+hg->nVtx))==ZOLTAN_MEMERR)
            MEMORY_ERROR;
    }
  }
    
/* Compute candidates' vertex visit order (selection). Random is default. */
  Zoltan_PHG_Vertex_Visit_Order (zz, hg, hgp, visit);
  
/* Loop processing ncandidates vertices per column each round.
 * Each loop has 3 phases, phase 3 may be repeated as necessary
 * Phase 1: send ncandidates vertices for global matching - horizontal comm
 * Phase 2: sum  inner products, find best in column - vertical communication
 * Phase 3: return best match selections        - horizontal communication
 * Phase 4: finalize decision/conflict resolution */

  MACRO_TIMER_STOP (0);
  for (i=0; i< hg->nVtx; ++i) {
    visited[i] = 0;
    lhead[i] = i;
    if (hgp->UsePrefPart) lheadpref[i] = hg->pref_part[i];
    match[i] = VTX_LNO_TO_GNO(hg, i);
  }
  for (i=0; i<VtxDim; ++i)
    tw[i] = 0.0;
  if (hg->vwgt) {
    for (i=0; i< hg->nVtx; ++i)
      for (j=0; j<VtxDim; ++j)
        tw[j] += (cw[i*VtxDim+j] = hg->vwgt[i*VtxDim+j]);
  } else {
    tw[0] += hg->nVtx;
    for (i=0; i< hg->nVtx; ++i)                 
      cw[i] = 1.0;
  }
  MPI_Allreduce(tw, maxw, VtxDim, MPI_FLOAT, MPI_SUM, hgc->row_comm);
  for (j=0; j<VtxDim; ++j)
    maxw[j] /= 4.0; /* we don't allow a cluster to grow more than half of
                       a part */

                    
  ZOLTAN_FREE(&tw);

  vindex = 0;                        /* marks current position in visit array */
  nRounds = 1;
  for (round = 0; nRounds; ++round) {
    MACRO_TIMER_START (1, "matching phase 1", 0);
    
    /************************ PHASE 1: ****************************************/

    /* Select upto nCandidates unmatched vertices to globally match. */

    for (locCandCnt = 0; locCandCnt < nCandidates && vindex < hg->nVtx; ++vindex) {
      int v=visit[vindex];          
      if (!visited[v])  {         /* unmatched */
        locCandidates [locCandCnt++] = v;    /* select it as a candidate */
        visited[v] = 1;
        match[v] = -match[v] - 1;
      }
    }
/*    uprintf(hgc, "Round %d:   locCandCnt=%zd\n", round, locCandCnt); */

   /* Send buffer: offset-to-ints, offset-to-floats, all gnos, all ints, all floats */

    int_size = float_size = gno_size = 0;

    for (i = 0; i < locCandCnt; i++)  {
      lno = locCandidates[i];
      /* UVCUVC: CHECK if it is possible to use sparse communication
         current code only works if all candidates have been communicated */
/*      if (hg->vindex[lno+1] > hg->vindex[lno]) */
      int_size += hg->vindex[lno+1] - hg->vindex[lno];   /* edge lnos */
    }

    int_size += (locCandCnt * (2 + (hgp->UsePrefPart ? 1 : 0))); /* index, count, pref part */
    float_size += locCandCnt * VtxDim;                           /* weights */
    gno_size = locCandCnt;                                       /* candidate gnos */

    int_size *= sizeof(int);
    gno_size *= sizeof(ZOLTAN_GNO_TYPE);
    float_size *= sizeof(float);

    sendsize = gno_size + int_size + float_size;

    if (sendsize > nSend){
      MACRO_RESIZE (1.2 * sendsize, nSend, sendbuf);    /* resize send buffer */    
    }

    s = sendbuf;

    for (i = 0; i < locCandCnt; i++)   {
      lno = locCandidates[i];
/*      if (hg->vindex[lno+1] > hg->vindex[lno]) { */  /* UVCUVC CHECK sparse comm? */

        gnoptr = (ZOLTAN_GNO_TYPE *)s;
        intptr = (int *)(gnoptr + 1);

        *gnoptr = VTX_LNO_TO_GNO(hg, lno);                  /* gno of candidate */
        *intptr++ = candIdx[hgc->myProc_x] + i;               /* candidate index */
        
        if (hgp->UsePrefPart)          
          *intptr++ = hg->pref_part[lno];       /* pref partition info */

        *intptr++ = hg->vindex[lno+1] - hg->vindex[lno];            /* edge count */

        for (j = hg->vindex[lno]; j < hg->vindex[lno+1]; j++) { 
          *intptr++ = hg->vedge[j];                                   /* edge lno */
        }

        floatptr = (float *)intptr;

        for (j = 0, fptr = cw + (lno*VtxDim); j < VtxDim; j++)  {
          *floatptr++ = *fptr++;
        }

        s = (char *)floatptr;
       
/*      } */
    }

    /* communication to determine global size of rec buffer */
    MPI_Allgather (&sendsize, 1, MPI_INT, size, 1, MPI_INT, hgc->row_comm);
    
    /* determine size of the rec buffer & reallocate bigger iff necessary */
    recsize = 0;
    for (i = 0; i < hgc->nProc_x; i++)
      recsize += size[i];          /* compute total size of edgebuf in ints */
    if (recsize > nEdgebuf){
      MACRO_RESIZE (1.2 * recsize, nEdgebuf, edgebuf);  /* enlarge edgebuf */
    }
    
    /* setup displacement array necessary for MPI_Allgatherv */
    dest[0] = 0;
    for (i = 1; i < hgc->nProc_x; i++)
      dest[i] = dest[i-1] + size[i-1];

    /* communicate vertices & their edges to all row neighbors */
    MPI_Allgatherv(sendbuf, sendsize, MPI_CHAR, edgebuf, size, dest, MPI_CHAR, hgc->row_comm);

    /* Communication has grouped candidates by processor, rescramble!     */
    /* Otherwise all candidates from proc column 0 will be matched first, */
    for (i = 0; i < total_nCandidates; i++)
      candvisit[i] = i;
    Zoltan_Rand_Perm_Int (candvisit, total_nCandidates, &(hgc->RNGState_col));

    for (i = 0; i < total_nCandidates; i++)
      idxptr[i] = -1;                 /* to flag missing sparse candidates */

    r = edgebuf;
    n = 0;

    while (r < edgebuf + recsize){

      gnoptr = (ZOLTAN_GNO_TYPE *)r;
      intptr = (int *)(gnoptr + 1);

      candidate_index = intptr[0];
      if (hgp->UsePrefPart){
        count = intptr[2];
      }
      else{
        count = intptr[1];
      }

      floatptr = (float *)(intptr + header_size + count);

      gno_locs[n] = (char *)gnoptr - edgebuf;

      idxptr[candidate_index] = n;

      memcpy(candw + candidate_index*VtxDim, floatptr, sizeof(float)*VtxDim);

      r = (char *)(floatptr + VtxDim);
      n++;
    }

    MACRO_TIMER_STOP (1);

    /************************ PHASE 2: ***************************************/
      
    /* for each candidate vertex, compute all local partial inner products */    

    kstart = old_kstart = 0;         /* next candidate (of nTotal) to process */
    while (kstart < total_nCandidates) {
      MACRO_TIMER_START (2, "Matching kstart A", 0);
      for (i = 0; i < hgc->nProc_y; i++)
        rows[i] = -1;                  /* to flag data not found for that row */

      sendsize = 0;                    /* position in send buffer */
      sendcnt  = 0;                    /* count of messages in send buffer */

      for (k = kstart; k < total_nCandidates; k++)  {

        n = idxptr[candvisit[k]];

        if (n == -1) 
          continue;                /* don't have this sparse candidate locally */

        gnoptr = (ZOLTAN_GNO_TYPE *)(edgebuf + gno_locs[n]);

        intptr = (int *)(gnoptr + 1);

        candidate_gno   = *gnoptr;          /* gno of candidate vertex */
        candidate_index = *intptr++;          /* candidate_index of vertex */
        if (hgp->UsePrefPart)          
          pref = *intptr++;     /* pref vertex information */          
    
        count = *intptr++;      /* intptr points to edge lnos */

        /* now compute the row's nVtx inner products for kth candidate */
        m = 0;
        if (hg->ewgt != NULL) 
          AGG_INNER_PRODUCT1(hg->ewgt[*intptr])
        else 
          AGG_INNER_PRODUCT1(1.0)

           /* if local vtx, remove self inner product (useless maximum) */
        if (VTX_TO_PROC_X(hg, candidate_gno) == hgc->myProc_x)
          sums[VTX_GNO_TO_LNO(hg, candidate_gno)] = 0.0;
        
        /* if it is partitioning with preferred parts and/or fixed vertices
           check if matches are OK also eliminate sending value 0.0*/
        for (count=0; count<m; ) 
            if (sums[aux[count]]>PSUM_THRESHOLD
                && AGG_MATCH_OK(hgp, hg, pref, lheadpref[aux[count]]))
                ++count;
            else {
                sums[aux[count]] = 0.0;
                aux[count] = aux[--m];
            }
        if (count == 0)
          continue;         /* no partial sums to append to message */       

        /* iff necessary, resize send buffer to fit at least first message */

        msgsize = sizeof(ZOLTAN_GNO_TYPE) + ((2 + count) * sizeof(int)) + (count * sizeof(float));

        if (sendcnt == 0 && (msgsize > nSend)) {
          MACRO_RESIZE (1.2 * msgsize, nSend, sendbuf);  /* increase buffer size */
        }

        if (sendsize + msgsize <= nSend)  {
          /* flag first data in each destination row for merging */
          if (rows[candidate_gno % hgc->nProc_y] != 1)  {
            rows[candidate_gno % hgc->nProc_y] = 1;
            candidate_index = -candidate_index -1;
          } 

          /* current partial sums fit, so put them into the send buffer */

          gnoptr = (ZOLTAN_GNO_TYPE *)(sendbuf + sendsize);
          intptr = (int *)(gnoptr + 1);
          floatptr = (float *)(intptr + 2 + count);

          dest[sendcnt]   = candidate_gno % hgc->nProc_y;    /* destination */            
          size[sendcnt++] = msgsize;          /* size of message */
          sendsize       += msgsize;          /* cummulative size of message */
          *gnoptr++ = candidate_gno;
          *intptr++ = candidate_index;        
          *intptr++ = count;
          for (i = 0; i < count; i++)  {          
            *intptr++ = aux[i];                          /* lno of partial sum */
            *floatptr++ = sums[aux[i]];                      /* partial sum */           
            sums[aux[i]] = 0.0;
          }          
        }
        else  {           /* psum message doesn't fit into buffer */
          for (i = 0; i < count; i++)              
            sums[aux[i]] = 0.0;
/*          uprintf(hgc, "WARYNING: psum msg doesn't fit into buffer!!!! sub-phases!!!!!\n"); */
          break;   
        }  
      }                  /* DONE: loop over k */                    

      MACRO_TIMER_STOP (2);
      MACRO_TIMER_START (3, "Matching kstart B", 0);
    
      /* synchronize all rows in this column to next kstart value */
      old_kstart = kstart;
      MPI_Allreduce (&k, &kstart, 1, MPI_INT, MPI_MIN, hgc->col_comm);

      /* Send inner product data in send buffer to appropriate rows */
      ierr = communication_by_plan (zz, sendcnt, dest, size, 1, sendbuf, &reccnt, 
                                    &recsize, &nRec, &recvbuf, hgc->col_comm, IPM_TAG);
      if (ierr != ZOLTAN_OK)
        goto End;
    
      /* build index into receive buffer pointer for each proc's row of data */
      for (i = 0; i < hgc->nProc_y; i++)
        rows[i] = recsize;       /* sentinal in case no row's data available */
      k = 0;
      r = recvbuf;

      while (r < recvbuf + recsize){

        gnoptr = (ZOLTAN_GNO_TYPE *)r;
        intptr = (int *)(gnoptr + 1);
        count = intptr[1];
        floatptr = (float *)(intptr + 2 + count);
        r = (char *)(floatptr + count);

        if (*intptr < 0)  {
          *intptr = -(*intptr) - 1;           /* make sentinal candidate_index positive */
          rows[k++] = (char *)gnoptr - recvbuf;  /* points to gno */
        }
      }

      /* UVCUVC */
      if (k>hgc->nProc_y)
          errexit("k(%d)!=nProc_y(%d) recsize %d", k, hgc->nProc_y, recsize);
    
      /* merge partial i.p. sum data to compute total inner products */

      sendsize = 0;
      
      for (n = old_kstart; n < kstart; n++) {
        m = 0;       
        for (i = 0; i < hgc->nProc_y; i++) {
          if (rows[i] < recsize){
            gnoptr = (ZOLTAN_GNO_TYPE *)(recvbuf + rows[i]);
            intptr = (int *)(gnoptr + 1);

            if (intptr[0] == candvisit[n])  {

              candidate_index = intptr[0];
              candidate_gno   = gnoptr[0];
              count           = intptr[1];
              floatptr = (float *)(intptr + 2 + count);

              for (j = 0; j < count; j++)  {
                lno = intptr[2+j];
                if (sums[lno] == 0.0)   /* is this first time for this lno? */ 
                    aux[m++] = lno;     /* then save the lno */
                sums[lno] += *floatptr++;    /* sum the psums */
              }
              rows[i] = (char *)floatptr - recvbuf;
            }
          }
        }

        bestsum = -1.0;                        /* any negative value will do */
        bestlno = -1;                          /* any negative value will do */
/*
  if (m)
            uprintf(hgc, "for index=%d Candidate = %d\n", n, candidate_gno);
        else
            uprintf(hgc, "for index=%d no local partner\n", n);
*/
        for (i = 0; i < m; i++)  {
          float val, mcw=0.0;
          int  lno=aux[i];


          for (j=0; j<VtxDim; ++j)
              if (cw[lno*VtxDim+j]>mcw)
                  mcw = cw[lno*VtxDim+j];
          if (mcw==0.0)
              mcw = 1.0;
          val = sums[lno] / mcw;
/*          printf("[v=%d (lno=%zd) ip=%f mcw=%f val=%f] ", VTX_LNO_TO_GNO(hg, lno), lno, sums[lno], mcw, val); */
          if (val > bestsum && match[lno]>=0)  {
            bestsum = val;
            bestlno = lno;
          }
          sums[lno] = 0.0;  
        }

        /*
        if (m)
            printf("\n");
        */
        
        /* Put <candidate_gno, candidate_index, count, <lno, tsum>> into send */

        msgsize = sizeof(ZOLTAN_GNO_TYPE) + (2 * sizeof(int)) + sizeof(float);

        if (bestlno >= 0)  {
          if (sendsize + msgsize > nSend ) {
            uprintf(hgc, "resize with nSend=%d sendsize=%d m=%d\n", nSend, sendsize, m);
            MACRO_REALLOC (1.2*(sendsize + msgsize), nSend, sendbuf);
          }      
          gnoptr = (ZOLTAN_GNO_TYPE *)(sendbuf + sendsize);
          intptr = (int *)(gnoptr + 1);
          floatptr = (float *)(intptr + 2);

          *gnoptr = candidate_gno;
          *intptr++ = candidate_index;
          *intptr++ = bestlno;
          *floatptr = bestsum;
/*          uprintf(hgc, "cand_gno=%d partner_lno=%d with ip=%f\n", candidate_gno, bestlno, bestsum);*/

          sendsize += msgsize;
        }     
      }
    
      /* Communicate total inner product results to MASTER ROW */
      MPI_Gather(&sendsize, 1, MPI_INT, size, 1, MPI_INT, 0, hgc->col_comm);
    
      if (hgc->myProc_y == 0) {
        recsize = 0;
        for (i = 0; i < hgc->nProc_y; i++)
          recsize += size[i];        
          
        dest[0] = 0;
        for (i = 1; i < hgc->nProc_y; i++)
          dest[i] = dest[i-1] + size[i-1];
          
        if (recsize > nRec)
          MACRO_RESIZE (1.2 * recsize, nRec, recvbuf);      /* make rec buffer bigger */
      }
      
      MPI_Gatherv (sendbuf, sendsize, MPI_CHAR, recvbuf, size, dest, MPI_CHAR, 0, hgc->col_comm);

/*      uprintf(hgc, "recsize=%d\n", recsize);*/
      
      /* Determine best vertex and best sum for each candidate */
      if (hgc->myProc_y == 0) {   /* do following only if I am the MASTER ROW */
        for (r = recvbuf; r < recvbuf + recsize;)  {
          gnoptr = (ZOLTAN_GNO_TYPE *)r;
          intptr = (int *)(gnoptr + 1);
          floatptr = (float *)(intptr + 2);

          r = (char *)(floatptr + 1);

          candidate_gno   = *gnoptr;
          candidate_index = *intptr++;
          bestlno         = *intptr++;                    /* count of nonzero pairs */
          bestsum         = *floatptr;

/*          uprintf(hgc, "local best for cand %zd (idx=%d) is %d with ip=%f\n", candidate_gno, candidate_index, (bestlno<0) ? -1 : match[bestlno], bestsum);*/
          master_data[candidate_index].candidate = candidate_gno;
          master_data[candidate_index].partner = match[bestlno];          
          master_data[candidate_index].ip = bestsum;
          if (match[bestlno]<0)
              errexit("hey for bestlno: match[%d]=%zd\n", bestlno, match[bestlno]);
        }
      } 

      MACRO_TIMER_STOP (3);    
    }            /* DONE: kstart < max_nTotal loop */



    /************************ PHASE 3 & 4 ********************************/
    
    MACRO_TIMER_START (4, "Matching Phase 3&4", 1);

    replycnt=0;
    msgsize =  header_size * sizeof(int) + VtxDim * sizeof(float);

    if (hgc->myProc_y == 0) {      
      HG_Ptr = hg;

      MPI_Allreduce(master_data, global_best, total_nCandidates,
                    phasethreetype, phasethreeop, hgc->row_comm);

      /*
      uprintf(hgc, "Cand\tPartner\tIP\n");
      for (i = 0; i < total_nCandidates; i++)
        uprintf(hgc, "%d\t%d\t%.3f\n",  global_best[i].candidate,
                global_best[i].partner, global_best[i].ip);
      */

      /* Reinitialize master_data for next round */
      for (i = 0; i < total_nCandidates; i++) {        
        master_data[i].candidate = -1;
        master_data[i].partner   = -1;
        master_data[i].ip        = -1.0;
      }

      if (total_nCandidates * msgsize > nSend) 
        MACRO_RESIZE (total_nCandidates * msgsize, nSend, sendbuf);  /* increase buffer size */

      s = sendbuf;         
      
      for (i = 0; i < total_nCandidates; i++) {
        int cproc, pproc;
        candidate_gno = global_best[i].candidate;
        
        if (candidate_gno == -1)
          continue;
        
        partner_gno = global_best[i].partner;
        cproc = VTX_TO_PROC_X(hg, candidate_gno);
        pproc = VTX_TO_PROC_X(hg, partner_gno);
        if (pproc == hgc->myProc_x)   { /* partner is mine */
          /* check weight constraints */
          
          int plno = VTX_GNO_TO_LNO(hg, partner_gno);
          for (j=0; j<VtxDim; ++j)
            if (cw[plno*VtxDim+j]+candw[i*VtxDim+j]>maxw[j])
              break;
          dest[replycnt++] = cproc;

          intptr = (int *)s;
          floatptr = (float *)(intptr + header_size);
          s = (char *)(floatptr + VtxDim);

          *intptr++ = i;
          *intptr++ = plno;

          if (hgp->UsePrefPart)
              *intptr++ = lheadpref[plno];

          if (j<VtxDim) { /* reject due to weight constraint*/
            *floatptr = -1.0; /* negative means rejected */
/*            uprintf(hgc, "I'm rejecting (%zd, %zd)\n", candidate_gno, partner_gno); */

          } else { /* accept */ 
            for (j=0; j<VtxDim; ++j) { /* modify weight immediately */
              cw[plno*VtxDim+j] += candw[i*VtxDim+j];
              *floatptr++ = cw[plno*VtxDim+j];
            }
            visited[plno] = 1;
/* this printf only works if all vertices are local
   uprintf(hgc, "I'm acceptiong (%d [%d], %d[%d])\n", candidate_gno, lheadpref[candidate_gno], partner_gno, lheadpref[partner_gno]);
*/
            /* was partner a candidate ?*/
            if (match[plno]<0) 
              errexit("HEY HEY HEY  match[%d(gno=%zd)]=%zd\n", plno, partner_gno, match[plno]);
          }
        }
      }
    }

    for (i=0; i<locCandCnt; ++i) {
        int lno=locCandidates[i];
        if (match[lno]<0)
            match[lno] = -match[lno]-1;
        else
            errexit("hey hey hey match[%d]=%zd\n", lno, match[lno]);
    }

    /* bcast accepted match to column so that they can
       set visited array for local partners and also set cw properly */   
    MPI_Bcast (&replycnt, 1, MPI_INT, 0, hgc->col_comm);
    
    if (hgc->myProc_y!=0 && (replycnt*msgsize)>nSend)
      MACRO_REALLOC (replycnt*msgsize, nSend, sendbuf);  /* increase buffer size */
    
    MPI_Bcast (sendbuf, replycnt*msgsize, MPI_CHAR, 0, hgc->col_comm);
    if (hgc->myProc_y!=0) {
      int plno;
      
      s = sendbuf;
      for (i=0; i<replycnt; ++i) {

        intptr = (int *)s;
        floatptr = (float *)(intptr + header_size);
        s = (char *)(floatptr + VtxDim); 

        intptr++; /* skip candidate_index*/

        plno=*intptr++;

        if (*floatptr < 0.0) { /* reject due to weight constraint*/
        } else {
/*        uprintf(hgc, "Set visited flag of %d (gno=%d)\n", plno, VTX_LNO_TO_GNO(hg, plno)); */
          visited[plno] = 1; 

          memcpy(&cw[plno*VtxDim], floatptr, sizeof(float)*VtxDim); 
        }
      }      
    }


    if (hgc->myProc_y == 0) {    
      /* send accept/reject message */
      communication_by_plan(zz, replycnt, dest, NULL, msgsize, sendbuf,
                            &reccnt, &recsize, &nRec, &recvbuf, hgc->row_comm, CONFLICT_TAG);

      msgsize = sizeof(ZOLTAN_GNO_TYPE) + header_size * sizeof(int) + VtxDim * sizeof(float);

      if (reccnt*msgsize > nSend) 
        MACRO_RESIZE (reccnt*msgsize, nSend, sendbuf);  /* increase buffer size */
      s = sendbuf;         
      for (r = recvbuf; r < recvbuf + recsize;) {
        int ci, lheadno, pref;
        ZOLTAN_GNO_TYPE partner;

        intptr = (int *)r;
        floatptr = (float *)(intptr + header_size);
        r = (char *)(floatptr + VtxDim);

        ci = *intptr++;

        lno=VTX_GNO_TO_LNO(hg, global_best[ci].candidate);    /* global_best is not one of my columns */

        partner=global_best[ci].partner;
        pref = -1;

        ++intptr;                     /* skip plno */

        if (hgp->UsePrefPart)
            pref = *intptr++;

        if (*floatptr < 0.0) { /* rejected */
          lheadno = -1;
/*          uprintf(hgc, "(%zd, %zd) has been rejected\n", global_best[ci].candidate, global_best[ci].partner);*/
        } else { /* accepted */
          lheadno = Zoltan_KVHash_Insert(&hash, partner, lno); 

          for (j=0; j<VtxDim; ++j) 
            cw[lheadno*VtxDim+j] = floatptr[j];
          if (hgp->UsePrefPart)
              lheadpref[lno] = pref;
          lhead[lno] = lheadno;
          match[lno] = partner;
          
/*          uprintf(hgc, "(%zd, %zd) has been accepted\n", global_best[ci].candidate, global_best[ci].partner);*/
        }

        gnoptr = (ZOLTAN_GNO_TYPE *)s;
        intptr = (int *)(gnoptr + 1);
        floatptr = (float *)(intptr + header_size);
        s = (char *)(floatptr + VtxDim);
        
        *gnoptr = partner;  
        *intptr++ = lno;
        *intptr++ = lheadno;
        if (hgp->UsePrefPart)
            *intptr++ = pref;
        if (lheadno!=-1){
          memcpy(floatptr, &cw[lheadno*VtxDim], sizeof(float)*VtxDim);
        }
      }
      recsize = s - sendbuf;
    }

    MPI_Bcast (&recsize, 1, MPI_INT, 0, hgc->col_comm); /* bcase nSend */
    if (recsize>nSend) /* for procs other than 0; that might be true */
      MACRO_RESIZE (recsize, nSend, sendbuf);  /* increase buffer size */

    MPI_Bcast (sendbuf, recsize, MPI_CHAR, 0, hgc->col_comm);

    if (hgc->myProc_y !=0) { /* master row already done this */
      int lno, lheadno, partner, pref;
      for (s = sendbuf; s < sendbuf + recsize; ) {

        gnoptr = (ZOLTAN_GNO_TYPE *)s;
        intptr = (int *)(gnoptr + 1);
        floatptr = (float *)(intptr + header_size);
        s = (char *)(floatptr + VtxDim);

        lno     = *intptr++;
        lheadno = *intptr++;
        partner = *gnoptr;
        pref = 0;
        
        if (hgp->UsePrefPart)
            pref = *intptr++;

        if (lheadno!=-1) { 
          lhead[lno] = lheadno;
          match[lno] = partner;
          memcpy(&cw[lheadno*VtxDim], floatptr, sizeof(float)*VtxDim); 

          if (hgp->UsePrefPart)
              lheadpref[lheadno] = pref;
        }
      }
    }

    for (; vindex < hg->nVtx && visited[visit[vindex]]; ++vindex);
    
    i = vindex < hg->nVtx;
    MPI_Allreduce(&i, &nRounds, 1, MPI_INT, MPI_SUM, hgc->row_comm);
    MACRO_TIMER_STOP (5);                       /* end of phase 4 */
  }                                             /* DONE: loop over rounds */

 End:
  MPI_Op_free(&phasethreeop);
  MPI_Type_free(&phasethreetype);
  ZOLTAN_FREE(&global_best);
  ZOLTAN_FREE(&gno_locs);

  if (hgc->myProc_y==0 && total_nCandidates)
    Zoltan_KVHash_Destroy(&hash);

  
ZOLTAN_FREE(&candIdx);
ZOLTAN_FREE(&cw);
ZOLTAN_FREE(&tw);
ZOLTAN_FREE(&maxw);
ZOLTAN_FREE(&candw);
ZOLTAN_FREE(&lhead);
ZOLTAN_FREE(&lheadpref);
ZOLTAN_FREE(&visit);
ZOLTAN_FREE(&visited);
ZOLTAN_FREE(&sums);
ZOLTAN_FREE(&sendbuf);
ZOLTAN_FREE(&dest);
ZOLTAN_FREE(&size);
ZOLTAN_FREE(&recvbuf);
ZOLTAN_FREE(&aux);
ZOLTAN_FREE(&idxptr);
ZOLTAN_FREE(&candvisit);
ZOLTAN_FREE(&edgebuf);
ZOLTAN_FREE(&locCandidates);
ZOLTAN_FREE(&rows);
ZOLTAN_FREE(&master_data);
ZOLTAN_FREE(&master_procs);
    
/*
  Zoltan_Multifree (__FILE__, __LINE__, 22, &candIdx, &cw, &tw, &maxw, &candw, &lhead, &lheadpref,
                    &visit, &visited, &sums, &sendbuf, &dest, &size, &recvbuf, &aux, &idxptr, &candvisit,
                    &edgebuf, &locCandidates, &rows, &master_data, &master_procs);
*/
  ZOLTAN_TRACE_EXIT(zz, yo);
    
  return ierr;
}


/****************************************************************************/
static int pmatching_geom (ZZ *zz,
			  HGraph *hg,
			  ZOLTAN_GNO_TYPE *match, /* the matching array */
                          PHGPartParams *hgp)
{
  int ierr = ZOLTAN_OK;

  char *yo = "pmatching_geom";
  ZZ *zz2 = Zoltan_Create(hg->comm->Communicator);
  int i;
  char s[8]; 
  int changes, num_gid_entries, num_lid_entries, local_vtx;
  ZOLTAN_GNO_TYPE *procmatch;
  MPI_Datatype zoltan_gno_mpi_type;

  /* --LB_Partition arguments-- */
  int num_import;                         /* Returned */
  ZOLTAN_ID_PTR import_global_ids = NULL; /* Returned */
  ZOLTAN_ID_PTR import_local_ids = NULL;  /* Returned */
  int *import_procs = NULL;               /* Returned */
  int *import_to_part = NULL;             /* Returned */
  int num_export;                  /* Returned: number of input GIDs */
  ZOLTAN_ID_PTR candidate_ids = NULL; /* Returned:  candidates for each
                                             input GID */
  ZOLTAN_ID_PTR export_local_ids = NULL;  /* Not computed */
  int *export_procs = NULL;               /* Not computed */
  int *export_to_part = NULL;             /* Not computed */
  /* ----------------- */

  /* Register new geometric callbacks and parameters */
  if (Zoltan_Set_Fn(zz2, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) geometric_get_num_obj,
		    (void *) hg) == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Fn()\n");
    goto End;
  }

  if (Zoltan_Set_Fn(zz2, ZOLTAN_OBJ_LIST_FN_TYPE,(void (*)()) geometric_get_obj_list,
		    (void *) hg) == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Fn()\n");
    goto End;
  }

  if (Zoltan_Set_Fn(zz2, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) geometric_get_num_geom,
		    (void *) hg) == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Fn()\n");
    goto End;
  }

  if (Zoltan_Set_Fn(zz2, ZOLTAN_GEOM_MULTI_FN_TYPE, (void (*)()) geometric_get_geom_multi,
		    (void *) hg) == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Fn()\n");
    goto End;
  }
    
  sprintf(s, "%d", 1);
  if (Zoltan_Set_Param(zz2, "NUM_GID_ENTRIES", s) == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }

  if (Zoltan_Set_Param(zz2, "NUM_LID_ENTRIES", s) == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }

  if (Zoltan_Set_Param(zz2, "DEBUG_LEVEL", "0") == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }

  if ((hgp->geometric_red <= 0) || (hgp->geometric_red > 1)) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid hybrid reduction factor. Using default value 0.1.");
      hgp->geometric_red = 0.1;
  }
  /* Parts are reduced by a factor, result should not be 0 */  
  local_vtx = (int)(hgp->geometric_red * (double) hg->nVtx
                                       / (double) hg->comm->nProc_y);
  sprintf(s, "%d", (local_vtx == 0) ? 1 : local_vtx);
  if (Zoltan_Set_Param(zz2, "NUM_LOCAL_PARTS", s) == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }

  if (Zoltan_Set_Param(zz2, "OBJ_WEIGHT_DIM", "1") == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }

  if (Zoltan_Set_Param(zz2, "REMAP", "0") == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }

  if (Zoltan_Set_Param(zz2, "CHECK_GEOM", "0")== ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }

  if (Zoltan_Set_Param(zz2, "RETURN_LISTS", "CANDIDATE_LISTS")== ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }
	
  if (Zoltan_Set_Param(zz2, "LB_METHOD", hgp->redm_str) == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "fatal: error returned from Zoltan_Set_Param()\n");
    goto End;
  }

  /* using the wrapper (just like any other RCB/RIB user) */
  ierr = Zoltan_LB_Partition(zz2, &changes, &num_gid_entries, &num_lid_entries,
                    &num_import, &import_global_ids,
		    &import_local_ids, &import_procs, &import_to_part,
		    &num_export, &candidate_ids, &export_local_ids,
		    &export_procs, &export_to_part);

  if(!(procmatch = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC(hg->nVtx, sizeof(ZOLTAN_GNO_TYPE))))
    MEMORY_ERROR;

  for (i = 0; i < hg->nVtx; i++)
    match[i] = 0;

  for (i = 0; i < num_import; i++)
    procmatch[import_local_ids[i]] = candidate_ids[i];

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();
  MPI_Allreduce(procmatch, match, hg->nVtx, zoltan_gno_mpi_type, MPI_SUM, hg->comm->col_comm);

#ifdef KDDKDD_DEBUG
 {/* KDDKDD */
   int kdd;
   printf("%d KDDKDD Sanity Check %d == %d ? \n", zz->Proc, num_import, num_export);
   for (kdd = 0; kdd < num_import; kdd++) { 
     printf("%d KDDKDD Input (%d %d)  (%f %f)  wgt %f Candidate %d\n", zz->Proc, import_global_ids[kdd], import_local_ids[kdd], hg->coor[2*import_local_ids[kdd]], hg->coor[2*import_local_ids[kdd]+1], hg->vwgt[import_local_ids[kdd]], candidate_ids[kdd]); 
   } 
 /* KDDKDD */

 } 
#endif    

  /* 
   * Perform geometric matching only on geometric_levels levels.  If done,
   * switch to agglomerative matching.
   */
  if (hg->info+1 >= hgp->geometric_levels) {
    hgp->matching = pmatching_agg_ipm;
    sprintf(hgp->redm_str, "agg");
  }
  if (zz->Proc == 0) printf("KDDKDD GEOM_MATCHING %d %d %x\n", hg->info, hgp->geometric_levels, hgp->matching);

 End:
  Zoltan_Destroy(&zz2);
  ZOLTAN_FREE(&procmatch);
  ZOLTAN_FREE(&import_global_ids);
  ZOLTAN_FREE(&import_local_ids);
  ZOLTAN_FREE(&import_procs);
  ZOLTAN_FREE(&import_to_part);
  ZOLTAN_FREE(&candidate_ids);

  ZOLTAN_TRACE_EXIT(zz, yo);
  
  return ierr;

}

static int geometric_get_num_obj(void *data, int *ierr) {
  /* Return number of vertices / number of processes on the y-axis */
  HGraph *hg;
  int count, diff;
  
  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  hg = ((HGraph*) data);
  *ierr = ZOLTAN_OK;

  count = hg->nVtx / hg->comm->nProc_y;
  diff  = hg->nVtx % hg->comm->nProc_y;
  
/* Use processor ID to select which will take extra processes */
  return hg->comm->myProc_y < diff? count + 1 : count;
}

static void geometric_get_obj_list(void *data, int num_gid, int num_lid,
			     ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
			     int wdim, float *wgt, int *ierr) {
  /* Return GNOs, LNOs, and weights */
  
  HGraph *hg;
  int local_obj;  /* Number of objects endemic to this processor */
  int count, diff, idx;
  int i, j;
  
  if (data == NULL) { 
    *ierr = ZOLTAN_FATAL; 
    return;
  } 

  hg = ((HGraph*) data);
  *ierr = ZOLTAN_OK;

  count = hg->nVtx / hg->comm->nProc_y;
  diff  = hg->nVtx % hg->comm->nProc_y;
  local_obj =  hg->comm->myProc_y < diff? count + 1 : count;
  idx = (count * hg->comm->myProc_y) + (hg->comm->myProc_y > diff ?
					diff : hg->comm->myProc_y);
  
  for (i = 0; i < local_obj; i++) {
    global_id[i] = VTX_LNO_TO_GNO(hg, (idx + i));
    local_id[i] = idx + i;
  }

  if (wdim > 0)
    for (i = 0; i < local_obj; i++)
      for (j = 0; j < wdim; j++)
	wgt[(i * wdim) + j] = hg->vwgt[local_id[i] * wdim + j];
}

static int geometric_get_num_geom(void *data, int *ierr) {
/* Return number of values needed to express the geometry of the HG */

  HGraph *hg;

  if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return 0;
  }

  hg = ((HGraph*) data);
  *ierr = ZOLTAN_OK;

  return hg->nDim;
}

static void geometric_get_geom_multi(void *data, int num_obj, int num_gid,
			       int num_lid, ZOLTAN_ID_PTR global_id,
			       ZOLTAN_ID_PTR local_id, int num_dim,
			       double *coor, int *ierr) {
  /* Return coords by GNOs of object list */

    HGraph *hg;
    int local_obj;
    ZOLTAN_ID_TYPE glob;
    int i, j;
    
  if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return;
  }

  hg = ((HGraph*) data);
  *ierr = ZOLTAN_OK;

  local_obj = geometric_get_num_obj(hg, ierr);
				    
  for (i = 0; i < local_obj; i++) { 
    glob = local_id[i];
    for (j = 0; j < hg->nDim; j++)
      coor[i * hg->nDim + j] = hg->coor[glob * hg->nDim + j]; 
  }
}

#undef MACRO_REALLOC
#undef MACRO_TIMER_START
#undef MACRO_TIMER_STOP
#undef INNER_PRODUCT
#undef INNER_PRODUCT2
#undef ROUNDS_CONSTANT
#undef IPM_TAG
#undef HEADER_SIZE
#undef PSUM_THRESHOLD
#undef TSUM_THRESHOLD

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif


