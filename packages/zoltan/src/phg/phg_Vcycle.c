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
#include "phg_distrib.h"
#include <limits.h>

  /*
#define PROCESSOR_REDUCTION
  */

    /*
#define _DEBUG
    */

/* maximum amount of memory that will be allocated by
   Zoltan_PHG_Compute_ConCut function */
#define MAXMEMORYALLOC 4*1024*1024

#define COMM_TAG 23973


typedef struct tagVCycle {
    HGraph           *hg;         /* for finer==NULL, hg and Part contains   */
    Partition         Part;       /* original hg and Part, don't delete them */  
    int              *vdest;      /* necessary to unredistribute             */
                                  /* vdest size = hg->nVtx
				     vdest[i] is the dest proc of vtx i in the
				     unredistributed communicators. */
    int              *vlno;       /* vlno size = hg->nVtx 
                                     vlno[i] is the local vertex number of vtx
				     i on proc vdest[i]. */
    int              *LevelMap;   /* necessary to uncoarsen                  */
                                  /* LevelMap size = hg->nVtx 
                                     LevelMap[i] is the vtx number of the
                                     coarse vertex containing fine vtx i 
                                     on the next level.
                                     LevelMap[i] = j >= 0 if local coarse vtx
                                     j is on the same processor as i;
                                     LevelMap[i] = -gno -1 < 0 if 
                                     coarse vtx gno is on a different proc
                                     from i. */
    int               LevelCnt;   /* 2 * count of external vertices matched to
                                     vertices owned by this proc. */
                                  /* # of negative values in LevelMap. 
                                     Number of external vertices being
                                     combined into coarse vertices on this
                                     processor */
    int               LevelSndCnt; /* number of vertices being returned by
                                      the Zoltan_comm_do_reverse(). Used to
                                      establish the receive buffer size. 
                                      Number of vertices I own that are being
                                      combined into a coarse vertex on another
                                      processor. */
    int              *LevelData;  /* buffer for external vertex information  */
                                  /* LevelData size  = LevelCnt
                                     LevelCnt/2 pairs of (external_lno, my_lno)
                                     describing matches made across procs.
                                     Proc owning my_lno will have the 
                                     coarse vtx resulting from the match
                                     and, thus, will have to send part
                                     assignment to external_lno when 
                                     uncoarsening.  */
    struct tagVCycle *finer; 
    struct Zoltan_Comm_Obj  *comm_plan;    
    struct Zoltan_Timer *timer;   /* timer for this level of the V-cycle */
    int timer_match;              /* Timers for various stages */
    int timer_coarse;             /* Declared static so we can accumulate */
    int timer_refine;             /* times over calls to Zoltan_PHG_Partition */
    int timer_project;
} VCycle; 



/****************************************************************************/
/* Routine to set function pointers corresponding to input-string options. */
int Zoltan_PHG_Set_Part_Options (ZZ *zz, PHGPartParams *hgp)
{
  int err;
  char *yo = "Zoltan_PHG_Set_Part_Options";  

  if (hgp->bal_tol < 1.0)  {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_BALANCE_TOLERANCE.");
    return ZOLTAN_FATAL;
  }

  /* Set coarsening method. */
  hgp->matching = NULL;
  if (!(Zoltan_PHG_Set_Matching_Fn (hgp)))  {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_COARSENING_METHOD.");
    return ZOLTAN_FATAL;
  }

  /* Set (serial) coarse partitioning method.  NOTE: May need parallel
   * partitioning method later if reduction to 1 proc fails              */
  hgp->CoarsePartition = Zoltan_PHG_Set_CoarsePartition_Fn(hgp, &err);
  if (err != ZOLTAN_OK)  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_COARSEPARTITION_METHOD.");
      return ZOLTAN_FATAL;
  }

  /* Set refinement method. */
  if (!(hgp->Refinement = Zoltan_PHG_Set_Refinement_Fn(hgp->refinement_str)))  {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_REFINEMENT_METHOD.");
    return ZOLTAN_FATAL;
  }
  return ZOLTAN_OK;
}

/******************************************************************************/



static int allocVCycle(VCycle *v)
{
  if (!v->hg || !v->hg->nVtx)
    return ZOLTAN_OK;
    
  if (!v->Part && !(v->Part = (int*) ZOLTAN_CALLOC (v->hg->nVtx, sizeof(int))))
    return ZOLTAN_MEMERR;
 
  if (!v->LevelMap 
   && !(v->LevelMap = (int*) ZOLTAN_CALLOC (v->hg->nVtx, sizeof(int))))
     return ZOLTAN_MEMERR;

  return ZOLTAN_OK;
}

/*****************************************************************************/



static VCycle *newVCycle(ZZ *zz, HGraph *hg, Partition part, VCycle *finer,
                         int vcycle_timing)
{
  VCycle *vcycle;
    
  if (!(vcycle = (VCycle*) ZOLTAN_MALLOC (sizeof(VCycle)))) 
    return NULL;
        
  vcycle->finer    = finer;
  vcycle->Part     = part;
  vcycle->vdest    = NULL;
  vcycle->vlno     = NULL;
  vcycle->LevelMap = NULL;
  vcycle->LevelData = NULL;
  vcycle->LevelCnt = 0;
  vcycle->LevelSndCnt = 0;
  vcycle->comm_plan = NULL;
  if (vcycle_timing)
    vcycle->timer = Zoltan_Timer_Create(zz->Timer);
  else
    vcycle->timer = NULL;
  vcycle->timer_match = -1;
  vcycle->timer_coarse = -1;
  vcycle->timer_refine = -1;
  vcycle->timer_project = -1;
  vcycle->hg       = hg ? hg : (HGraph*) ZOLTAN_MALLOC (sizeof(HGraph));
  if (!vcycle->hg)  {
    ZOLTAN_FREE (&vcycle);
    return NULL;
  }

  if (hg && (allocVCycle(vcycle) != ZOLTAN_OK))  {
    ZOLTAN_FREE (&vcycle->hg);
    ZOLTAN_FREE (&vcycle);
  }
  return vcycle;
}

/****************************************************************************/



/*  Main partitioning function for hypergraph partitioning. */
int Zoltan_PHG_Partition (
  ZZ *zz,               /* Zoltan data structure */
  HGraph *hg,           /* Input hypergraph to be partitioned */
  int p,                /* Input:  number partitions to be generated */
  float *part_sizes,    /* Input:  array of length p containing percentages
                           of work to be assigned to each partition */
  Partition parts,      /* Input:  initial partition #s; aligned with vtx 
                           arrays. 
                           Output:  computed partition #s */
  PHGPartParams *hgp)   /* Input:  parameters for hgraph partitioning. */
{

  PHGComm *hgc = hg->comm;
  VCycle  *vcycle=NULL, *del=NULL;
  int  i, err = ZOLTAN_OK, middle;
  int  origVcnt     = hg->dist_x[hgc->nProc_x];   /* for processor */
  int  origVedgecnt = hg->dist_y[hgc->nProc_y];   /* reduction test */
  int  prevVcnt     = 2*hg->dist_x[hgc->nProc_x]; /* initialized so that the */
  int  prevVedgecnt = 2*hg->dist_y[hgc->nProc_y]; /* while loop will be entered
						     before any coarsening */
  char *yo = "Zoltan_PHG_Partition";
  static int timer_match = -1,    /* Timers for various stages */
             timer_coarse = -1,   /* Declared static so we can accumulate */
             timer_refine = -1,   /* times over calls to Zoltan_PHG_Partition */
             timer_coarsepart = -1,
             timer_project = -1,
             timer_vcycle = -1;   /* times everything in Vcycle not included
                                     in above timers */
  int do_timing = (hgp->use_timers > 1);
  int vcycle_timing = (hgp->use_timers > 4);

  ZOLTAN_TRACE_ENTER(zz, yo);
    
  if (do_timing) {
    if (timer_vcycle < 0) 
      timer_vcycle = Zoltan_Timer_Init(zz->ZTime, 0, "Vcycle");
    if (timer_match < 0) 
      timer_match = Zoltan_Timer_Init(zz->ZTime, 1, "Matching");
    if (timer_coarse < 0) 
      timer_coarse = Zoltan_Timer_Init(zz->ZTime, 1, "Coarsening");
    if (timer_coarsepart < 0)
      timer_coarsepart = Zoltan_Timer_Init(zz->ZTime, 1,
                                           "Coarse_Partition");
    if (timer_refine < 0) 
      timer_refine = Zoltan_Timer_Init(zz->ZTime, 1, "Refinement");
    if (timer_project < 0) 
      timer_project = Zoltan_Timer_Init(zz->ZTime, 1, "Project_Up");

    ZOLTAN_TIMER_START(zz->ZTime, timer_vcycle, hgc->Communicator);
  }

  if (!(vcycle = newVCycle(zz, hg, parts, NULL, vcycle_timing))) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "VCycle is NULL.");
    return ZOLTAN_MEMERR;
  }

  /****** Coarsening ******/    
#define COARSEN_FRACTION_LIMIT 0.9  /* Stop if we don't make much progress */
  while ((hg->redl>0) && (hg->dist_x[hgc->nProc_x] > hg->redl)
	 && ((hg->dist_x[hgc->nProc_x] < (int) (COARSEN_FRACTION_LIMIT * prevVcnt + 0.5)) /* prevVcnt initialized to 2*hg->dist_x[hgc->nProc_x] */
	     || (hg->dist_y[hgc->nProc_y] < (int) (COARSEN_FRACTION_LIMIT * prevVedgecnt + 0.5))) /* prevVedgecnt initialized to 2*hg->dist_y[hgc->nProc_y] */
    && hg->dist_y[hgc->nProc_y] && hgp->matching) {
      int *match = NULL;
      VCycle *coarser=NULL, *redistributed=NULL;
        
      prevVcnt     = hg->dist_x[hgc->nProc_x];
      prevVedgecnt = hg->dist_y[hgc->nProc_y];

#ifdef _DEBUG      
      /* UVC: load balance stats */
      Zoltan_PHG_LoadBalStat(zz, hg);
#endif
      
      if (hgp->output_level >= PHG_DEBUG_LIST) {
          uprintf(hgc,
                  "START %3d |V|=%6d |E|=%6d #pins=%6d %d/%s/%s/%s p=%d...\n",
                  hg->info, hg->nVtx, hg->nEdge, hg->nPins, hg->redl, 
                  hgp->redm_str,
                  hgp->coarsepartition_str, hgp->refinement_str, p);
          if (hgp->output_level > PHG_DEBUG_LIST) {
              err = Zoltan_HG_Info(zz, hg);
              if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
                  goto End;
          }
      }
      if (hgp->output_level >= PHG_DEBUG_PLOT)
        Zoltan_PHG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, NULL,
         "coarsening plot");

      if (do_timing) {
        ZOLTAN_TIMER_STOP(zz->ZTime, timer_vcycle, hgc->Communicator);
        ZOLTAN_TIMER_START(zz->ZTime, timer_match, hgc->Communicator);
      }
      if (vcycle_timing) {
        if (vcycle->timer_match < 0) {
          char str[80];
          sprintf(str, "VC Matching %d", hg->info);
          vcycle->timer_match = Zoltan_Timer_Init(vcycle->timer, 0, str);
        }
        ZOLTAN_TIMER_START(vcycle->timer, vcycle->timer_match,
                           hgc->Communicator);
      }

      /* Allocate and initialize Matching Array */
      if (hg->nVtx && !(match = (int*) ZOLTAN_MALLOC (hg->nVtx*sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory: Matching array");
        return ZOLTAN_MEMERR;
      }
      for (i = 0; i < hg->nVtx; i++)
        match[i] = i;
        
      /* Calculate matching (packing or grouping) */
      err = Zoltan_PHG_Matching (zz, hg, match, hgp);
      if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
        ZOLTAN_FREE ((void**) &match);
        goto End;
      }
      if (vcycle_timing)
        ZOLTAN_TIMER_STOP(vcycle->timer, vcycle->timer_match,
                          hgc->Communicator);

      if (do_timing) {
        ZOLTAN_TIMER_STOP(zz->ZTime, timer_match, hgc->Communicator);
        ZOLTAN_TIMER_START(zz->ZTime, timer_coarse, hgc->Communicator);
      }

      if (vcycle_timing) {
        if (vcycle->timer_coarse < 0) {
          char str[80];
          sprintf(str, "VC Coarsening %d", hg->info);
          vcycle->timer_coarse = Zoltan_Timer_Init(vcycle->timer, 0, str);
        }
        ZOLTAN_TIMER_START(vcycle->timer, vcycle->timer_coarse,
                           hgc->Communicator);
      }
            
      if (!(coarser = newVCycle(zz, NULL, NULL, vcycle, vcycle_timing))) {
        ZOLTAN_FREE ((void**) &match);
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "coarser is NULL.");
        goto End;
      }

      /* Construct coarse hypergraph and LevelMap */
      err = Zoltan_PHG_Coarsening (zz, hg, match, coarser->hg, vcycle->LevelMap,
       &vcycle->LevelCnt, &vcycle->LevelSndCnt, &vcycle->LevelData, 
       &vcycle->comm_plan, hgp);
      if (err != ZOLTAN_OK && err != ZOLTAN_WARN) 
        goto End;

      if (vcycle_timing)
        ZOLTAN_TIMER_STOP(vcycle->timer, vcycle->timer_coarse,
                          hgc->Communicator);
        
      if (do_timing) {
        ZOLTAN_TIMER_STOP(zz->ZTime, timer_coarse, hgc->Communicator);
        ZOLTAN_TIMER_START(zz->ZTime, timer_vcycle, hgc->Communicator);
      }

      ZOLTAN_FREE ((void**) &match);

      if ((err=allocVCycle(coarser))!= ZOLTAN_OK)
        goto End;
      vcycle = coarser;
      hg = vcycle->hg;

      if (hgc->nProc>1 &&
	  (hg->dist_x[hgc->nProc_x] < (int) (0.5 * origVcnt + 0.5) ||
	   hg->dist_y[hgc->nProc_y] < (int) (0.5 * origVedgecnt + 0.5))) {
	/* redistribute to half the processors */
	origVcnt     = hg->dist_x[hgc->nProc_x];   /* update for processor */
	origVedgecnt = hg->dist_y[hgc->nProc_y];   /* reduction test */

	if (!(hg->vmap = (int *) ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))) {
	  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory: hg->vmap");
	  return ZOLTAN_MEMERR;
	}

	for (i = 0; i < hg->nVtx; i++)
	  hg->vmap[i] = i;

	middle = (int)((float) (hgc->nProc-1) * 0.5);

	if (hgp->nProc_x_req!=1 && hgp->nProc_y_req!=1) { /* Want 2D decomp */
	  if ((middle+1) > SMALL_PRIME && Zoltan_PHG_isPrime(middle+1))
	    --middle; /* if it was prime just use one less #procs (since
			 it should be bigger than 7 it is safe to decrement) */
	}

#ifdef PROCESSOR_REDUCTION
	if (!(hgc = (PHGComm*) ZOLTAN_MALLOC (sizeof(PHGComm)))) {
	  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory: PHGComm");
	  return ZOLTAN_MEMERR;
	}

	if (!(redistributed = newVCycle(zz,NULL,NULL,vcycle,vcycle_timing))) {
	  ZOLTAN_FREE (&hgc);
	  ZOLTAN_PRINT_ERROR (zz->Proc, yo, "redistributed is NULL.");
	  goto End;
	}

	Zoltan_PHG_Redistribute(zz, hgp, hg, 0, middle, hgc, redistributed->hg,
				&vcycle->vlno, &vcycle->vdest);

	if ((err=allocVCycle(redistributed))!= ZOLTAN_OK)
	  goto End;
	vcycle = redistributed;

	if (hgc->myProc < 0)
	  /* I'm not in the redistributed part so I should go to uncoarsening
	     refinement and wait */
	  goto Refine;

	hg = vcycle->hg;
	hg->redl = hgp->redl; /* not set with hg creation */
#endif
      }
  }

  if (hgp->output_level >= PHG_DEBUG_LIST) {
    uprintf(hgc, "START %3d |V|=%6d |E|=%6d #pins=%6d %d/%s/%s/%s p=%d...\n",
     hg->info, hg->nVtx, hg->nEdge, hg->nPins, hg->redl, 
     hgp->redm_str, hgp->coarsepartition_str, hgp->refinement_str, p);
    if (hgp->output_level > PHG_DEBUG_LIST) {
      err = Zoltan_HG_Info(zz, hg);
      if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        goto End;
    }
  }
  if (hgp->output_level >= PHG_DEBUG_PLOT)
    Zoltan_PHG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, NULL,
     "coarsening plot");

  /* free array that may have been allocated in matching */
  if (hgp->vtx_scal) ZOLTAN_FREE(&(hgp->vtx_scal));

  if (do_timing) {
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_vcycle, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer_coarsepart, hgc->Communicator);
  }

  /****** Coarse Partitioning ******/
  err = Zoltan_PHG_CoarsePartition (zz, hg, p, part_sizes, vcycle->Part, hgp);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
    goto End;

  if (do_timing) {
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_coarsepart, hgc->Communicator);
    ZOLTAN_TIMER_START(zz->ZTime, timer_vcycle, hgc->Communicator);
  }

Refine:
  del = vcycle;
  /****** Uncoarsening/Refinement ******/
  while (vcycle) {
    VCycle *finer = vcycle->finer;
    hg = vcycle->hg;

    if (hgc->myProc >= 0) {
      if (do_timing) {
	ZOLTAN_TIMER_STOP(zz->ZTime, timer_vcycle, hgc->Communicator);
	ZOLTAN_TIMER_START(zz->ZTime, timer_refine, hgc->Communicator);
      }
      if (vcycle_timing) {
	if (vcycle->timer_refine < 0) {
	  char str[80];
	  sprintf(str, "VC Refinement %d", hg->info);
	  vcycle->timer_refine = Zoltan_Timer_Init(vcycle->timer, 0, str);
	}
	ZOLTAN_TIMER_START(vcycle->timer, vcycle->timer_refine,
			   hgc->Communicator);
      }

      err = Zoltan_PHG_Refinement (zz, hg, p, part_sizes, vcycle->Part, hgp);
        
      if (do_timing) {
	ZOLTAN_TIMER_STOP(zz->ZTime, timer_refine, hgc->Communicator);
	ZOLTAN_TIMER_START(zz->ZTime, timer_vcycle, hgc->Communicator);
      }
      if (vcycle_timing)
	ZOLTAN_TIMER_STOP(vcycle->timer, vcycle->timer_refine,
			  hgc->Communicator);

                          
      if (hgp->output_level >= PHG_DEBUG_LIST)     
	uprintf(hgc, 
		"FINAL %3d |V|=%6d |E|=%6d #pins=%6d %d/%s/%s/%s p=%d bal=%.2f cutl=%.2f\n",
		hg->info, hg->nVtx, hg->nEdge, hg->nPins, hg->redl, 
		hgp->redm_str,
		hgp->coarsepartition_str, hgp->refinement_str, p,
		Zoltan_PHG_Compute_Balance(zz, hg, part_sizes, p, vcycle->Part),
		Zoltan_PHG_Compute_ConCut(hgc, hg, vcycle->Part, p, &err));

      if (hgp->output_level >= PHG_DEBUG_PLOT)
	Zoltan_PHG_Plot(zz->Proc, hg->nVtx, p, hg->vindex, hg->vedge, vcycle->Part,
			"partitioned plot");
        
      if (do_timing) {
	ZOLTAN_TIMER_STOP(zz->ZTime, timer_vcycle, hgc->Communicator);
	ZOLTAN_TIMER_START(zz->ZTime, timer_project, hgc->Communicator);
      }
      if (vcycle_timing) {
	if (vcycle->timer_project < 0) {
	  char str[80];
	  sprintf(str, "VC Project Up %d", hg->info);
	  vcycle->timer_project = Zoltan_Timer_Init(vcycle->timer, 0, str);
	}
	ZOLTAN_TIMER_START(vcycle->timer, vcycle->timer_project,
			   hgc->Communicator);
      }
    }

    if (finer)  {
      int *rbuffer;
            
      /* Project coarse partition to fine partition */
      if (finer->comm_plan) {
	/* easy to assign partitions to internal matches */
	for (i = 0; i < finer->hg->nVtx; i++)
	  if (finer->LevelMap[i] >= 0)   /* if considers only the local vertices */
	    finer->Part[i] = vcycle->Part[finer->LevelMap[i]];
          
	/* now that the course partition assignments have been propagated */
	/* upward to the finer level for the local vertices, we need to  */    
	/* fill the LevelData (matched pairs of a local vertex with a    */
	/* off processor vertex) with the partition assignment of the    */
	/* local vertex - can be done totally in the finer level!        */    
	for (i = 0; i < finer->LevelCnt; i++)  {
	  ++i;          /* skip over off processor lno */
	  finer->LevelData[i] = finer->Part[finer->LevelData[i]]; 
	}
            
	/* allocate rec buffer to exchange LevelData information */
	rbuffer = NULL;
	if (finer->LevelSndCnt > 0)  {
	  rbuffer = (int*) ZOLTAN_MALLOC (2 * finer->LevelSndCnt * sizeof(int));
	  if (!rbuffer)    {
	    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
	    return ZOLTAN_MEMERR;
	  }
	}       
      
	/* get partition assignments from owners of externally matched vtxs */  
	Zoltan_Comm_Resize (finer->comm_plan, NULL, COMM_TAG, &i);
	Zoltan_Comm_Do_Reverse (finer->comm_plan, COMM_TAG+1, 
         (char*) finer->LevelData, 2 * sizeof(int), NULL, (char*) rbuffer);

	/* process data to assign partitions to expernal matches */
	for (i = 0; i < 2 * finer->LevelSndCnt;)  {
	  int lno, partition;
	  lno       = rbuffer[i++];
	  partition = rbuffer[i++];      
	  finer->Part[lno] = partition;         
	}

	ZOLTAN_FREE (&rbuffer);                  
	Zoltan_Comm_Destroy (&finer->comm_plan);                   
      } else {
	/* ints local and partition numbers */
	int *sendbuf = NULL, size = finer->vlno ? 2 * hg->nVtx : 0;
	if (finer->vlno) {
	  sendbuf = (int*) ZOLTAN_MALLOC (2 * hg->nVtx * sizeof(int));
	  if (!sendbuf) {
	    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
	    return ZOLTAN_MEMERR;
	  }

	  for (i = 0; i < hg->nVtx; ++i) {
	    sendbuf[2 * i] = finer->vlno[i];     /* assign local numbers */
	    sendbuf[2 * i + 1] = vcycle->Part[i];/* assign partition numbers */
	  }
	}

	hgc = finer->hg->comm; /* updating hgc is required when the processors
				   change */
	/* Create comm plan to unredistributed processors */
	Zoltan_Comm_Create(&finer->comm_plan, finer->vlno ? hg->nVtx : 0,
			   finer->vdest, hgc->Communicator, COMM_TAG+2, &size);
	/* allocate rec buffer to exchange sendbuf information */
	rbuffer = NULL;
	rbuffer = (int*) ZOLTAN_MALLOC (2 * finer->hg->nVtx * sizeof(int));

	if (!rbuffer) {
	  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
	  return ZOLTAN_MEMERR;
	}

	/* Use plan to send partitions to the unredistributed processors */

	Zoltan_Comm_Do(finer->comm_plan, COMM_TAG+3, (char *) sendbuf,
		       2*sizeof(int), (char *) rbuffer);

	MPI_Bcast(rbuffer, 2*finer->hg->nVtx, MPI_INT, 0, hgc->col_comm);
	
	/* process data to assign partitions to unredistributed processors */
	for (i = 0; i < 2 * finer->hg->nVtx;) {
	  int lno, partition;
	  lno       = rbuffer[i++];
	  partition = rbuffer[i++];
	  finer->Part[lno] = partition;
	}

	if (finer->vlno)
	  ZOLTAN_FREE (&sendbuf);

	ZOLTAN_FREE (&rbuffer);
	Zoltan_Comm_Destroy (&finer->comm_plan);
      }
    }

    if (do_timing) {
      ZOLTAN_TIMER_STOP(zz->ZTime, timer_project, hgc->Communicator);
      ZOLTAN_TIMER_START(zz->ZTime, timer_vcycle, hgc->Communicator);
    }
    if (vcycle_timing)
      ZOLTAN_TIMER_STOP(vcycle->timer, vcycle->timer_project,
                        hgc->Communicator);

    vcycle = finer;
  }       /* while (vcycle) */
    
End:
  vcycle = del;
  while (vcycle) {
    if (vcycle_timing) {
      Zoltan_Timer_PrintAll(vcycle->timer, 0, hgc->Communicator, stdout);
      Zoltan_Timer_Destroy(&vcycle->timer);
    }
    if (vcycle->finer) {   /* cleanup by level */
      Zoltan_HG_HGraph_Free (vcycle->hg);

      if (vcycle->LevelData)
	Zoltan_Multifree (__FILE__, __LINE__, 4, &vcycle->Part,
			  &vcycle->LevelMap, &vcycle->LevelData, &vcycle->hg);
      else if (vcycle->vlno)
	Zoltan_Multifree (__FILE__, __LINE__, 5, &vcycle->Part, &vcycle->vdest,
			  &vcycle->vlno, &vcycle->LevelMap, &vcycle->hg);
      else
	Zoltan_Multifree (__FILE__, __LINE__, 3, &vcycle->Part,
			  &vcycle->LevelMap, &vcycle->hg);
    }
    else                   /* cleanup top level */
      Zoltan_Multifree (__FILE__, __LINE__, 2, &vcycle->LevelMap,
                        &vcycle->LevelData);
    del = vcycle;
    vcycle = vcycle->finer;
    ZOLTAN_FREE(&del);
  }

  if (do_timing)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_vcycle, hgc->Communicator);
  ZOLTAN_TRACE_EXIT(zz, yo) ;
  return err;
}
    
/****************************************************************************/

double Zoltan_PHG_Compute_NetCut(
  PHGComm *hgc,
  HGraph *hg,
  Partition part,
  int p
)
{
/* Calculates the cutsize of a partition by summing the weight of all edges
 * which span more than one part. Time O(|H|). 
 * Results are returned on all processors of hgc->Communicator. 
 */
  int i, j, *netpart = NULL, *allparts = NULL;    
  double cut = 0.0, totalcut=0.0;
  char *yo = "Zoltan_PHG_Compute_NetCut";

  if (hg->nEdge && !(netpart = (int*) ZOLTAN_CALLOC (hg->nEdge, sizeof(int)))) {
    ZOLTAN_PRINT_ERROR (hgc->myProc, yo, "Memory error.");
    return ZOLTAN_MEMERR;
  }

  if (!hgc->myProc_x)
    if (hg->nEdge && 
       !(allparts = (int*) ZOLTAN_CALLOC(hgc->nProc_x*hg->nEdge,sizeof(int)))) {
      ZOLTAN_PRINT_ERROR (hgc->myProc, yo, "Memory error.");
      return ZOLTAN_MEMERR;
    }

  for (i = 0; i < hg->nEdge; ++i) 
    if (hg->hindex[i] >= hg->hindex[i+1])
       netpart[i] = -1;
    else  {
       j = hg->hindex[i];
       netpart[i] = part[hg->hvertex[j]];
       for (++j; j < hg->hindex[i+1]  &&  part[hg->hvertex[j]] == netpart[i]; ++j);
         if (j != hg->hindex[i+1])
           netpart[i] = -2;
    }

  if (hg->nEdge)
    MPI_Gather(netpart, hg->nEdge, MPI_INT, allparts, hg->nEdge, MPI_INT, 0, 
               hgc->row_comm);
  ZOLTAN_FREE (&netpart);

  if (!hgc->myProc_x) { 
    for (i = 0; i < hg->nEdge; ++i) {
      int p=-1;
      for (j = 0; j < hgc->nProc_x; ++j)
        if (allparts[j*hg->nEdge+i] == -2)
          break;
        else if (allparts[j*hg->nEdge+i] >= 0) {
          if (p == -1)
            p = allparts[j*hg->nEdge+i];
          else if (p != allparts[j*hg->nEdge+i])
            break;
        }            
      if (j < hgc->nProc_x)
        cut += (hg->ewgt ? hg->ewgt[i] : 1.0);
    }
        
    ZOLTAN_FREE (&allparts);
    MPI_Reduce (&cut, &totalcut, 1, MPI_DOUBLE, MPI_SUM, 0, hgc->col_comm);
  }

  MPI_Bcast (&totalcut, 1, MPI_DOUBLE, 0, hgc->Communicator);
  return totalcut;    
}

/****************************************************************************/



/******************************************************************************/
double Zoltan_PHG_Compute_ConCut(
  PHGComm *hgc,
  HGraph *hg,
  Partition part, 
  int p,
  int *ierr)
{
/* Calculates the cutsize of a partition. For each edge it calculates the number
 * of parts it spans across. This value minus one is the cutsize of this edge
 * and the total cutsize is the sum of the single cutsizes. Time O(|H|). 
 * NOTE:  non-zero result is returned ONLY on processor (0,0) of the 2D
 * decomposition (i.e., processor zero overall).
 */
    double cut = 0.0, totalcut=0.0;
    char *yo = "Zoltan_PHG_Compute_ConCut";
    
    if (hg->nEdge) {
        int i, j, *cuts=NULL, *rescuts=NULL, *parts, nEdge, start;
            
        nEdge = MIN(MAXMEMORYALLOC / (2*sizeof(int)*p), hg->nEdge);

        if (!(cuts = (int*) ZOLTAN_MALLOC (p * nEdge * sizeof(int)))) {
            ZOLTAN_PRINT_ERROR(hgc->myProc, yo, "Memory error.");
            *ierr = ZOLTAN_MEMERR;
            goto End;
        }   
        if (!hgc->myProc_x)
            if (!(rescuts = (int*) ZOLTAN_MALLOC (p * nEdge * sizeof(int)))) {
                ZOLTAN_PRINT_ERROR(hgc->myProc, yo, "Memory error.");
                *ierr = ZOLTAN_MEMERR;
                goto End;
            }

        start = 0;
        do {
            int nparts, end=MIN(hg->nEdge, start+nEdge);

            memset(cuts, 0, sizeof(int)*p*(end-start));
            for (i = start; i < end; ++i) {
                parts = &cuts[(i-start)*p];
                for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j) 
                    ++parts[part[hg->hvertex[j]]];
            }

            MPI_Reduce(cuts, rescuts, p*(end-start), MPI_INT, MPI_SUM, 0, 
                       hgc->row_comm);
            
            if (!hgc->myProc_x) {
                for (i = start; i < end; ++i) {
                    parts = &rescuts[(i-start)*p];
                    for (j = nparts = 0; j< p; ++j)
                        if (parts[j])
                            ++nparts;
                    if (nparts>1)
                        cut +=  ((nparts-1) * (hg->ewgt ? hg->ewgt[i] : 1.0));
                    else if (nparts==0) {
                        char msg[80];
                        sprintf(msg, "Hyperedge %d has no vertices!\n", i);
                        ZOLTAN_PRINT_ERROR(hgc->myProc, yo, msg);
                        *ierr = ZOLTAN_FATAL;
                        goto End;
                    }
                }        
            }
            start += nEdge;            
        } while (start<hg->nEdge);
        
End:

        ZOLTAN_FREE (&cuts);
        ZOLTAN_FREE (&rescuts);
    }
    
    if (!hgc->myProc_x) 
        MPI_Reduce (&cut, &totalcut, 1, MPI_DOUBLE, MPI_SUM, 0, hgc->col_comm);

    MPI_Bcast (&totalcut, 1, MPI_DOUBLE, 0, hgc->Communicator);
    return totalcut;
}

/****************************************************************************/


double Zoltan_PHG_Compute_Balance (
  ZZ *zz,
  HGraph *hg,
  float *part_sizes,
  int p,
  Partition part
)
{
  int i;
  double *lsize_w, *size_w, max_imbal, tot_w;
  char *yo = "Zoltan_PHG_Compute_Balance";
  
  if (!hg || !hg->comm || !hg->comm->row_comm)  {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to compute balance");
    return 1.0;
  }  
  
  if (!(lsize_w = (double*) ZOLTAN_CALLOC (2*p, sizeof(double)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
  }
  size_w = lsize_w + p;
  
  if (hg->vwgt)
    for (i = 0; i < hg->nVtx; i++)
      lsize_w[part[i]] += hg->vwgt[i];
  else
    for (i = 0; i < hg->nVtx; i++)
      lsize_w[part[i]]++;
        
  MPI_Allreduce(lsize_w, size_w, p, MPI_DOUBLE, MPI_SUM, hg->comm->row_comm);
  
  max_imbal = tot_w = 0.0;
  for (i = 0; i < p; i++) 
      tot_w += size_w[i];
  if (tot_w) {
      for (i = 0; i < p; i++)
          if (part_sizes[i]) {
              double ib= (size_w[i]-part_sizes[i]*tot_w)/(part_sizes[i]*tot_w);
              if (ib>max_imbal)
                  max_imbal = ib;
          }
  }

  ZOLTAN_FREE (&lsize_w);

  return  1.0+max_imbal;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
