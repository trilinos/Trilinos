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


#include "zz_const.h"
#include "hypergraph.h"
#include <limits.h>

  /* Generate serial HG.
   * Two cases:  
   * Case 1:  HE = [OldProc, NewProc] 
   *          with HEwgt = #obj in OldProc that are now assigned to NewProc.
   * Case 2:  HE = [OldPart, NewPart] 
   *          with HEwgt = #obj in OldPart that are now assigned to NewPart.
   */

#define HEINFO_ENTRIES 3

static int gather_and_remap(ZZ *, int, int *);
static int set_remap_flag(ZZ *);
static int malloc_HEinfo(ZZ *, int, int **);
static int do_match(ZZ*, HGraph *, int *, int);


/******************************************************************************/
int Zoltan_LB_Remap_From_Import_Lists(
  ZZ *zz,
  int nkeep,           /* # objs whose partition and processor are unchanged */
  int *keep_part,      /* partition assignments for kept objs:  old == new */
  int nimp,            /* # objs imported to this processor */
  int *imp_old_proc,   /* processors from which objs are imported */
  int *imp_old_part,   /* partitions from which objs are imported */
  int *imp_new_part    /* partitions to which objs are imported */
)
{
/*  Routine to remap partitions (to new processors or new partition numbers)
 *  to reduce data movement.
 *  This routine assumes the load-balancing algorithm built import lists.
 *  Objects described are those that ENDED UP on my_proc due to load balancing.
 *  For all these objects, new_proc == my_proc.
 */
char *yo = "Zoltan_LB_Remap_From_Import_Lists";
int ierr = ZOLTAN_OK;
int i, cnt, tmp;

int old_size;                 /* # of old entries to remap to. If remapping
                                 parts to processors, old_size = Num_Procs;
                                 if renumbering partitions, old_size = old 
                                 num parts. */
int fp;                       /* First partition on this processor in new
                                 decomposition. */
int np;                       /* # of partitions on this processor in new
                                 decomposition. */
int my_proc = zz->Proc;       /* This processor's rank. */
int minp, maxp;               /* Lowest and highest partition numbers on this
                                 processor in old decomposition;
                                 partition numbers are assumed to be dense,
                                 but no particular distribution is assumed. */
int HEwgt_size;               /* # of HE weights allocated. */
int *HEwgt = NULL;            /* Array of HE weights.  Initially includes
                                 zero weights; later zero-weights are removed.*/
int HEcnt;                    /* # of HEs allocated. */
int *HEinfo = NULL;           /* Array of HE info; for each HE, two pins and 
                                 one edge weight. Stored as a single vector
                                 to minimize communication calls.  */


  /* Determine type of remapping that is appropriate */
  ierr = set_remap_flag(zz);

  if (zz->LB.Remap_Flag == ZOLTAN_LB_REMAP_NONE) 
    return ierr;

  Zoltan_LB_Proc_To_Part(zz, my_proc, &np, &fp);
  if (zz->LB.Remap_Flag == ZOLTAN_LB_REMAP_PROCESSORS) {

    /* Renumber new processors to minimize changes in proc assignment. */

    HEwgt_size = zz->Num_Proc;
    HEwgt = (int *) ZOLTAN_CALLOC(HEwgt_size, sizeof(int));
    if (!HEwgt) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < nkeep; i++) 
      HEwgt[my_proc]++;
    for (i = 0; i < nimp; i++) 
      HEwgt[imp_old_proc[i]]++;

    HEcnt = 0;
    for (i = 0; i < HEwgt_size; i++)
      if (HEwgt[i] != 0) HEcnt++;
 
    ierr = malloc_HEinfo(zz, HEcnt, &HEinfo);
    if (ierr < 0) 
      goto End;

    cnt = 0;
    for (i = 0; i < HEwgt_size; i++) {
      if (HEwgt[i] != 0) {
        tmp = cnt * HEINFO_ENTRIES;
        HEinfo[tmp] = i;               /* Old processor number */
        HEinfo[tmp+1] = my_proc;       /* New processor number */
        HEinfo[tmp+2] = HEwgt[i];         /* shift non-zero weights down. */
        cnt++;
      }
    }
  }

  else {  /* ZOLTAN_LB_REMAP_PARTITIONS */

    /* Renumber new partitions to minimize changes in partition assignment */

    for (minp = INT_MAX, maxp = 0, i = 0; i < nkeep; i++) {
      if (keep_part[i] < minp) minp = keep_part[i];
      if (keep_part[i] > maxp) maxp = keep_part[i];
    }
    for (i = 0; i < nimp; i++) {
      if (imp_old_part[i] < minp) minp = imp_old_part[i];
      if (imp_old_part[i] > maxp) maxp = imp_old_part[i];
    }
    old_size = maxp - minp + 1; 

    HEwgt_size = np * old_size;
    HEwgt = (int *) ZOLTAN_CALLOC(HEwgt_size, sizeof(int));
    if (HEwgt_size && !HEwgt) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < nkeep; i++) {
      tmp = (keep_part[i]-fp) * old_size;
      HEwgt[tmp + (keep_part[i]-minp)]++;
    }

    for (i = 0; i < nimp; i++) {
      tmp = (imp_new_part[i] - fp) * old_size;
      HEwgt[tmp + (imp_old_part[i]-minp)]++;
    }

    HEcnt = 0;
    for (i = 0; i < HEwgt_size; i++)
      if (HEwgt[i] != 0) HEcnt++;
   
    ierr = malloc_HEinfo(zz, HEcnt, &HEinfo);
    if (ierr < 0) 
      goto End;

    cnt = 0;
    for (i = 0; i < HEwgt_size; i++) {
      if (HEwgt[i] != 0) {
        tmp = cnt * HEINFO_ENTRIES;
        HEinfo[tmp] = i%old_size + minp;  /* Old partition number */
        HEinfo[tmp+1] = i/old_size + fp; /* New partition number */
        HEinfo[tmp+2] = HEwgt[i];           /* shift non-zero weights down. */
        cnt++;
      }
    }
  }

  ZOLTAN_FREE(&HEwgt);

  ierr = gather_and_remap(zz, HEcnt, HEinfo);

End:
  
  if (HEinfo) ZOLTAN_FREE(&HEinfo);
  if (HEwgt) ZOLTAN_FREE(&HEwgt);

  return ierr;
}    

/******************************************************************************/
int Zoltan_LB_Remap_From_Export_Lists(
  ZZ *zz,
  int nkeep,           /* # objs whose partition and processor are unchanged */
  int *keep_part,      /* partition assignments for kept objs:  old == new */
  int nexp,            /* # objs exported from this processor */
  int *exp_old_part,   /* partitions from which objs are exported */
  int *exp_new_proc,   /* processors to which objs are exported */
  int *exp_new_part    /* partitions to which objs are exported */
) 
{
/*  Routine to remap partitions (to new processors or new partition numbers)
 *  to reduce data movement.
 *  This routine assumes the load-balancing algorithm built export lists.
 *  Objects described are those that STARTED on zz->Proc due to load balancing.
 *  For all these objects, old_proc == zz->Proc.
 */
char *yo = "Zoltan_LB_Remap_From_Export_Lists";
int ierr = ZOLTAN_OK;
int i, cnt, tmp;
int my_proc = zz->Proc;       /* This processor's rank. */

int nimp = 0;
int *imp_old_proc = NULL,     /* Temporary arrays if inversion of export to */
    *imp_old_part = NULL,     /* import lists is needed. */
    *imp_new_part = NULL;

int HEwgt_size;               /* # of HE weights allocated. */
int *HEwgt = NULL;            /* Array of HE weights.  Initially includes
                                 zero weights; later zero-weights are removed.*/
int HEcnt;                    /* # of HEs allocated. */
int *HEinfo = NULL;           /* Array of HE info; for each HE, two pins and 
                                 one edge weight. Stored as a single vector
                                 to minimize communication calls.  */

  /* Determine type of remapping that is appropriate */
  ierr = set_remap_flag(zz);

  if (zz->LB.Remap_Flag == ZOLTAN_LB_REMAP_NONE) 
    /* Not doing remapping */
    return ierr;

  if (zz->LB.Remap_Flag == ZOLTAN_LB_REMAP_PROCESSORS) {
    /* Build HEs based on processor assignment.
     * We know the old processor for all objects we are keeping and all
     * export objects -- it is my_proc!
     * We also know the new processor number for all objects initially on
     * my_proc (since we built export lists.)
     * This case is a special case of partition remapping; it is easy to 
     * build the hyperedges in this special case.
     */

    HEwgt_size = zz->Num_Proc;
    HEwgt = (int *) ZOLTAN_CALLOC(HEwgt_size, sizeof(int));
    if (!HEwgt) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < nkeep; i++) 
      HEwgt[my_proc]++;
    for (i = 0; i < nexp; i++)
      HEwgt[exp_new_proc[i]]++;

    HEcnt = 0;
    for (i = 0; i < HEwgt_size; i++)
      if (HEwgt[i] != 0) HEcnt++;
   
    ierr = malloc_HEinfo(zz, HEcnt, &HEinfo);
    if (ierr < 0) 
      goto End;

    cnt = 0;
    for (i = 0; i < HEwgt_size; i++) {
      if (HEwgt[i] != 0) {
        tmp = cnt * HEINFO_ENTRIES;
        HEinfo[tmp] = my_proc;    /* Old processor number */
        HEinfo[tmp+1] = i;        /* New processor number */
        HEinfo[tmp+2] = HEwgt[i];    /* shift non-zero weights down. */
        cnt++;
      }
    }
    ZOLTAN_FREE(&HEwgt);
  
    ierr = gather_and_remap(zz, HEcnt, HEinfo);
  }

  else {  /* ZOLTAN_LB_REMAP_PARTITIONS */
    /* Cannot renumber partitions given export lists without summing HE weights
     * across processors.  This summation is not straightforward.  Also, a 
     * potentially large number of HEs may exist 
     * (max_old_partition_number * zz->Num_Global_Parts).  Rather than build
     * this large matrix, just compute import lists from the export lists
     * and run the import-list algorithm.
     */
    ZOLTAN_COMM_OBJ *plan;
    int msg_tag = 22345;

    ierr = Zoltan_Comm_Create(&plan, nexp, exp_new_proc, zz->Communicator,
                              msg_tag, &nimp);

    if (nimp > 0) {
      imp_old_proc = (int *) ZOLTAN_MALLOC(3 * nimp * sizeof(int));
      imp_old_part = imp_old_proc + nimp;
      imp_new_part = imp_old_part + nimp;
      if (!imp_old_proc) {
        ierr = ZOLTAN_MEMERR;
        ZOLTAN_PRINT_ERROR(my_proc, yo, "Memory error.");
        goto End;
      }
    }

    ierr = Zoltan_Comm_Info(plan, NULL, NULL, NULL, NULL, NULL, NULL, 
                            NULL, NULL, NULL, NULL, NULL, imp_old_proc, NULL);

    msg_tag++;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) exp_old_part, sizeof(int),
                          (char *) imp_old_part);

    msg_tag++;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) exp_new_part, sizeof(int),
                          (char *) imp_new_part);

    Zoltan_Comm_Destroy(&plan);

    ierr = Zoltan_LB_Remap_From_Import_Lists(zz, nkeep, keep_part, nimp, 
                                      imp_old_proc, imp_old_part, imp_new_part);
  }

End:
  if (HEwgt) ZOLTAN_FREE(&HEwgt);
  if (HEinfo) ZOLTAN_FREE(&HEinfo);
  if (imp_old_proc) ZOLTAN_FREE(&imp_old_proc);

  return ierr;
}

/******************************************************************************/
static int set_remap_flag(
  ZZ *zz
)
{
char *yo = "set_remap_flag";
int ierr = ZOLTAN_OK;

/* Set remap flag based on distribution of partitions to processors. */

  if (zz->LB.PartDist == NULL) {
    /* # Partitions == # Processors, uniformly distributed; remap processors */
    zz->LB.Remap_Flag = ZOLTAN_LB_REMAP_PROCESSORS;
  }
  else if (!(zz->LB.Single_Proc_Per_Part)) {
    /* Partitions spread across >1 processor; remapping not supported. */
    zz->LB.Remap_Flag = ZOLTAN_LB_REMAP_NONE;
    ZOLTAN_PRINT_WARN(zz->Proc, yo, 
                      "Remapping is not available when computed partitions "
                      "are spread across multiple processors.");
    ierr = ZOLTAN_WARN;
  }
  else {
    /* # Partitions != # processors, or partitions not uniformly distributed */
    zz->LB.Remap_Flag = ZOLTAN_LB_REMAP_PARTITIONS;
  }

  return ierr;
}

/******************************************************************************/
static int do_match(
  ZZ *zz,
  HGraph *hg,   /* Hypergraph data structure on which to do the matching. */
  int *match,   /* Matching array -- output */
  int limit     /* max number of matches that are allowed */
)
{
/* Temporary function; will be replace by a real matching function later. */
int ierr = ZOLTAN_OK;
int i;
  
  /* Default initialization -- no change in mapping */
  for (i = 0; i < hg->nVtx; i++)
    match[i] = i;

  /* Silly hand-made matching -- reverse partition assignments. */
  for (i = limit; i < hg->nVtx; i++) {
    int part = i - limit;
    int tmp = 2 * limit - i - 1;
    if (part%2) continue;   /* Don't match odd-numbered parts */
    match[i] = tmp;
    match[tmp] = i;
  }


  /* Diagnostics -- print matching vector */
  if (zz->Proc == zz->Debug_Proc)
    for (i = 0; i < hg->nVtx; i++)
      printf("%d  MATCH %d<->%d\n", zz->Proc, i, match[i]);

  return ierr;
}

/******************************************************************************/
static int malloc_HEinfo(
  ZZ *zz,
  int HEcnt,    /* Number of HEs to allocate */
  int **HEinfo  /* Array of HE info; for each HE, two pins and 
                   one edge weight. Stored as a single vector
                   to minimize communication calls.  */
)
{
/* Routine for allocating HEs to use in remap's matching routine. */
char *yo = "malloc_HEinfo";
int ierr = ZOLTAN_OK;

  if (HEcnt) {
    *HEinfo = (int *) ZOLTAN_MALLOC(HEINFO_ENTRIES * HEcnt * sizeof(int));
    if (*HEinfo == NULL)  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
    }
  }
  else 
    *HEinfo = NULL;

  return ierr;
}

/******************************************************************************/
static int gather_and_remap(
  ZZ *zz, 
  int HEcnt,                  /* # of HEs allocated. */
  int *HEinfo                 /* Array of HE info; for each HE, two pins and 
                                 one edge weight. Stored as a single vector
                                 to minimize communication calls.  */
)
{
char *yo = "gather_and_remap";
int ierr = ZOLTAN_OK;
int i, uidx, tmp;
int num_move, gnum_move = 0;
int *each_size = NULL;        /* sizes (# HEs * HEINFO_ENTRIES) for each proc */
int *recvbuf = NULL;          /* Receive buffer for gatherv */
int *displs = NULL;           /* Displacement buffer for gatherv */
int send_size;                /* Local # HEs * HEINFO_ENTRIES */
int total_size;               /* Total # ints in gatherv */
int total_HEcnt;              /* Total (across all procs) number of HEs. */
int max0, max1;               /* Max values of pin 0 and pin 1 for each HE. */
int *match = NULL;            /* Vector describing the matching. 
                                 match[i] = j ==> match[j] = i ==> 
                                 vertices i and j are matched. */
int *used = NULL;             /* Vector indicating which partitions are used
                                 in the matching. */
int limit;                    /* Maximum number of matches that are allowed */
HGraph hg;                    /* Hypergraph for matching */

  /* Diagnostics:  put into high Debug_Level later */
  for (num_move = 0, i = 0; i < HEcnt; i++) {
    tmp = i * HEINFO_ENTRIES;
    if (HEinfo[tmp] != HEinfo[tmp+1]) num_move += HEinfo[tmp+2];
  }
  MPI_Allreduce(&num_move, &gnum_move, 1, MPI_INT, MPI_SUM, zz->Communicator);
  if (zz->Proc == zz->Debug_Proc)
    printf("%d REMAP--BEFORE: TOTAL NUM MOVE = %d\n\n", zz->Proc, gnum_move);

  /* Gather HEs from each processor into a local complete HG. */

  each_size = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
  if (!each_size) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  send_size = HEcnt * HEINFO_ENTRIES;
  MPI_Allgather(&send_size, 1, MPI_INT, each_size, 1, MPI_INT,
                zz->Communicator);

  for (total_size = 0, i = 0; i < zz->Num_Proc; i++) {
    total_size += each_size[i];
  }

  recvbuf = (int *) ZOLTAN_MALLOC((zz->Num_Proc + total_size) * sizeof(int));
  displs = recvbuf + total_size;
  if (!recvbuf) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  displs[0] = 0;
  for (i = 1; i < zz->Num_Proc; i++)
    displs[i] = displs[i-1] + each_size[i-1];

  MPI_Allgatherv(HEinfo, send_size, MPI_INT, 
                 recvbuf, each_size, displs, MPI_INT, zz->Communicator);

  total_HEcnt = total_size / HEINFO_ENTRIES;
  for (max0 = -1, max1 = -1, i = 0; i < total_HEcnt; i++) {
    tmp = i * HEINFO_ENTRIES;
    if (recvbuf[tmp] > max0) max0 = recvbuf[tmp];
    if (recvbuf[tmp+1] > max1) max1 = recvbuf[tmp+1];
  }
  /* Increment max0 and max1 so that they are the maximum number of unique
     pin values for pin0 and pin1 respectively; i.e., allow pin value == 0. */
  max0++;
  max1++;
  
  /* Sanity check */
  if (max1 != zz->LB.Num_Global_Parts) 
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unexpected value for max1.");

  /* Set up global HG */

  Zoltan_HG_HGraph_Init(&hg);
  if (total_HEcnt) {
    hg.nVtx = max0 + max1;  
    hg.nEdge = total_HEcnt;
    hg.nInput = total_HEcnt * 2;   /* two pins per HE */
    hg.EdgeWeightDim = 1;
    hg.ewgt = (float *) ZOLTAN_MALLOC(total_HEcnt * sizeof(float));
    hg.hindex = (int *) ZOLTAN_MALLOC((total_HEcnt + 1) * sizeof(int));
    hg.hvertex = (int *) ZOLTAN_MALLOC((hg.nInput) * sizeof(int));
    if (!hg.ewgt || !hg.hindex || !hg.hvertex) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < total_HEcnt; i++) {
      tmp = i * HEINFO_ENTRIES;
      hg.hindex[i] = i+i; 
      hg.hvertex[i+i] = recvbuf[tmp];
      hg.hvertex[i+i+1] = recvbuf[tmp+1]+max0;
      hg.ewgt[i] = recvbuf[tmp+2];
    }
    hg.hindex[total_HEcnt] = total_HEcnt + total_HEcnt;

    ierr = Zoltan_HG_Create_Mirror(zz, &hg);
    if (ierr < 0) goto End;

    /* Call serial matching routine. */
    /* Store remap vector in zz->LB. */
  
    if (zz->Proc == zz->Debug_Proc) Zoltan_HG_Print(zz, &hg, stdout);
  }

  /* Do matching */

  match = (int *) ZOLTAN_CALLOC(hg.nVtx + max1, sizeof(int));
  used = match + hg.nVtx;
  if (hg.nVtx && !match) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  limit = max0;   /* Max # matches allowed */
  do_match(zz, &hg, match, limit);

  /* Build remapping vector */

  zz->LB.Remap = (int *) ZOLTAN_MALLOC(zz->LB.Num_Global_Parts * sizeof(int));
  if (!(zz->LB.Remap)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }


  /* First, process all parts that were matched.  Mark matched parts as used. */
  for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
    zz->LB.Remap[i] = -1;  /* FOR SANITY CHECK ONLY, REMOVE LATER */
    tmp = match[i+max0];
    if (tmp != i+max0) {
      zz->LB.Remap[i] = tmp;
      used[tmp] = 1;
    }
  }

  /* Second, process all unmatched parts; assign them to unused partitions.  */

  for (uidx = 0, i = 0; i < zz->LB.Num_Global_Parts; i++) {
    tmp = match[i+max0];
    if (tmp == i+max0) {
      while (used[uidx]) uidx++;   /* Find next unused partition */
      zz->LB.Remap[i] = uidx;
      used[uidx] = 1;
    }
  }

  ZOLTAN_FREE(&match);


  /* Diagnostics:  put into high Debug_Level later */
  if (zz->Proc == zz->Debug_Proc) 
    for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
      printf("%d REMAP Part %d to Part %d\n", zz->Proc, i, zz->LB.Remap[i]);
    }

  for (gnum_move = 0, i = 0; i < hg.nEdge; i++) {
    tmp = i + i;
    if (hg.hvertex[tmp] != zz->LB.Remap[hg.hvertex[tmp+1]-max0]) {
      gnum_move += hg.ewgt[i];
    }
  }
  if (zz->Proc == zz->Debug_Proc)
    printf("%d REMAP--AFTER: TOTAL NUM MOVE = %d\n\n", zz->Proc, gnum_move);

End:
  ZOLTAN_FREE(&each_size);
  ZOLTAN_FREE(&recvbuf);
  Zoltan_HG_HGraph_Free(&hg);
  return ierr;
}

/******************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
