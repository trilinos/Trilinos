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
#include "hg_hypergraph.h"
#include "hg_util.h"
#include <limits.h>

/*
 * Values indicating how partition remapping should be done.
 */
#define ZOLTAN_LB_REMAP_NONE 0
#define ZOLTAN_LB_REMAP_PROCESSORS 1
#define ZOLTAN_LB_REMAP_PARTITIONS 2

#define HEINFO_ENTRIES 3

static int gather_and_build_remap(ZZ *, int *, int, int *);
static int set_remap_type(ZZ *, int *);
static int malloc_HEinfo(ZZ *, int, int **);
static int do_match(ZZ*, HGraph *, int *, int);
static int matching_pgm(ZZ *, HGraph *, int *, int *);
static int local_HEs_from_import_lists(ZZ *, int, int, int *, int *, int *,
  int *, int **);
static int local_HEs_from_export_lists(ZZ *, int, int, int *, int *, int *,
  int *, int **);
static float measure_stays(ZZ *, HGraph *, int, char *);

/******************************************************************************/

int Zoltan_LB_Remap(
  ZZ *zz,
  int *new_map,        /* Upon return, flag indicating whether part or proc
                          assignments actually changed due to remapping. */
  int nobj,            /* # objs the processor knows about after partitioning */
  int *proc,           /* processors for the objs; 
                          if export_list_flag == 1, 
                             proc contains new proc assignment
                          else
                             proc contains old proc assignment
                          Upon return, proc contains remapped new proc 
                          assignment regardless of export_list_flag's value. */
  int *old_part,       /* old partition assignments for the objs */
  int *new_part,       /* new partition assignments for the objs.
                          Upon return, new_part contains remapped new
                          partition assignments */
  int export_list_flag /* Flag indicating whether the algorithm computes
                          export lists or import lists. The HG for matching
                          is built differently depending on whether 
                          the algorithm knows export or import info.  */
)
{
char *yo = "Zoltan_LB_Remap";
int ierr = ZOLTAN_OK;
int i;
int remap_type;               /* Type of remapping to be done: 
                                 Procs, Parts, or None */
int HEcnt = 0;                /* Number of local hyperedges */
int *HEinfo = NULL;           /* Array of HE info; for each HE, two pins and 
                                 one edge weight. Stored as a single vector
                                 to minimize communication calls.  */

  *new_map = 0;

  /* Determine type of remapping that is appropriate */
  ierr = set_remap_type(zz, &remap_type);

  if (remap_type != ZOLTAN_LB_REMAP_NONE) {
    /* Build local hyperedges */
    if (export_list_flag) 
      ierr = local_HEs_from_export_lists(zz, remap_type,
                                         nobj, proc, old_part, new_part,
                                         &HEcnt, &HEinfo);
    else 
      ierr = local_HEs_from_import_lists(zz, remap_type,
                                         nobj, proc, old_part, new_part,
                                         &HEcnt, &HEinfo);

    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building local HEs");
      goto End;
    }

    /* Gather local hyperedges to each processor; build remap vector */
    ierr = gather_and_build_remap(zz, new_map, HEcnt, HEinfo);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from gather_and_build_remap.");
      goto End;
    }
  
    if (*new_map) {
      /* Update partition and processor information for algorithms */
      for (i = 0; i < nobj; i++) {
        new_part[i] = zz->LB.Remap[new_part[i]];
        proc[i] = Zoltan_LB_Part_To_Proc(zz, new_part[i], NULL);
      }
    }
  }

End:

  ZOLTAN_FREE(&HEinfo);
  return(ierr);
}


/******************************************************************************/
static int local_HEs_from_import_lists(
  ZZ *zz,
  int remap_type,      /* type of remapping to do:  parts, procs, or none. */
  int nobj,            /* # objs the processor knows about (keep + imports) */
  int *proc,           /* On input, old processor assignment for each obj; 
                          Upon return, remapped new proc assignment for
                          each obj. */
  int *old_part,       /* old partition assignments for each objs */
  int *new_part,       /* On input, new partition assignments for each objs.
                          Upon return, remapped new partition assignments */
  int *HEcnt,          /* # of HEs allocated. */
  int **HEinfo         /* Array of HE info; for each HE, two pins and 
                          one edge weight. Stored as a single vector
                          to minimize communication calls.  */
)
{
/*  Routine to remap partitions (to new processors or new partition numbers)
 *  to reduce data movement.
 *  This routine assumes the load-balancing algorithm built import lists.
 *  Objects described are those that ENDED UP on my_proc due to load balancing.
 *  For all these objects, new_proc == my_proc.
 */
char *yo = "local_HEs_from_import_lists";
int ierr = ZOLTAN_OK;
int i, cnt, tmp;
int *tmp_HEinfo;
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


  if (remap_type == ZOLTAN_LB_REMAP_PROCESSORS) {

    /* Renumber new processors to minimize changes in proc assignment. */

    HEwgt_size = zz->Num_Proc;
    HEwgt = (int *) ZOLTAN_CALLOC(HEwgt_size, sizeof(int));
    if (!HEwgt) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < nobj; i++) 
      HEwgt[proc[i]]++;    /* At this point, proc has old proc assignments */

    *HEcnt = 0;
    for (i = 0; i < HEwgt_size; i++)
      if (HEwgt[i] != 0) (*HEcnt)++;
 
    ierr = malloc_HEinfo(zz, *HEcnt, HEinfo);
    if (ierr < 0) 
      goto End;
    tmp_HEinfo = *HEinfo;

    cnt = 0;
    for (i = 0; i < HEwgt_size; i++) {
      if (HEwgt[i] != 0) {
        tmp = cnt * HEINFO_ENTRIES;
        tmp_HEinfo[tmp] = i;               /* Old processor number */
        tmp_HEinfo[tmp+1] = my_proc;       /* New processor number */
        tmp_HEinfo[tmp+2] = HEwgt[i];      /* shift non-zero weights down. */
        cnt++;
      }
    }
  }

  else {  /* ZOLTAN_LB_REMAP_PARTITIONS */

    /* Renumber new partitions to minimize changes in partition assignment */

    for (minp = INT_MAX, maxp = 0, i = 0; i < nobj; i++) {
      if (old_part[i] < minp) minp = old_part[i];
      if (old_part[i] > maxp) maxp = old_part[i];
    }

    /* Don't include old partition numbers that are greater than 
     * zz->LB.Num_Global_Parts - 1; they are not valid values for 
     * remapping of new partition numbers. 
     */
    if (minp >= zz->LB.Num_Global_Parts)
      minp = zz->LB.Num_Global_Parts-1;
    if (maxp >= zz->LB.Num_Global_Parts)
      maxp = zz->LB.Num_Global_Parts-1;

    old_size = maxp - minp + 1; 

    Zoltan_LB_Proc_To_Part(zz, my_proc, &np, &fp);
    HEwgt_size = np * old_size;
    HEwgt = (int *) ZOLTAN_CALLOC(HEwgt_size, sizeof(int));
    if (HEwgt_size && !HEwgt) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i = 0; i < nobj; i++) {
      if (old_part[i] < zz->LB.Num_Global_Parts) {  
        /* Include only HEs to old partitions numbered 
         * 0 to zz->LB.Num_Global_Parts-1; these are the only valid
         * remapping values for the new partition numbers.
         */
        tmp = (new_part[i]-fp) * old_size;
        HEwgt[tmp + (old_part[i]-minp)]++;
      }
    }

    *HEcnt = 0;
    for (i = 0; i < HEwgt_size; i++)
      if (HEwgt[i] != 0) (*HEcnt)++;
   
    ierr = malloc_HEinfo(zz, *HEcnt, HEinfo);
    if (ierr < 0) 
      goto End;
    tmp_HEinfo = *HEinfo;

    cnt = 0;
    for (i = 0; i < HEwgt_size; i++) {
      if (HEwgt[i] != 0) {
        tmp = cnt * HEINFO_ENTRIES;
        tmp_HEinfo[tmp] = i%old_size + minp;  /* Old partition number */
        tmp_HEinfo[tmp+1] = i/old_size + fp;  /* New partition number */
        tmp_HEinfo[tmp+2] = HEwgt[i];         /* shift non-zero weights down. */
        cnt++;
      }
    }
  }

End:
  
  if (HEwgt) ZOLTAN_FREE(&HEwgt);

  return ierr;
}    

/******************************************************************************/
static int local_HEs_from_export_lists(
  ZZ *zz,
  int remap_type,      /* type of remapping to do:  parts, procs, or none. */
  int nobj,            /* # objs the processor knows about (keep + exports) */
  int *new_proc,       /* On input, new processor assignment for each obj; 
                          Upon return, remapped new proc assignment for
                          each obj. */
  int *old_part,       /* old partition assignments for each objs */
  int *new_part,       /* On input, new partition assignments for each objs.
                          Upon return, remapped new partition assignments */
  int *HEcnt,          /* # of HEs allocated. */
  int **HEinfo         /* Array of HE info; for each HE, two pins and 
                          one edge weight. Stored as a single vector
                          to minimize communication calls.  */
) 
{
/*  Routine to remap partitions (to new processors or new partition numbers)
 *  to reduce data movement.
 *  This routine assumes the load-balancing algorithm built export lists.
 *  Objects described are those that STARTED on zz->Proc due to load balancing.
 *  For all these objects, old_proc == zz->Proc.
 */
char *yo = "local_HEs_from_export_lists";
int ierr = ZOLTAN_OK;
int i, cnt, tmp;
int *tmp_HEinfo;
int my_proc = zz->Proc;       /* This processor's rank. */

int nimp = 0;
int *imp_proc = NULL,         /* Temporary arrays if inversion of export to */
    *imp_old_part = NULL,     /* import lists is needed. */
    *imp_new_part = NULL;

int HEwgt_size;               /* # of HE weights allocated. */
int *HEwgt = NULL;            /* Array of HE weights.  Initially includes
                                 zero weights; later zero-weights are removed.*/

  if (remap_type == ZOLTAN_LB_REMAP_PROCESSORS) {
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

    for (i = 0; i < nobj; i++) 
      HEwgt[new_proc[i]]++;

    *HEcnt = 0;
    for (i = 0; i < HEwgt_size; i++)
      if (HEwgt[i] != 0) (*HEcnt)++;
   
    ierr = malloc_HEinfo(zz, *HEcnt, HEinfo);
    if (ierr < 0) 
      goto End;
    tmp_HEinfo = *HEinfo;

    cnt = 0;
    for (i = 0; i < HEwgt_size; i++) {
      if (HEwgt[i] != 0) {
        tmp = cnt * HEINFO_ENTRIES;
        tmp_HEinfo[tmp] = my_proc;    /* Old processor number */
        tmp_HEinfo[tmp+1] = i;        /* New processor number */
        tmp_HEinfo[tmp+2] = HEwgt[i]; /* shift non-zero weights down. */
        cnt++;
      }
    }
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

    ierr = Zoltan_Comm_Create(&plan, nobj, new_proc, zz->Communicator,
                              msg_tag, &nimp);

    if (nimp > 0) {
      imp_proc = (int *) ZOLTAN_MALLOC(3 * nimp * sizeof(int));
      imp_old_part = imp_proc + nimp;
      imp_new_part = imp_old_part + nimp;
      if (!imp_proc) {
        ierr = ZOLTAN_MEMERR;
        ZOLTAN_PRINT_ERROR(my_proc, yo, "Memory error.");
        goto End;
      }
    }

    ierr = Zoltan_Comm_Info(plan, NULL, NULL, NULL, NULL, NULL, NULL, 
                            NULL, NULL, NULL, NULL, NULL, imp_proc, NULL);

    msg_tag++;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) old_part, sizeof(int),
                          (char *) imp_old_part);

    msg_tag++;
    ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) new_part, sizeof(int),
                          (char *) imp_new_part);

    Zoltan_Comm_Destroy(&plan);

    ierr = local_HEs_from_import_lists(zz, remap_type, nimp, imp_proc,
                                       imp_old_part, imp_new_part,
                                       HEcnt, HEinfo);
  }

End:

  if (HEwgt) ZOLTAN_FREE(&HEwgt);
  if (imp_proc) ZOLTAN_FREE(&imp_proc);

  return ierr;
}

/******************************************************************************/
static int set_remap_type(
  ZZ *zz, 
  int *remap_type
)
{
int ierr = ZOLTAN_OK;

/* Set remap type based on distribution of partitions to processors. */

  if (zz->LB.Remap_Flag == 0) {
    /* No remapping requested */
    *remap_type = ZOLTAN_LB_REMAP_NONE;
  }
  else if (!(zz->LB.Uniform_Parts)) {
    /* Remapping does not respect requested non-uniform partition sizes;
       no remapping done. */
    *remap_type = ZOLTAN_LB_REMAP_NONE;
    ierr = ZOLTAN_WARN;
  }
  else if (!(zz->LB.Single_Proc_Per_Part)) {
    /* Some partitions spread across >1 processor; remapping not supported. */
    *remap_type = ZOLTAN_LB_REMAP_NONE;
    ierr = ZOLTAN_WARN;
  }
  else if (zz->LB.PartDist == NULL) {
    /* # Partitions == # Processors, uniformly distributed; remap processors */
    *remap_type = ZOLTAN_LB_REMAP_PROCESSORS;
  }
  else {
    /* # Partitions != # processors, or partitions not uniformly distributed */
    *remap_type = ZOLTAN_LB_REMAP_PARTITIONS;
  }

  return ierr;
}

/******************************************************************************/
static int do_match(
  ZZ *zz,
  HGraph *hg,   /* Hypergraph data structure on which to do the matching. */
  int *match,   /* Matching array -- output */
  int limit    /* max number of matches that are allowed */
)
{
/* Temporary function; will be replace by a real matching function later. */
int ierr = ZOLTAN_OK;
int i;
  
  /* Default initialization -- no change in mapping */
  for (i = 0; i < hg->nVtx; i++)
    match[i] = i;

  ierr = matching_pgm(zz, hg, match, &limit);

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
static int gather_and_build_remap(
  ZZ *zz, 
  int *new_map,               /* Upon return, flag indicating whether parts
                                 assignments were changed due to remap. */
  int HEcnt,                  /* # of HEs allocated. */
  int *HEinfo                 /* Array of HE info; for each HE, two pins and 
                                 one edge weight. Stored as a single vector
                                 to minimize communication calls.  */
)
{
char *yo = "gather_and_remap";
int ierr = ZOLTAN_OK;
int i, uidx, tmp;
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
float before,                 /* Amount of data that overlaps between old and */
      after;                  /* new decomposition before and after remapping, 
                                 respectively. */


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
  /* Ideally, max1 should equal LB.Num_Global_Parts, but ParMETIS3 sometimes
   * does not return the correct number of non-empty partitions, allowing
   * max1 to be less than LB.Num_Global_Parts. 
   * (e.g., ewgt.adaptive-partlocal1-v3.4.?).
   */
  if (max1 > zz->LB.Num_Global_Parts) 
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unexpected value for max1.");

  /* Set up global HG */

  Zoltan_HG_HGraph_Init(&hg);
  if (total_HEcnt) {
    hg.nVtx = max0 + zz->LB.Num_Global_Parts;  
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
  }

  before = measure_stays(zz, &hg, max0, "BEFORE");

  /* Do matching */

  match = (int *) ZOLTAN_CALLOC(hg.nVtx + zz->LB.Num_Global_Parts, sizeof(int));
  used = match + hg.nVtx;
  if (hg.nVtx && !match) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* Max # matches allowed */
  limit = (max0 < zz->LB.Num_Global_Parts ? max0 : zz->LB.Num_Global_Parts); 
  do_match(zz, &hg, match, limit);

      
  /* Build remapping vector, if non-trivial matching was returned. */

  *new_map = 0;
  for (i = 0; i < zz->LB.Num_Global_Parts; i++) 
    if (match[i+max0] != i+max0) {
      *new_map = 1;
      break;
    }

  if (*new_map) {

    zz->LB.Remap = (int *) ZOLTAN_MALLOC(zz->LB.Num_Global_Parts * sizeof(int));
    if (!(zz->LB.Remap)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }


    /* First, process all parts that were matched. Mark matched parts as used.*/

    for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
      zz->LB.Remap[i] = -1; 
      tmp = match[i+max0];
      if (tmp != i+max0) {
        zz->LB.Remap[i] = tmp;
        used[tmp] = 1;
      }
    }

    /* Second, process unmatched parts; if possible, keep same part number. */

    for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
      if (zz->LB.Remap[i] > -1) continue;  /* Already processed part i */
      /* match[i+max0] == i+max0 */
      if (!used[i]) {  /* Keep the same part number if it is not used */
        zz->LB.Remap[i] = i;
        used[i] = 1;
      }
    }
  
    /* Third, process remaining unmatched parts; assign them to 
       unused partitions.*/
  
    for (uidx = 0, i = 0; i < zz->LB.Num_Global_Parts; i++) {
      if (zz->LB.Remap[i] > -1) continue;  /* Already processed part i */
      /* match[i+max0] == i+max0 */
      while (used[uidx]) uidx++;   /* Find next unused partition */
      zz->LB.Remap[i] = uidx;
      used[uidx] = 1;
    }
  }

  if (*new_map) {
    after = measure_stays(zz, &hg, max0, "AFTER ");

    if (before >= after) {
      /* No benefit from remapping; don't keep it! */
      ZOLTAN_FREE(&zz->LB.Remap);
      *new_map = 0;
    }

    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL && zz->Proc == zz->Debug_Proc &&
        zz->LB.Remap) 
      for (i = 0; i < zz->LB.Num_Global_Parts; i++) 
        printf("%d REMAP Part %d to Part %d\n", zz->Proc, i, zz->LB.Remap[i]);
  }

End:
  ZOLTAN_FREE(&match);
  ZOLTAN_FREE(&each_size);
  ZOLTAN_FREE(&recvbuf);
  Zoltan_HG_HGraph_Free(&hg);
  return ierr;
}

/******************************************************************************/
static float measure_stays(
  ZZ *zz,
  HGraph *hg, 
  int max0,
  char *when
)
{
/* Routine that measures and prints the amount of data that doesn't move
 * as described by the hypergraph. 
 */

float stay = 0.;
int tmp, i;

  for (i = 0; i < hg->nEdge; i++) {
    tmp = i + i;
    if (zz->LB.Remap) {
      if (hg->hvertex[tmp] == zz->LB.Remap[hg->hvertex[tmp+1]-max0]) 
        stay += hg->ewgt[i];
    }
    else {
      if (hg->hvertex[tmp] == (hg->hvertex[tmp+1]-max0))
        stay += hg->ewgt[i];
    }
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL && zz->Proc == zz->Debug_Proc)
    printf("%d REMAP--%s: TOTAL AMT STAY = %g\n\n",
            zz->Proc, when, stay);

  return(stay);
}


/******************************************************************************/

/* path growing matching, hypergraph version */
static int matching_pgm (ZZ *zz, HGraph *hg, int *match, int *limit)
{
int i, j, k, side = 0, edge, vertex, *Match[2] = {NULL, NULL};
int limits[2], neighbor, next_vertex, pins;
double w[2]={0.0,0.0}, weight, max_weight, *sims = NULL;
char  *yo = "matching_pgm";

  limits[0] = limits[1] = *limit;
  Match[0] = match;
  if (hg->nVtx) {
    if (!(Match[1] = (int*)   ZOLTAN_MALLOC (hg->nVtx * sizeof(int)))
     || !(sims     = (double*) ZOLTAN_CALLOC (hg->nVtx,  sizeof(double))) ) {
      Zoltan_Multifree (__FILE__, __LINE__, 2, &Match[1], &sims);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      return ZOLTAN_MEMERR;
    }
  }

  for (i = 0; i < hg->nVtx; i++)
    Match[1][i] = i;

  for (i = 0; i < hg->nVtx  &&  limits[side] > 0; i++) {
    if (Match[0][i] == i && Match[1][i] == i) {
      vertex = i;
      while (vertex > 0 && limits[side] > 0) {
        max_weight = 0.0;
        next_vertex = -1;

        for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
          edge = hg->vedge[j];
          pins = hg->hindex[edge+1] - hg->hindex[edge];
          weight = 2.0 / ((pins-1)*pins);
          if (hg->ewgt)
            weight *= hg->ewgt[edge];
          for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; k++) {
            neighbor = hg->hvertex[k];
            if (neighbor != vertex && Match[0][neighbor] == neighbor && 
                Match[1][neighbor]==neighbor)
               sims[neighbor] += weight;
          }
        }
        for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++) {
          edge = hg->vedge[j];
          for (k = hg->hindex[edge]; k < hg->hindex[edge+1]; k++) {
            neighbor = hg->hvertex[k];
            if (sims[neighbor] > 0.0) {
              if (sims[neighbor] > max_weight) {
                max_weight = sims[neighbor];
                next_vertex = neighbor;
              }
              sims[neighbor] = 0.0;
            }
          }
        }

        if (next_vertex >= 0) {
          Match[side][vertex] = next_vertex;
          Match[side][next_vertex] = vertex;
          limits[side]--;
          w[side] += max_weight;
          side = 1-side;
        }
        vertex = next_vertex;
      }
    }
  }

  if (w[0] < w[1]) {
    for (i = 0; i < hg->nVtx; i++)
      match[i] = Match[1][i];
    *limit = limits[1];
  }
  else
    *limit = limits[0];

  Zoltan_Multifree (__FILE__, __LINE__, 2, &Match[1], &sims);
  return ZOLTAN_OK;
}

/******************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
