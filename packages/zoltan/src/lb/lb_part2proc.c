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

/* Routines to handle mapping of partitions to processors. */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_LB_Build_ProcDist(ZZ *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Part_To_Proc(ZZ *zz, int part)
{
/* Routine that maps partitions to processors.
 * If a partition is entirely within a processor, that processor's rank is
 * returned.
 * If a partition is spread across several processors, find the range of
 * processors.  If zz->Proc is one of them, return zz->Proc.  Otherwise,
 * return a random processor in the range of processors.
 * (This random feature should be replaced by something smarter to reduce
 * data movement. KDD)
 */
char *yo = "Zoltan_LB_Part_To_Proc";
int proc;
int *pdist = zz->LB.PartDist;    /* Temporary variable */
static int first_time = 1;
int num_procs_for_part;

  if (zz->LB.PartDist == NULL) {
    /*  number of parts == number of procs, uniformly distributed. 
     *  return input part. */
    proc = part;
  }
  else {
    /*  number of parts != number of procs or 
     *  non-uniform distribution of parts     */
    if (part >= 0 && part < zz->LB.Num_Global_Parts) {
      num_procs_for_part = pdist[part+1] - pdist[part];
      if (zz->LB.Single_Proc_Per_Part || num_procs_for_part <= 1)
        proc = pdist[part];
      else if (zz->Proc >= pdist[part] && zz->Proc < pdist[part+1])
        proc = zz->Proc;
      else {
        if (first_time) {
          srand(zz->Proc);
          first_time = 0;
        }
        proc = (rand() % num_procs_for_part) + pdist[part];
      }
    }
    else {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid partition number.");
      proc = -1;
    }
  }
  return proc;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Proc_To_Part(
  ZZ *zz, 
  int proc,       /* Input: processor number */
  int *nparts,    /* Output: Number of partitions on processor proc (>= 0) */
  int *fpart      /* Output: Partition number of first partition on proc. */
)
{
/* Routine that returns the number of partitions and the partition number
 * of the lowest-numbered partition on a given processor.
 * If there are no partitions on a processor, nparts = 0 and fpart = -1.
 */
int *partdist = zz->LB.PartDist;
int *procdist;
int ierr = ZOLTAN_OK;
int tmp;
  
  if (partdist == NULL) {
    *nparts = 1;
    *fpart = proc;
  }

  else {
    if (zz->LB.ProcDist == NULL) {
      ierr = Zoltan_LB_Build_ProcDist(zz);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        *nparts = 0;
        *fpart = -1;
        goto End;
      }
    }
    procdist = zz->LB.ProcDist;
    if (procdist[proc] == -1) {
      *nparts = 0;
      *fpart = -1;
    }
    else {
      tmp = proc+1;
      while (procdist[tmp] == -1) tmp++;
      *nparts = procdist[tmp] - procdist[proc];
      *nparts = ((*nparts < 1) ? 1 : *nparts);
      *fpart = procdist[proc];
    }
  }

End:
  return ierr;
}

/*****************************************************************************/

static int Zoltan_LB_Build_ProcDist(
  ZZ *zz
)
{
/* Routine that computes the inverse of array LB.PartDist.
 * Builds array LB.ProcDist that maps processors to partitions.
 * Entry i of LB.ProcDist is the lowest partition number on processor i. 
 * If processor i has no partitions, ProcDist[i] = -1.
 */
char *yo = "Zoltan_LB_Build_ProcDist";
int ierr = ZOLTAN_OK;
int *partdist = zz->LB.PartDist;
int *procdist;
int i, j;
   
  if (partdist != NULL) {
    procdist = zz->LB.ProcDist 
             = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    if (procdist == NULL) {
      ierr = ZOLTAN_MEMERR;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      goto End;
    }

    for (j = 0, i = 0; i < zz->Num_Proc; i++) {
      if (partdist[j] == i) {
        /* Partition j is on processor i */
        procdist[i] = j;
        while (partdist[j] == i) j++;
      }
      else if (!zz->LB.Single_Proc_Per_Part)
        /* processor i has continuation of previous processor's partition */
        procdist[i] = procdist[i-1];
      else
        /* processor i has no partitions */
        procdist[i] = -1;
    }
    procdist[zz->Num_Proc] = zz->LB.Num_Global_Parts;
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    printf("%d LB.ProcDist: ", zz->Proc);
    for (i = 0; i <= zz->Num_Proc; i++)
      printf("%d ", zz->LB.ProcDist[i]);
    printf("\n");
  }
End:
  return ierr;
}

/*****************************************************************************/


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
