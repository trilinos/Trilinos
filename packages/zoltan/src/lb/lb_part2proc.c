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
#include "zz_util_const.h"

/* Routines to handle mapping of partitions to processors. */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_LB_Build_ProcDist(ZZ *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Part_To_Proc(ZZ *zz, int part, ZOLTAN_ID_PTR gid)
{
/* Routine that maps partitions to processors.
 * If a partition is entirely within a processor, that processor's rank is
 * returned.
 * If a partition is spread across several processors, find the range of its
 * processors.  If zz->Proc is one of them, return zz->Proc.  Otherwise,
 * hash the input gid to a processor within the range of processors.
 * NOTE:  The special case of returning zz->Proc when it is within range 
 * reduces data movement, but can result in different processor assignments 
 * for the same gid on different processors.
 * If all processors must map a gid to the same processor, this special
 * case must be removed.
 */
char *yo = "Zoltan_LB_Part_To_Proc";
int proc;
int *pdist = zz->LB.PartDist;    /* Temporary variable */
int num_procs_for_part;
int hash_value;

  if (zz->LB.PartDist == NULL) {
    /*  number of parts == number of procs, uniformly distributed. 
     *  return input part. */
    proc = part;
  }
  else if (part >= 0 && part < zz->LB.Num_Global_Parts) {
    /*  number of parts != number of procs or 
     *  non-uniform distribution of parts     */
    num_procs_for_part = pdist[part+1] - pdist[part];
    if (zz->LB.Single_Proc_Per_Part || num_procs_for_part <= 1)
      proc = pdist[part];
    else if (zz->Proc >= pdist[part] && zz->Proc < pdist[part+1])
      /* zz->Proc is in range of procs holding part; return zz->Proc
       * to prevent data movement for exported items.  */
      proc = zz->Proc;
    else {
      /* Map the gid to a processor within range for the partition.
       * Use Zoltan_Hash to attempt to evenly distribute the gids to
       * processors holding the partition. */
      if (gid != NULL) 
        hash_value = Zoltan_Hash(gid, zz->Num_GID, num_procs_for_part);
      else 
        hash_value = 0;
      proc = pdist[part] + hash_value;
    }
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid partition number.");
    proc = -1;
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
char *yo = "Zoltan_LB_Proc_To_Part";
int *partdist = zz->LB.PartDist;
int *procdist;
int ierr = ZOLTAN_OK;
int tmp;
  
  if (proc < 0 || proc >= zz->Num_Proc) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input proc is out of range.");
    ierr = ZOLTAN_FATAL;
    *nparts = 0;
    *fpart = -1;
    goto End;
  }

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
