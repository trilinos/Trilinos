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

static int part_compare(const void *, const void *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Part_To_Proc(ZZ *zz, int part)
{
/* Routine that maps partitions to processors.
 * Assumption:  Given Num_Local_Parts for each processor, assign 
 * 0 to Num_Local_Parts(0)-1 to Proc 0,
 * Num_Local_Parts(0) to Sum[i=0,1:Num_Local_Parts(i)]-1 to Proc 1,
 * ...
 * Sum[i=0,Num_Proc-2:Num_Local_Parts(i)] to
 *     Sum[i=0,Num_Proc-1:Num_Local_Parts(i)]-1 to Proc Num_Proc-1.
 */
char *yo = "Zoltan_LB_Part_To_Proc";
int proc;
int *parray = zz->LB.PartDist;    /* Temporary variable */

  if (zz->LB.Num_Global_Parts == zz->Num_Proc) {
    /*  number of parts == number of procs. return input part. */
    proc = part;
  }
  else {
    /*  number of parts != number of procs */
    if (part < 0 || part >= zz->LB.Num_Global_Parts) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid partition number.");
      proc = -1;
    }
    else {

      /* Check if part is on this proc; often likely; avoid the bsearch */
      if (part >= parray[zz->Proc] && part < parray[zz->Proc+1])
        proc = zz->Proc;
      else
        proc = ((int *)bsearch((void *)&part, parray, zz->Num_Proc+1, 
                                sizeof(int), part_compare)) - parray;
    }
  }
  return proc;
}

/*****************************************************************************/
static int part_compare(const void *key, const void *arr)
{
  if (*((int *)key) < *((int *)arr))
    return -1;
  if (*((int *)key) >= *((int *)(arr)+1))
    return 1;
  return 0;
}

/*****************************************************************************/


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
