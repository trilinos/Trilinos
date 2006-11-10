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


#include <stdio.h>
#include <memory.h>
#include "zz_const.h"
#include "zz_rand.h"
#include "params_const.h"

/* Simple partitioning method. */

/*****************************************************************************/
/*  Parameters structure for Simple method. */

static PARAM_VARS Simple_params[] = {
                  { NULL, NULL, NULL, 0 } };
/*****************************************************************************/

int Zoltan_Simple(
  ZZ *zz,                       /* The Zoltan structure.                     */
  float *part_sizes,            /* Input:  Array of size zz->LB.Num_Global_Parts
                                   * zz->Obj_Weight_Dim
                                   containing the percentage of work to be
                                   assigned to each partition.               */
  int *num_import,              /* Return -1. Random uses only export lists. */
  ZOLTAN_ID_PTR *import_global_ids, /* Not used. */
  ZOLTAN_ID_PTR *import_local_ids,  /* Not used. */
  int **import_procs,           /* Not used. */
  int **import_to_part,         /* Not used. */
  int *num_export,              /* Output: Number of objects to export. */
  ZOLTAN_ID_PTR *export_global_ids, /* Output: GIDs to export. */
  ZOLTAN_ID_PTR *export_local_ids,  /* Output: LIDs to export. */
  int **export_procs,           /* Output: Processsors to export to. */
  int **export_to_part          /* Output: Partitions to export to. */
)
{

  /* No import lists computed. */
  *num_import = -1;

  /* Dummy method: No change, so nothing to export. */
  *num_export = 0;

  return ZOLTAN_OK;
}

