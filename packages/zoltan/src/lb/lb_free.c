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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for freeing arrays allocated by Zoltan and
 *  returned to the application; these functions are all callable by the 
 *  application.  
 *
 *  Also includes routine for freeing memory in zz->LB (LB_Struct).
 *  This routine should be called only by Zoltan.
 */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_LB_Free_Part(
  ZOLTAN_ID_PTR *global_ids, /* Array of global IDs */
  ZOLTAN_ID_PTR *local_ids,  /* Array of local IDs */
  int **procs,               /* Array of processor IDs */
  int **to_part              /* Array of partition assignments */
)
{
/*
 *  Routine to free the arrays returning the results of the load balancing.
 */

  ZOLTAN_FREE(global_ids);
  ZOLTAN_FREE(local_ids);
  ZOLTAN_FREE(procs);
  ZOLTAN_FREE(to_part);

  return (ZOLTAN_OK);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Free_Data(
  ZOLTAN_ID_PTR *import_global_ids, /* Array of global IDs for non-local objects 
                                    assigned to this processor in the new
                                    decomposition.                           */
  ZOLTAN_ID_PTR *import_local_ids,  /* Array of local IDs for non-local objects
                                    assigned to the processor in the new
                                    decomposition.                           */
  int **import_procs,           /* Array of processor IDs of processors owning
                                   the non-local objects that are assigned to
                                   this processor in the new decomposition.  */
  ZOLTAN_ID_PTR *export_global_ids, /* Array of global IDs of
                                   objects to be exported to other processors
                                   to establish the new decomposition.       */
  ZOLTAN_ID_PTR *export_local_ids,  /* Array of local IDs of
                                   objects to be exported to other processors
                                   to establish the new decomposition.       */
  int **export_procs            /* Array of processor IDs
                                   to which objects will be exported 
                                   to establish the new decomposition.       */
)
{
/*
 *  Routine to free the arrays returning the results of the load balancing.
 */

  Zoltan_LB_Free_Part(import_global_ids, import_local_ids, import_procs, NULL);
  Zoltan_LB_Free_Part(export_global_ids, export_local_ids, export_procs, NULL);

  return (ZOLTAN_OK);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_LB_Free_Struct(struct Zoltan_LB_Struct *lb)
{
  ZOLTAN_FREE((void **) &(lb->PartDist));
  ZOLTAN_FREE((void **) &(lb->ProcDist));
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
