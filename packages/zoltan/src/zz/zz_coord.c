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
 * This function gets a list of coordinates one way or the other,
 * i.e., by calling either Get_Geom_Multi or Get_Geom for each object.
 * Degenerate geometry detection and transformations may be added here later.
 */

int Zoltan_Get_Coordinates(
  ZZ *zz, 
  int num_obj,               /* Input:  number of objects */
  ZOLTAN_ID_PTR global_ids,  /* Input:  global IDs of objects */
  ZOLTAN_ID_PTR local_ids,   /* Input:  local IDs of objects; may be NULL. */
  int *num_dim,              /* Output: dimension of coordinates */
  double **coords            /* Output: array of coordinates; malloc'ed by
                                        fn if NULL upon input. */
)
{
  char *yo = "Zoltan_Get_Coordinates";
  int i;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int alloced_coords = 0;
  ZOLTAN_ID_PTR lid;   /* Temporary pointers to local IDs; used to pass 
                          NULL to query functions when NUM_LID_ENTRIES == 0. */
  int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Error check -- make sure needed callback functions are registered. */

  if (zz->Get_Num_Geom == NULL || 
     (zz->Get_Geom == NULL && zz->Get_Geom_Multi == NULL)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register ZOLTAN_NUM_GEOM_FN and "
                       "either ZOLTAN_GEOM_MULTI_FN or ZOLTAN_GEOM_FN");
    goto End;
  }

  /* Get problem dimension. */

  *num_dim = zz->Get_Num_Geom(zz->Get_Num_Geom_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Error returned from ZOLTAN_GET_NUM_GEOM_FN");
    goto End;
  }
  if (*num_dim < 0 || *num_dim > 3) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Invalid dimension returned from ZOLTAN_NUM_GEOM_FN");
    goto End;
  }

  /* Get coordinates for object; allocate memory if not already provided. */

  if (*num_dim > 0 && num_obj > 0) {
    if (*coords == NULL) {
      alloced_coords = 1;
      *coords = (double *) ZOLTAN_MALLOC(num_obj * (*num_dim) * sizeof(double));
      if (*coords == NULL) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error");
        goto End;
      }
    }

    if (zz->Get_Geom_Multi != NULL) {
      zz->Get_Geom_Multi(zz->Get_Geom_Multi_Data, zz->Num_GID, zz->Num_LID,
                         num_obj, global_ids, local_ids, *num_dim, *coords,
                         &ierr);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                           "Error returned from ZOLTAN_GET_GEOM_MULTI_FN");
        goto End;
      }
    }
    else {
      for (i = 0; i < num_obj; i++) {
        lid = (num_lid_entries ? &(local_ids[i*num_lid_entries]) : NULL);
        zz->Get_Geom(zz->Get_Geom_Data, num_gid_entries, num_lid_entries,
                     global_ids + i*num_gid_entries, lid, 
                     (*coords) + i*(*num_dim), &ierr);
        if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                             "Error returned from ZOLTAN_GET_GEOM_FN");
          goto End;
        }
      }
    }
  }

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error found; no coordinates returned.");
    if (alloced_coords) ZOLTAN_FREE(coords);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
