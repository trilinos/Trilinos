// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <ctype.h>
#include "zz_const.h"
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* 
 * This function gets a list of objects one way or the other,
 * i.e., by calling either Get_Obj_List or Get_First_Obj+Get_Next_Obj.
 * It also retrieves object weights and initial partition assignments. 
 *
 * It is called through either Zoltan_Get_Obj_List or
 * Zoltan_Get_Obj_List_Special_Malloc.  In the second case, the
 * special malloc functions required by the Fortran interface are
 * used in creating the global_id, local_id and parts arrays.
 */

static int _Zoltan_Get_Obj_List(
  ZZ *zz, 
  int *num_obj,
  ZOLTAN_ID_PTR *global_ids, 
  ZOLTAN_ID_PTR *local_ids, 
  int wdim, 
  float **objwgts,
  int **parts,
  int special_malloc
)
{
  char *yo = "_Zoltan_Get_Obj_List";
  int i, n;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int alloced_gids = 0, alloced_lids = 0;
  int gid_off, lid_off;
  ZOLTAN_ID_PTR lid, next_lid; /* Temporary pointers to local IDs; used to pass 
                                  NULL to query functions when 
                                  NUM_LID_ENTRIES == 0. */
  float *next_objwgt;          /* Temporarry pointer to an object weight; used
                                  to pass NULL to query functions when wdim=0 */
  int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  *num_obj = 0;
  *objwgts = NULL;

  if (zz->Get_Num_Obj != NULL) {
    *num_obj = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_Obj.");
      goto End;
    }
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register ZOLTAN_NUM_OBJ_FN.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (*num_obj > 0) {

    /* 
     * Test global_ids and local_ids for NULL.  
     * Should be NULL for doing partitioning.
     * Should not be NULL if doing ordering.
     */
    if (special_malloc){
      if (*global_ids == NULL) {
        Zoltan_Special_Malloc(zz, (void **)global_ids,
                          *num_obj, ZOLTAN_SPECIAL_MALLOC_GID);
        alloced_gids = 1;
      }
      if (*local_ids == NULL) {
        Zoltan_Special_Malloc(zz, (void **)local_ids,
                          *num_obj, ZOLTAN_SPECIAL_MALLOC_LID);
        alloced_lids = 1;
      }

      Zoltan_Special_Malloc(zz, (void **)parts,
                          *num_obj, ZOLTAN_SPECIAL_MALLOC_INT);
    }
    else{
      if (*global_ids == NULL) {
        *global_ids = ZOLTAN_MALLOC_GID_ARRAY(zz, *num_obj);
        alloced_gids = 1;
      }
      if (*local_ids == NULL) {
        *local_ids  = ZOLTAN_MALLOC_LID_ARRAY(zz, *num_obj);
        alloced_lids = 1;
      }
      *parts = (int *) ZOLTAN_MALLOC(*num_obj * sizeof(int));
    }

    if (wdim > 0)
      *objwgts  = (float*) ZOLTAN_MALLOC (sizeof(float) * *num_obj * wdim);

    if ((*global_ids == NULL) || (num_lid_entries > 0 && *local_ids == NULL) ||
        (*parts == NULL) ||
        (wdim > 0 && *objwgts == NULL)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory Error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    if (zz->Get_Obj_List != NULL){
      /* Get object list directly */
      zz->Get_Obj_List(zz->Get_Obj_List_Data, 
                       num_gid_entries, num_lid_entries,
                       *global_ids, *local_ids, 
                       wdim, *objwgts, &ierr);
    }
    else if ((zz->Get_First_Obj != NULL) && (zz->Get_Next_Obj != NULL)){
      /* Use iterator functions to loop through object list */
      if (zz->Get_First_Obj(zz->Get_First_Obj_Data, 
                            num_gid_entries, num_lid_entries, 
                            *global_ids, *local_ids, 
                            wdim, *objwgts, &ierr)){
        n = *num_obj;
        i = 0;
        while (!ierr && (i<n-1)){ 
          gid_off = i * num_gid_entries;
          lid_off = i * num_lid_entries;
          lid = (num_lid_entries ? &((*local_ids)[lid_off]) : NULL);
          next_lid = (num_lid_entries ? &((*local_ids)[lid_off+num_lid_entries]) 
                                      : NULL);
          next_objwgt = (wdim ? (*objwgts) + (i+1)*wdim : NULL);
          zz->Get_Next_Obj(zz->Get_Next_Obj_Data, 
                           num_gid_entries, num_lid_entries, 
                           &((*global_ids)[gid_off]), lid, 
                           &((*global_ids)[gid_off+num_gid_entries]),
                           next_lid,
                           wdim, next_objwgt, &ierr);
          i++;
        }
      }
    }
    else { /* No way to get objects */
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register ZOLTAN_OBJ_LIST_FN or "
                         "ZOLTAN_FIRST_OBJ_FN/ZOLTAN_NEXT_OBJ_FN.");
      ierr = ZOLTAN_FATAL;
      goto End;
    }

    /* Get partition information for objects. */
    /* Call user-callback if provided; otherwise, all parts == zz->Proc */
    
    if (zz->Get_Part == NULL && zz->Get_Part_Multi == NULL) {
      for (i = 0; i < *num_obj; i++) 
        (*parts)[i] = zz->Proc;
    }
    else if (zz->Get_Part_Multi) {
      zz->Get_Part_Multi(zz->Get_Part_Multi_Data, 
                              num_gid_entries, num_lid_entries, *num_obj,
                              *global_ids, *local_ids, *parts, &ierr);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                           "Error returned from ZOLTAN_PART_MULTI_FN");
        goto End;
      }
    }
    else {
      for (i = 0; i < *num_obj; i++) {
        lid = (num_lid_entries ? &((*local_ids)[i*num_lid_entries]) : NULL);
        (*parts)[i] = zz->Get_Part(zz->Get_Part_Data,
                              num_gid_entries, num_lid_entries,
                              &((*global_ids)[i*num_gid_entries]), lid, &ierr);
        if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                             "Error returned from ZOLTAN_PART_FN");
          goto End;
        }
      }
    }
  }

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error found; no lists returned.");
    if (special_malloc){
      if (alloced_gids) 
        Zoltan_Special_Free(zz, (void **)global_ids, ZOLTAN_SPECIAL_MALLOC_GID);
      if (alloced_lids) 
        Zoltan_Special_Free(zz, (void **)local_ids, ZOLTAN_SPECIAL_MALLOC_LID);
      Zoltan_Special_Free(zz, (void **)parts, ZOLTAN_SPECIAL_MALLOC_INT);
    }
    else{
      if (alloced_gids) ZOLTAN_FREE(global_ids);
      if (alloced_lids) ZOLTAN_FREE(local_ids);
      ZOLTAN_FREE(parts);
    }
    ZOLTAN_FREE(objwgts);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}
int Zoltan_Get_Obj_List(
  ZZ *zz, 
  int *num_obj,
  ZOLTAN_ID_PTR *global_ids, 
  ZOLTAN_ID_PTR *local_ids, 
  int wdim, 
  float **objwgts,
  int **parts
)
{
  char *yo = "Zoltan_Get_Obj_List";
  int rc=ZOLTAN_OK;
  ZOLTAN_TRACE_ENTER(zz, yo);

  rc = _Zoltan_Get_Obj_List(zz, num_obj, global_ids, local_ids,
               wdim, objwgts, parts, 0);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return rc;
}
int Zoltan_Get_Obj_List_Special_Malloc(
  ZZ *zz, 
  int *num_obj,
  ZOLTAN_ID_PTR *global_ids, 
  ZOLTAN_ID_PTR *local_ids, 
  int wdim, 
  float **objwgts,
  int **parts
)
{
  char *yo = "Zoltan_Get_Obj_List_Special_Malloc";
  int rc=ZOLTAN_OK;
  ZOLTAN_TRACE_ENTER(zz, yo);

  rc = _Zoltan_Get_Obj_List(zz, num_obj, global_ids, local_ids,
               wdim, objwgts, parts, 1);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return rc;
}
/*****************************************************************************/
/*****************************************************************************/

/* 
 * For debugging purposes, print out object information.
 */

int Zoltan_Print_Obj_List(
  ZZ *zz,
  ZOLTAN_ID_PTR Gids,
  ZOLTAN_ID_PTR Lids,
  int wdim,
  float *Weights,
  int *Parts,
  int howMany
)
{
  int len, i, ierr, j, k, rowSize;
  int num_obj=0;
  int lidSize = zz->Num_LID;
  int gidSize = zz->Num_GID;
 
  if (zz->Get_Num_Obj != NULL) {
    num_obj = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  }

  if (howMany <= 0){
    len = num_obj;
  }
  else if (howMany < num_obj){
    len = howMany;
  }
  else{
    len = num_obj;
  }

  if (len > 0){
    if (Gids && gidSize){
      printf("Global ID's \n");
      for(i=0; i<len; i++){
        ZOLTAN_PRINT_GID(zz, Gids);
        Gids += gidSize;
        if (i && (i%10==0)) printf("\n");
      }
      printf("\n");
    }
  
    if (Lids && lidSize){
      printf("Local ID's \n");
      for(i=0; i<len; i++){
        ZOLTAN_PRINT_LID(zz, Lids);
        Lids += lidSize;
        if (i && (i%10==0)) printf("\n");
      }
      printf("\n");
    }
  
    if (Weights && wdim){
      rowSize = 15 / wdim; 
      printf("Object weights\n");
      k=0;
      for(i=0; i<len; i++){
        printf("(");
        for (j=0; j<wdim; j++){
          printf("%f ", Weights[k++]);
        }
        printf(") ");
        if (i && (i%rowSize==0)) printf("\n");
      }
      printf("\n");
    }
  
    if (Parts){
      printf("Partitions\n");
      for(i=0; i<len; i++){
        printf("%d ", Parts[i]); 
        if (i && (i%10==0)) printf("\n");
      }
      printf("\n");
    }
  }

  return len;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
