// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "zz_const.h"
#include "rcb.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RCB_Build_Structure(ZZ *zz, int *num_obj, int *max_obj, int wgtflag,
                               double overalloc, int use_ids, int gen_tree)
{
/*
 *  Function to build the geometry-based data structures for 
 *  Steve Plimpton's RCB implementation.
 */
char *yo = "Zoltan_RCB_Build_Structure";
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */
struct rcb_tree *treeptr;
int i, ierr = 0;

  /*
   * Allocate an RCB data structure for this Zoltan structure.
   * If the previous data structure is still there, free the Dots and IDs first;
   * the other fields can be reused.
   */

  if (zz->LB.Data_Structure == NULL) {
    rcb = (RCB_STRUCT *) ZOLTAN_MALLOC(sizeof(RCB_STRUCT));
    if (rcb == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    zz->LB.Data_Structure = (void *) rcb;
    rcb->Tree_Ptr = NULL;
    rcb->Box = NULL;
    rcb->Global_IDs = NULL;
    rcb->Local_IDs = NULL;
    memset(&(rcb->Dots), 0, sizeof(struct Dot_Struct));

    Zoltan_Initialize_Transformation(&(rcb->Tran));

    rcb->Box = (struct rcb_box *) ZOLTAN_MALLOC(sizeof(struct rcb_box));
    if (rcb->Box == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_RCB_Free_Structure(zz);
      return(ZOLTAN_MEMERR);
    }
    if (gen_tree) {
      rcb->Tree_Ptr = (struct rcb_tree *)
        ZOLTAN_CALLOC(zz->LB.Num_Global_Parts, sizeof(struct rcb_tree));
      if (rcb->Tree_Ptr == NULL) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        Zoltan_RCB_Free_Structure(zz);
        return(ZOLTAN_MEMERR);
      }
      /* initialize Tree_Ptr */
      for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
         treeptr = &(rcb->Tree_Ptr[i]);
         /* initialize dim to -1 to prevent use of cut */
         treeptr->dim = -1;
         treeptr->cut = 0.0;
         treeptr->parent = treeptr->left_leaf = treeptr->right_leaf = 0;
      }
    }
  }
  else {
    rcb = (RCB_STRUCT *) zz->LB.Data_Structure;
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    Zoltan_Free_And_Reset_Dot_Structure(&rcb->Dots);
  }

  ierr = Zoltan_RB_Build_Structure(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                               &(rcb->Dots), num_obj, max_obj,
                               &(rcb->Num_Dim),
                               wgtflag, overalloc, use_ids, 1);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Error returned from Zoltan_RB_Build_Structure.");
    Zoltan_RCB_Free_Structure(zz);
    return(ierr);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_RCB_Free_Structure(ZZ *zz)
{
/*
 * Deallocate the persistent RCB data structures in zz->Structure.
 */
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */

  rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);

  if (rcb != NULL) {
    ZOLTAN_FREE(&(rcb->Tree_Ptr));
    ZOLTAN_FREE(&(rcb->Box));
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    Zoltan_Free_And_Reset_Dot_Structure(&rcb->Dots);
    ZOLTAN_FREE(&(zz->LB.Data_Structure));

  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#define COPY_BUFFER(buf, type, num) \
  if (from->buf) { \
    to->buf = (type *)ZOLTAN_MALLOC((num) * sizeof(type)); \
    if (!to->buf) { \
      Zoltan_RCB_Free_Structure(toZZ); \
      ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory."); \
      return(ZOLTAN_MEMERR); \
    } \
    memcpy(to->buf, from->buf, (num) * sizeof(type)); \
  } \
  else { \
    to->buf = NULL; \
  }

int Zoltan_RCB_Copy_Structure(ZZ *toZZ, ZZ const *fromZZ)
{
  char *yo = "Zoltan_RCB_Copy_Structure";
  RCB_STRUCT *to;
  RCB_STRUCT const *from;

  from = (RCB_STRUCT const *)fromZZ->LB.Data_Structure;
  Zoltan_RCB_Free_Structure(toZZ);

  if (!from){
    return(ZOLTAN_OK);
  }

  to = (RCB_STRUCT *)ZOLTAN_MALLOC(sizeof(RCB_STRUCT));
  if (to == NULL) {
    ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  toZZ->LB.Data_Structure = (void *)to;
  *to = *from;

  COPY_BUFFER(Tree_Ptr, struct rcb_tree, fromZZ->LB.Num_Global_Parts);

  COPY_BUFFER(Box, struct rcb_box, 1);

  return ZOLTAN_OK;
}

size_t Zoltan_RCB_Serialize_Structure_Size(ZZ const *zz) {
  size_t size = sizeof(int);
  size += sizeof(struct rcb_tree) * zz->LB.Num_Global_Parts;
  size += sizeof(struct rcb_box);
  size += sizeof(ZZ_Transform);
  return size;
}

void Zoltan_RCB_Serialize_Structure(ZZ const *zz, char **buf)
{
  char *yo = "Zoltan_RCB_Serialize_Structure";
  RCB_STRUCT const *zzrcb = (RCB_STRUCT const *) zz->LB.Data_Structure;

  if (!zzrcb)
    return;

  /* Need only the tree structure for Point_Assign and Box_Assign */
  char *bufptr = *buf;

  *((int *) bufptr) = zzrcb->Num_Dim;
  bufptr += sizeof(int);

  size_t copysize = sizeof(struct rcb_tree) * zz->LB.Num_Global_Parts;
  memcpy(bufptr, (void *)(zzrcb->Tree_Ptr), copysize);
  bufptr += copysize;

  memcpy(bufptr, (void *)(zzrcb->Box), sizeof(struct rcb_box));
  bufptr += sizeof(struct rcb_box);

  memcpy(bufptr, (void *) &(zzrcb->Tran), sizeof(ZZ_Transform));
  bufptr += sizeof(ZZ_Transform);

  *buf = bufptr;
}

void Zoltan_RCB_Deserialize_Structure(ZZ *zz, char **buf)
{
  char *yo = "Zoltan_RCB_Serialize_Structure";
  RCB_STRUCT *zzrcb = (RCB_STRUCT *) ZOLTAN_MALLOC(sizeof(RCB_STRUCT));
  zz->LB.Data_Structure = zzrcb;

  /* initialize as in Zoltan_RCB_Build_Structure */
  zzrcb->Global_IDs = NULL;
  zzrcb->Local_IDs = NULL;
  memset(&(zzrcb->Dots), 0, sizeof(struct Dot_Struct));

  /* Need only the tree structure for Point_Assign and Box_Assign */
  char *bufptr = *buf;

  zzrcb->Num_Dim = *((int *) bufptr);
  bufptr += sizeof(int);

  /* Copy the tree structure */
  size_t copysize = sizeof(struct rcb_tree) * zz->LB.Num_Global_Parts;
  zzrcb->Tree_Ptr = (struct rcb_tree *) ZOLTAN_MALLOC(copysize);
  memcpy((void *)(zzrcb->Tree_Ptr), bufptr, copysize);
  bufptr += copysize;

  /* Copy the box */
  zzrcb->Box = (struct rcb_box *) ZOLTAN_MALLOC(sizeof(struct rcb_box));
  memcpy((void *)(zzrcb->Box), bufptr, sizeof(struct rcb_box));
  bufptr += sizeof(struct rcb_box);

  /* Copy the transformation */
  memcpy((void *) &(zzrcb->Tran), bufptr, sizeof(ZZ_Transform));
  bufptr += sizeof(ZZ_Transform);

  *buf = bufptr;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
