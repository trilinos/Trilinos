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
#include "rcb.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RCB_Build_Structure(ZZ *zz, int *num_obj, int *max_obj, int wgtflag,
                           int use_ids)
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
    rcb->Dots = NULL;

    rcb->Tree_Ptr = (struct rcb_tree *)
      ZOLTAN_MALLOC(zz->LB.Num_Global_Parts * sizeof(struct rcb_tree));
    rcb->Box = (struct rcb_box *) ZOLTAN_MALLOC(sizeof(struct rcb_box));
    if (rcb->Tree_Ptr == NULL || rcb->Box == NULL) {
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
  else {
    rcb = (RCB_STRUCT *) zz->LB.Data_Structure;
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    ZOLTAN_FREE(&(rcb->Dots));
  }

  ierr = Zoltan_RB_Build_Structure(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                               &(rcb->Dots), num_obj, max_obj,
                               &(rcb->Num_Dim),
                               wgtflag, use_ids);
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
    ZOLTAN_FREE(&(rcb->Dots));
    ZOLTAN_FREE(&(zz->LB.Data_Structure));
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#define COPY_BUFFER(buf, type, allocLen, copyLen) \
  if (from->buf) { \
    to->buf = (type *)ZOLTAN_MALLOC((allocLen) * sizeof(type)); \
    if (!to->buf) { \
      Zoltan_RCB_Free_Structure(toZZ); \
      ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory."); \
      return(ZOLTAN_MEMERR); \
    } \
    if (copyLen > 0){ \
      memcpy(to->buf, from->buf, (copyLen) * sizeof(type)); \
    } else { \
      memset(to->buf, 0, allocLen); \
    } \
  } \
  else { \
    to->buf = NULL; \
  }


int Zoltan_RCB_Copy_Structure(ZZ *toZZ, ZZ *fromZZ)
{
  char *yo = "Zoltan_RCB_Copy_Structure";
  int num_obj, max_obj, rc;
  RCB_STRUCT *to, *from;
  ZOLTAN_ID_PTR gids, lids;

  from = (RCB_STRUCT *)fromZZ->LB.Data_Structure;
  Zoltan_RCB_Free_Structure(toZZ);

  if (!from){
    return(ZOLTAN_OK);
  }

  rc = Zoltan_Copy_Obj_List(fromZZ, from->Global_IDs, &gids,
    from->Local_IDs, &lids,  0, NULL, NULL, NULL, NULL, &num_obj);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    return rc;
  }

  to = (RCB_STRUCT *)ZOLTAN_MALLOC(sizeof(RCB_STRUCT));
  if (to == NULL) {
    ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&gids);
    ZOLTAN_FREE(&lids);
    return(ZOLTAN_MEMERR);
  }

  memset(to, 0, sizeof(RCB_STRUCT));
  toZZ->LB.Data_Structure = (void *)to;

  max_obj = (int)(1.5 * num_obj) + 1;

  to->Global_IDs = ZOLTAN_REALLOC_GID_ARRAY(fromZZ, gids, max_obj);
  to->Local_IDs = ZOLTAN_REALLOC_LID_ARRAY(fromZZ, lids, max_obj);

  COPY_BUFFER(Dots, struct Dot_Struct, max_obj, num_obj);

  COPY_BUFFER(Tree_Ptr, struct rcb_tree, 
              fromZZ->LB.Num_Global_Parts, fromZZ->LB.Num_Global_Parts);

  COPY_BUFFER(Box, struct rcb_box, 1, 1);

  to->Num_Dim = from->Num_Dim;

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
