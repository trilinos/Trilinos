/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lbi_const.h"
#include "lb_const.h"
#include "rcb_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */

static int initialize_dot(LB *, struct rcb_dot *, LB_GID, LB_LID, int, float);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RCB_Build_Structure(LB *lb, int *num_obj, int *max_obj, 
                                 int wgtflag)
{
/*
 *  Function to build the geometry-based data structures for 
 *  Steve Plimpton's RCB implementation.
 */
char *yo = "LB_RCB_Build_Structure";
char msg[256];
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */
LB_GID *objs_global;                  /* Array of global IDs returned by the
                                         application.                        */
LB_LID *objs_local;                   /* Array of local IDs returned by the
                                         application.                        */
float  *objs_wgt;                     /* Array of object weights returned by 
                                         the application.                    */
LB_GID obj_global_id;                 /* Global ID returned by application.  */
LB_LID obj_local_id;                  /* Local ID returned by application.   */
int num_geom;                         /* # values per object used to describe
                                         the geometry.                       */
float wgt = 1.;
int found;
int i, ierr = 0;

  /*
   * Allocate an RCB data structure for this load balancing structure.
   * If the previous data structure is still there, free the Dots first;
   * the other fields can be reused.
   */

  if (lb->Data_Structure == NULL) {
    rcb = (RCB_STRUCT *) LB_MALLOC(sizeof(RCB_STRUCT));
    if (rcb == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
    }
    lb->Data_Structure = (void *) rcb;
    rcb->Tree_Ptr = NULL;
    rcb->Box = NULL;
    rcb->Dots = NULL;

    rcb->Tree_Ptr = (struct rcb_tree *)
      LB_MALLOC(lb->Num_Proc* sizeof(struct rcb_tree));
    rcb->Box = (struct rcb_box *) LB_MALLOC(sizeof(struct rcb_box));
    if (rcb->Tree_Ptr == NULL || rcb->Box == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_RCB_Free_Structure(lb);
      return(LB_MEMERR);
    }
    /* initialize dim to -1 to prevent use of cut */
    for (i = 0; i < lb->Num_Proc; i++)
       rcb->Tree_Ptr[i].dim = -1;
  }
  else {
    rcb = (RCB_STRUCT *) lb->Data_Structure;
    LB_FREE(&(rcb->Dots));
  }

  /*
   * Allocate space for objects in RCB data structure.  Allow extra space
   * for objects that are imported to the processor.
   */

  *num_obj = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Num_Obj.");
    LB_RCB_Free_Structure(lb);
    return(ierr);
  }
  *max_obj = (int)(1.5 * *num_obj) + 1;
  rcb->Dots = (struct rcb_dot *) LB_MALLOC((*max_obj)*sizeof(struct rcb_dot));
  if (rcb->Dots == NULL) {
    LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    LB_RCB_Free_Structure(lb);
    return(LB_MEMERR);
  }

  /*
   * Compute the number of geometry fields per object.  For RCB, this
   * value should be one, two or three, describing the x-, y-, and z-coords.
   */

  num_geom = lb->Get_Num_Geom(lb->Get_Num_Geom_Data, &ierr);
  if (num_geom > 3) {
    sprintf(msg, "Number of geometry fields %d is "
                  "too great for RCB; valid range is 1-3\n",
                  num_geom);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_RCB_Free_Structure(lb);
    exit(LB_FATAL);
  }
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Num_Geom.");
    LB_RCB_Free_Structure(lb);
    return(ierr);
  }

  /*
   *  Access objects based on the method provided by the application.
   */

  if (lb->Get_Obj_List != NULL) {
    if (*num_obj) {

      /*
       *  Call the application for the IDs of all objects and initialize the
       *  dot for each object.
       */

      objs_global = (LB_GID *) LB_MALLOC((*num_obj)*sizeof(LB_GID));
      objs_local  = (LB_LID *) LB_MALLOC((*num_obj)*sizeof(LB_LID));
      objs_wgt    = (float  *) LB_MALLOC((*num_obj)*sizeof(float));
      if (objs_global == NULL || objs_local == NULL || objs_wgt == NULL) {
        LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
        LB_FREE(&objs_global);
        LB_FREE(&objs_local);
        LB_FREE(&objs_wgt);
        LB_RCB_Free_Structure(lb);
        return(LB_MEMERR);
      }

      if (wgtflag == 0)
        for (i = 0; i < *num_obj; i++) objs_wgt[i] = 0.;

      lb->Get_Obj_List(lb->Get_Obj_List_Data, objs_global, objs_local, 
                       wgtflag, objs_wgt, &ierr);
      if (ierr == LB_FATAL || ierr == LB_MEMERR) {
        LB_PRINT_ERROR(lb->Proc, yo, 
                       "Error returned from user function Get_Obj_List.");
        LB_FREE(&objs_global);
        LB_FREE(&objs_local);
        LB_FREE(&objs_wgt);
        LB_RCB_Free_Structure(lb);
        return(ierr);
      }

      for (i = 0; i < *num_obj; i++) {
        ierr = initialize_dot(lb, &(rcb->Dots[i]), objs_global[i],
                              objs_local[i], wgtflag, objs_wgt[i]);
        if (ierr == LB_FATAL || ierr == LB_MEMERR) 
          break;
      }
      LB_FREE(&objs_global);
      LB_FREE(&objs_local);
      LB_FREE(&objs_wgt);
      if (ierr == LB_FATAL || ierr == LB_MEMERR) {
        LB_PRINT_ERROR(lb->Proc, yo, 
                       "Error returned from initialize_dot.");
        LB_RCB_Free_Structure(lb);
        return(ierr);
      }
    }
  }
  else if (lb->Get_First_Obj != NULL && lb->Get_Next_Obj != NULL) {

    /*
     *  Call the application for each object and initialize the dot for 
     *  that object.
     */

    i = 0;
    found = lb->Get_First_Obj(lb->Get_First_Obj_Data, &obj_global_id,
                              &obj_local_id, wgtflag, &wgt, &ierr);
    if (ierr == LB_FATAL || ierr == LB_MEMERR) {
      LB_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned from user function Get_First_Obj.");
      LB_RCB_Free_Structure(lb);
      return(ierr);
    }

    while (found) {
      ierr = initialize_dot(lb, &(rcb->Dots[i]), obj_global_id, obj_local_id,
                           wgtflag, wgt);
      if (ierr == LB_FATAL || ierr == LB_MEMERR) {
        LB_PRINT_ERROR(lb->Proc, yo, 
                       "Error returned from initialize_dot.");
        LB_RCB_Free_Structure(lb);
        return(ierr);
      }
      i++;
      found = lb->Get_Next_Obj(lb->Get_Next_Obj_Data, obj_global_id,
                               obj_local_id, &obj_global_id, &obj_local_id,
                               wgtflag, &wgt, &ierr);
      if (ierr == LB_FATAL || ierr == LB_MEMERR) {
        LB_PRINT_ERROR(lb->Proc, yo, 
                       "Error returned from user function Get_Next_Obj.");
        LB_RCB_Free_Structure(lb);
        return(ierr);
      }
    }
    if (i != *num_obj) {
      sprintf(msg, "Number of objects returned %d != "
                   "Number of objects declared %d;"
                   "Check implementation of LB_FIRST_OBJ_FN and LB_NEXT_OBJ_FN",
                   i, *num_obj);
      LB_PRINT_ERROR(lb->Proc, yo, msg);
      LB_RCB_Free_Structure(lb);
      return(LB_FATAL);
    }
  }
  else {
    LB_PRINT_ERROR(lb->Proc, yo, "Must define and register either "
                    "LB_OBJ_LIST_FN or LB_FIRST_OBJ_FN/LB_NEXT_OBJ_FN pair.");
    LB_RCB_Free_Structure(lb);
    return(LB_FATAL);
  }
  return(LB_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_RCB_Free_Structure(LB *lb)
{
/*
 * Deallocate the persistent RCB data structures in lb->Structure.
 */
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */

  rcb = (RCB_STRUCT *) lb->Data_Structure;

  if (rcb != NULL) {
    LB_FREE(&(rcb->Tree_Ptr));
    LB_FREE(&(rcb->Box));
    LB_FREE(&(rcb->Dots));
    LB_FREE(&(lb->Data_Structure));
  }
}



/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int initialize_dot(LB *lb, struct rcb_dot *dot, LB_GID global_id, 
                          LB_LID local_id, int wgtflag, float wgt)
{
/*
 *  Function that initializes the dot data structure for RCB.  It uses the 
 *  global ID, coordinates and weight provided by the application.  
 */
int ierr = LB_OK;
char *yo = "initialize_dot";

  LB_SET_GID(dot->Tag.Global_ID, global_id);
  LB_SET_LID(dot->Tag.Local_ID, local_id);
  dot->Tag.Proc = lb->Proc;
  dot->X[0] = dot->X[1] = dot->X[2] = 0.0;
  lb->Get_Geom(lb->Get_Geom_Data, global_id, local_id, dot->X, &ierr);
  if (ierr == LB_FATAL || ierr == LB_MEMERR) {
    LB_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user defined Get_Geom function.");
    return(ierr);
  }
  if (wgtflag)
     dot->Weight = wgt;
  return(ierr);
}
