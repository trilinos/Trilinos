/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_rcbutilc_id = "$Id$";
#endif


#include "lb_const.h"
#include "rcb_const.h"
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */

static void initialize_dot(LB *, struct rcb_dot *, LB_ID, LB_ID);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void rcb_build_data_structure(LB *lb, int *num_obj, int *max_obj)
{
/*
 *  Function to build the geometry-based data structures for 
 *  Steve Plimpton's RCB implementation.
 */
char *yo = "rcb_build_data_structure";
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */
LB_ID *objs_global;                   /* Array of global IDs returned by the
                                         application.                        */
LB_ID *objs_local;                    /* Array of local IDs returned by the
                                         application.                        */
LB_ID obj_global_id;                  /* Global ID returned by application.  */
LB_ID obj_local_id;                   /* Local ID returned by application.   */
int num_geom;                         /* # values per object used to describe
                                         the geometry.                       */
int found;
int i;

  /*
   * Allocate an RCB data structure for this load balancing object.
   * If the previous data structure is still there, free the Dots first;
   * the other fields can be reused.
   */

  if (lb->Data_Structure == NULL) {
    rcb = (RCB_STRUCT *) LB_SMALLOC(sizeof(RCB_STRUCT));
    lb->Data_Structure = (void *) rcb;
    rcb->Tree_Ptr = (struct rcb_tree *) LB_array_alloc(__FILE__, __LINE__,
                                                       1, LB_Num_Proc, 
                                                       sizeof(struct rcb_tree));
    rcb->Box = (struct rcb_box *) LB_SMALLOC(sizeof(struct rcb_box));
  }
  else {
    rcb = (RCB_STRUCT *) lb->Data_Structure;
    LB_safe_free((void **) &(rcb->Dots));
  }

  /*
   * Allocate space for objects in RCB data structure.  Allow extra space
   * for objects that are imported to the processor.
   */

  *num_obj = lb->Get_Num_Obj();
  *max_obj = 1.5 * *num_obj;
  rcb->Dots = (struct rcb_dot *) LB_array_alloc(__FILE__, __LINE__, 1, *max_obj,
                                                sizeof(struct rcb_dot));


  /*
   * Compute the number of geometry fields per object.  For RCB, this
   * value should be one, two or three, describing the x-, y-, and z-coords.
   */

  num_geom = lb->Get_Num_Geom();
  if (num_geom > 3) {
    fprintf(stderr, "Error in %s:  Number of geometry fields %d is too great "
                    "for RCB; valid range is 1-3\n", yo, num_geom);
    exit(-1);
  }

  /*
   *  Access objects based on the method provided by the application.
   */

  if (lb->Get_Obj_List != NULL) {

    /*
     *  Call the application for the IDs of all objects and initialize the
     *  dot for each object.
     */

    objs_global = (LB_ID *) LB_array_alloc(__FILE__, __LINE__, 1, 2 * *num_obj,
                                           sizeof(LB_ID));
    objs_local = (LB_ID *) (objs_global + *num_obj);
    lb->Get_Obj_List(objs_global, objs_local);

    for (i = 0; i < *num_obj; i++) {
      initialize_dot(lb, &(rcb->Dots[i]), objs_global[i], objs_local[i]);
    }
    LB_safe_free((void **) &objs_global);
  }
  else if (lb->Get_First_Obj != NULL && lb->Get_Next_Obj != NULL) {

    /*
     *  Call the application for each object and initialize the dot for 
     *  that object.
     */

    i = 0;
    found = lb->Get_First_Obj(&obj_global_id, &obj_local_id);
    while (found) {
      initialize_dot(lb, &(rcb->Dots[i]), obj_global_id, obj_local_id);
      i++;
      found = lb->Get_Next_Obj(obj_global_id, obj_local_id, 
                                     &obj_global_id, &obj_local_id);
    }
    if (i != *num_obj) {
      fprintf(stderr, "Error in %s:  Number of objects returned %d != "
                      "Number of objects declared %d\n", yo, i, *num_obj);
      fprintf(stderr, "Check implementation of LB_FIRST_OBJ_FN and "
                      "LB_NEXT_OBJ_FN \n");
      exit(-1);
    }
  }
  else {
    fprintf(stderr, "Error in %s:  Must define and register either "
                    "LB_OBJ_LIST_FN or LB_FIRST_OBJ_FN/LB_NEXT_OBJ_FN pair\n",
                     yo);
    fprintf(stderr, "Cannot perform RCB without one of these functions.\n");
    exit(-1);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void initialize_dot(LB *lb, struct rcb_dot *dot, LB_ID global_id, 
                           LB_ID local_id)
{
/*
 *  Function that initializes the dot data structure for RCB.  It uses the 
 *  global ID, coordinates and weight provided by the application.  
 */

  dot->Tag.Global_ID = global_id;
  dot->Tag.Local_ID = local_id;
  dot->Tag.Proc = LB_Proc;
  dot->X[0] = dot->X[1] = dot->X[2] = 0.0;
  lb->Get_Geom(global_id, local_id, dot->X);
  if (lb->Get_Obj_Weight != NULL) {
    dot->Weight = lb->Get_Obj_Weight(global_id, local_id);
  }
}
