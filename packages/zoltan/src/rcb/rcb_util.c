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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */

void initialize_dot(LB *, struct rcb_dot *, int, LB_ID);
void rcb_build_data_structure(LB *, int *, int *);

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
LB_ID *objs;                          /* Array of IDs returned by the appl.  */
LB_ID obj_id;                         /* ID returned by the application.     */
int object_type = lb->Object_Type;    /* Object type set by the user.        */
int num_geom;                         /* # values per object used to describe
                                         the geometry.                       */
int i;

  /*
   * Allocate an RCB data structure for this load balancing object.
   * If the previous data structure is still there, free the Dots first;
   * the other fields can be reused.
   */

  if (lb->Data_Structure == NULL) {
    rcb = (RCB_STRUCT *) LB_smalloc(sizeof(RCB_STRUCT));
    lb->Data_Structure = (void *) rcb;
    rcb->Tree_Ptr = (struct rcb_tree *) LB_array_alloc(1, LB_Num_Proc, 
                                                       sizeof(struct rcb_tree));
    rcb->Box = (struct rcb_box *) LB_smalloc(sizeof(struct rcb_box));
  }
  else {
    rcb = (RCB_STRUCT *) lb->Data_Structure;
    safe_free((void **) &(rcb->Dots));
  }

  /*
   * Allocate space for objects in RCB data structure.  Allow extra space
   * for objects that are imported to the processor.
   */

  *num_obj = lb->Get_Num_Local_Obj(object_type);
  *max_obj = 1.5 * *num_obj;
  rcb->Dots = (struct rcb_dot *) LB_array_alloc(1, *max_obj,
                                                sizeof(struct rcb_dot));


  /*
   * Compute the number of geometry fields per object.  For RCB, this
   * value should be one, two or three, describing the x-, y-, and z-coords.
   */

  num_geom = lb->Get_Num_Geom(object_type);
  if (num_geom > 3) {
    fprintf(stderr, "Error in %s:  Number of geometry fields %d is too great "
                    "for RCB; valid range is 1-3\n", yo, num_geom);
    exit(-1);
  }

  /*
   *  Access objects based on the method provided by the application.
   */

  if (lb->Get_All_Local_Objs != NULL) {

    /*
     *  Call the application for the IDs of all objects and initialize the
     *  dot for each object.
     */

    objs = (LB_ID *) LB_array_alloc(1, *num_obj, sizeof(LB_ID));
    lb->Get_All_Local_Objs(object_type, objs);

    for (i = 0; i < *num_obj; i++) {
      initialize_dot(lb, &(rcb->Dots[i]), i, objs[i]);
    }
    LB_safe_free((void **) &objs);
  }
  else if (lb->Get_Next_Local_Obj != NULL) {

    /*
     *  Call the application for each object and initialize the dot for 
     *  that object.
     */

    for (i = 0, obj_id = NULL; i < *num_obj; i++) {
      obj_id = lb->Get_Next_Local_Obj(obj_id, object_type);
      initialize_dot(lb, &(rcb->Dots[i]), i, obj_id);
    }
  }
  else {
    fprintf(stderr, "Error in %s:  Must define and register either function "
                    "Get_Next_Local_Obj or Get_All_Local_Objs\n", yo);
    fprintf(stderr, "Cannot perform RCB without one of these functions.\n");
    exit(-1);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void initialize_dot(LB *lb, struct rcb_dot *dot, int local_id, LB_ID global_id)
{
/*
 *  Function that initializes the dot data structure for RCB.  It uses the 
 *  global ID, coordinates and weight provided by the application.  
 */

int object_type = lb->Object_Type;

  dot->Tag.Local_ID = local_id;
  dot->Tag.Global_ID = global_id;
  dot->Tag.Proc = LB_Proc;
  dot->X[0] = dot->X[1] = dot->X[2] = 0.0;
  lb->Get_Obj_Geom(global_id, object_type, dot->X);
  if (lb->Get_Obj_Weight != NULL)
    dot->Weight = lb->Get_Obj_Weight(global_id, object_type);
}
