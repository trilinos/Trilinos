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
#include "irb_const.h"

/* PROTOTYPES */

static int initialize_dot(LB *, struct irb_dot *, LB_GID, LB_LID, int, float);

int LB_IRB_Build_Structure(LB *lb, int *num_obj, int *max_obj, int wgtflag)
{
/* Function to build the geometry-based data structures for IRB method. */
char           *yo = "LB_IRB_Build_Structure";
char           msg[256];
IRB_STRUCT     *irb;                  /* Data structure for IRB.             */
LB_GID         *objs_global;          /* Array of global IDs returned by the
                                         application.                        */
LB_LID         *objs_local;           /* Array of local IDs returned by the
                                         application.                        */
float          *objs_wgt;             /* Array of object weights returned by
                                         the application.                    */
LB_GID         obj_global_id;         /* Global ID returned by application.  */
LB_LID         obj_local_id;          /* Local ID returned by application.   */
float          wgt = 1.;
int            found;
int            i, ierr = 0;

     /* Allocate an IRB data structure for this load balancing structure.
        If the previous data structure is still there, free the Dots first;
        the other fields can be reused. */

     if (lb->Data_Structure == NULL) {
        irb = (IRB_STRUCT *) LB_MALLOC(sizeof(IRB_STRUCT));
        if (irb == NULL) {
           LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
           return(LB_MEMERR);
        }
        lb->Data_Structure = (void *) irb;
        irb->Tree_Ptr = NULL;
        irb->Dots = NULL;

        irb->Tree_Ptr = (struct irb_tree *)
                        LB_MALLOC(lb->Num_Proc* sizeof(struct irb_tree));
        if (irb->Tree_Ptr == NULL) {
           LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
           LB_IRB_Free_Structure(lb);
           return(LB_MEMERR);
        }
     }
     else {
        irb = (IRB_STRUCT *) lb->Data_Structure;
        LB_FREE(&(irb->Dots));
     }

     /* Allocate space for objects in IRB data structure.  Allow extra space
        for objects that are imported to the processor. */

     *num_obj = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
     if (ierr) {
        LB_PRINT_ERROR(lb->Proc, yo, 
                       "Error returned from user function Get_Num_Obj.");
        LB_IRB_Free_Structure(lb);
        return(ierr);
     }
     *max_obj = (int)(1.5 * *num_obj) + 1;
     irb->Dots = (struct irb_dot *)
                 LB_MALLOC((*max_obj)*sizeof(struct irb_dot));
     if (irb->Dots == NULL) {
        LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
        LB_IRB_Free_Structure(lb);
        return(LB_MEMERR);
     }

     /* Compute the number of geometry fields per object.  For IRB, this
        value should be one, two or three, describing the x-, y-, and
        z-coords. */

     irb->Num_Geom = lb->Get_Num_Geom(lb->Get_Num_Geom_Data, &ierr);
     if (irb->Num_Geom > 3 || irb->Num_Geom < 1) {
        sprintf(msg, "Number of geometry fields %d is "
                     "invalid for IRB; valid range is 1-3",
                      irb->Num_Geom);
        LB_PRINT_ERROR(lb->Proc, yo, msg);
        LB_IRB_Free_Structure(lb);
        exit(LB_FATAL);
     }
     if (ierr) {
        LB_PRINT_ERROR(lb->Proc, yo, 
                       "Error returned from user function Get_Num_Geom.");
        LB_IRB_Free_Structure(lb);
        return(ierr);
     }

     /* Access objects based on the method provided by the application. */

     if (lb->Get_Obj_List != NULL) {
        if (*num_obj) {

           /* Call the application for the IDs of all objects and initialize
              the dot for each object. */

           objs_global = (LB_GID *) LB_MALLOC((*num_obj)*sizeof(LB_GID));
           objs_local  = (LB_LID *) LB_MALLOC((*num_obj)*sizeof(LB_LID));
           objs_wgt    = (float  *) LB_MALLOC((*num_obj)*sizeof(float));
           if (objs_global == NULL || objs_local == NULL || objs_wgt == NULL) {
              LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
              LB_FREE(&objs_global);
              LB_FREE(&objs_local);
              LB_FREE(&objs_wgt);
              LB_IRB_Free_Structure(lb);
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
              LB_IRB_Free_Structure(lb);
              return(ierr);
           }

           for (i = 0; i < *num_obj; i++) {
              ierr = initialize_dot(lb, &(irb->Dots[i]), objs_global[i],
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
              LB_IRB_Free_Structure(lb);
              return(ierr);
           }
        }
     }
     else if (lb->Get_First_Obj != NULL && lb->Get_Next_Obj != NULL) {

        /* Call the application for each object and initialize the dot for
           that object. */

        i = 0;
        found = lb->Get_First_Obj(lb->Get_First_Obj_Data, &obj_global_id,
                                  &obj_local_id, wgtflag, &wgt, &ierr);
        if (ierr == LB_FATAL || ierr == LB_MEMERR) {
           LB_PRINT_ERROR(lb->Proc, yo,
                          "Error returned from user function Get_First_Obj.");
           LB_IRB_Free_Structure(lb);
           return(ierr);
        }

        while (found) {
           ierr = initialize_dot(lb, &(irb->Dots[i]), obj_global_id,
                                 obj_local_id, wgtflag, wgt);
           if (ierr == LB_FATAL || ierr == LB_MEMERR) {
              LB_PRINT_ERROR(lb->Proc, yo, 
                             "Error returned from initialize_dot.");
              LB_IRB_Free_Structure(lb);
              return(ierr);
           }
           i++;
           found = lb->Get_Next_Obj(lb->Get_Next_Obj_Data, obj_global_id,
                                    obj_local_id, &obj_global_id,
                                    &obj_local_id, wgtflag, &wgt, &ierr);
           if (ierr == LB_FATAL || ierr == LB_MEMERR) {
              LB_PRINT_ERROR(lb->Proc, yo, 
                             "Error returned from user function Get_Next_Obj.");
              LB_IRB_Free_Structure(lb);
              return(ierr);
           }
        }
        if (i != *num_obj) {
           sprintf(msg, "Number of objects returned %d != "
                        "Number of objects declared %d;"
                        "Check implementation of LB_FIRST_OBJ_FN and "
                        "LB_NEXT_OBJ_FN", 
                         i, *num_obj);
           LB_PRINT_ERROR(lb->Proc, yo, msg);
           LB_IRB_Free_Structure(lb);
           return(LB_FATAL);
        }
     }
     else {
        LB_PRINT_ERROR(lb->Proc, yo, "Must define and register either "
                       "LB_OBJ_LIST_FN or LB_FIRST_OBJ_FN/LB_NEXT_OBJ_FN pair");
        LB_IRB_Free_Structure(lb);
        return(LB_FATAL);
     }

     return(LB_OK);
}

/*****************************************************************************/

void LB_IRB_Free_Structure(LB *lb)
{
/* Deallocate the persistent IRB data structures in lb->Structure.  */
IRB_STRUCT    *irb;                   /* Data structure for IRB. */

     irb = (IRB_STRUCT *) lb->Data_Structure;

     if (irb != NULL) {
        LB_FREE(&(irb->Tree_Ptr));
        LB_FREE(&(irb->Dots));
        LB_FREE(&(lb->Data_Structure));
     }
}

/*****************************************************************************/

static int initialize_dot(LB *lb, struct irb_dot *dot, LB_GID global_id,
                          LB_LID local_id, int wgtflag, float wgt)
{
/* Function that initializes the dot data structure for IRB.  It uses the
   global ID, coordinates and weight provided by the application.  */
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
