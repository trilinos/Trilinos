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
static char *cvs_lbutilc_id = "$Id$";
#endif

#include "lb_const.h"
#include "lb_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_perform_error_checking(LB *lb)
{
/* 
 *  Routine to make sure required functions are defined for the given method.
 *  Num_Objs, comm rtns should be defined for all methods.  
 */

}

/* 
 * This function gets a list of objects one way or the other,
 * i.e., by calling either Get_Obj_List or Get_First_Obj+Get_Next_Obj.
 */

void LB_Get_Obj_List(LB *lb, LB_GID *global_ids, LB_LID *local_ids, 
     int wdim, float *objwgts, int *ierr)
{
  int i, n;

  *ierr = DLB_OK;
  if (lb->Get_Obj_List != NULL){
    /* Get object list directly */
    lb->Get_Obj_List(lb->Get_Obj_List_Data, global_ids, local_ids, 
                     wdim, objwgts, ierr);
  }
  else if ((lb->Get_First_Obj != NULL) && (lb->Get_Next_Obj != NULL)){
    /* Use iterator functions to loop through object list */
    if (lb->Get_First_Obj(lb->Get_First_Obj_Data, global_ids, local_ids, 
        wdim, objwgts, ierr)){
      /* Determine the number of objects since we don't trust the user
         to write the Get_Next_Obj query function in a safe way! */
      n = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, ierr);
      i = 0;
      while (!(*ierr) && (i<n-1)){ 
        lb->Get_Next_Obj(lb->Get_Next_Obj_Data, global_ids[i], 
          local_ids[i], &global_ids[i+1], &local_ids[i+1], 
          wdim, &objwgts[(i+1)*wdim], ierr);
        i++;
      }
    }
  }
  else { /* No way to get objects */
    *ierr = DLB_FATAL;
  }
}

