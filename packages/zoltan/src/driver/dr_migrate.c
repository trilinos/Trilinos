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
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_migrate = "$Id$";
#endif

/*--------------------------------------------------------------------------*/
/* Purpose: Call Zoltan to migrate elements.                                */
/*          Contains all of the callback functions that Zoltan needs        */
/*          for the migration.                                              */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Revision History:                                                        */
/*    10 May 1999:       Date of creation.                                  */
/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "dr_const.h"
#include "dr_err_const.h"
#include "dr_loadbal_const.h"

/*
 *  PROTOTYPES for load-balancer interface functions.
 */
LB_OBJ_SIZE_FN migrate_elem_size;
LB_PACK_OBJ_FN migrate_pack_elem;
LB_UNPACK_OBJ_FN migrate_unpack_elem;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int migrate_elements(
  int Proc,
  ELEM_INFO *elements[],
  struct LB_Struct *lb_obj,
  int num_imp,
  LB_GID *imp_gids,
  LB_LID *imp_lids,
  int *imp_procs,
  int num_exp,
  LB_GID *exp_gids,
  LB_LID *exp_lids,
  int *exp_procs)
{
/* Local declarations. */

/***************************** BEGIN EXECUTION ******************************/

  /*
   * register migration functions
   */
  if (LB_Set_Fn(lb_obj, LB_PRE_MIGRATE_FN_TYPE, NULL, NULL) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }
  if (LB_Set_Fn(lb_obj, LB_OBJ_SIZE_FN_TYPE, (void *) migrate_elem_size, NULL)
        == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }
  if (LB_Set_Fn(lb_obj, LB_PACK_OBJ_FN_TYPE, (void *) migrate_pack_elem,
                (void *) *elements) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }
  if (LB_Set_Fn(lb_obj, LB_UNPACK_OBJ_FN_TYPE, (void *) migrate_unpack_elem,
                (void *) elements) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }

  if (LB_Help_Migrate(lb_obj, num_imp, imp_gids, imp_lids, imp_procs,
                      num_exp, exp_gids, exp_lids, exp_procs) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Help_Migrate()\n");
    return 0;
  }

  return 1;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int migrate_elem_size(void *data, int *ierr)
/*
 * Function to return size of element information for a single element.
 *
 * For now, do not worry about connect tables or adjacency lists
 * or coordinates of the elements or any parts of ELEM_INFO that
 * must be allocated. Just worry about passing the basic struct
 * with the globalID and the rest will be added later.
 */
{
  *ierr = LB_OK;

  return (sizeof(ELEM_INFO));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void migrate_pack_elem(void *data, LB_GID elem_gid, LB_LID elem_lid,
                       int mig_proc, int elem_data_size, char *buf, int *ierr)
{
  ELEM_INFO *elem, *elem_mig;

  if (data == NULL) {
    *ierr = LB_FATAL;
    return;
  }

  elem = (ELEM_INFO *) data; /* this is the head of the element struct array */

  elem_mig = (ELEM_INFO *) buf; /* this is the element struct to be migrated */

  /*
   * for now just copy the easy things
   */
  elem_mig->border   = elem[elem_lid].border;
  elem_mig->globalID = elem[elem_lid].globalID;
  elem_mig->elem_blk = elem[elem_lid].elem_blk;
  elem_mig->cpu_wgt  = elem[elem_lid].cpu_wgt;
  elem_mig->mem_wgt  = elem[elem_lid].mem_wgt;

  /*
   * need to remove this entry from this procs list of elements
   * do so by setting the globalID to -1
   */
  elem[elem_lid].globalID = -1;

  /*
   * need to update the Mesh struct to reflect this element
   * being gone
   */
  Mesh.num_elems--;
  Mesh.eb_cnts[elem[elem_lid].elem_blk]--;
  /*
   * NOTE: it is not worth the effort to determine the change in the
   * number of nodes on this processor until all of the migration is
   * completed.
   */

  *ierr = LB_OK;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void migrate_unpack_elem(void *data, LB_GID elem_gid, int elem_data_size,
                         char *buf, int *ierr)
{
  ELEM_INFO **elem, *tmp, *elem_mig;
  int i;

  if (data == NULL) {
    *ierr = LB_FATAL;
    return;
  }

  /*
   * In this case, data is a pointer to the head of the element
   * struct array. This is so that the array can be reallocated
   * to contain this element.
   */
  elem = (ELEM_INFO **) data;

  elem_mig = (ELEM_INFO *) buf;

  /*
   * check if the element array has any space
   * if not, allocate some new space
   */
  if (Mesh.num_elems == Mesh.elem_array_len) {
    Mesh.elem_array_len += 5;
    tmp = (ELEM_INFO_PTR) realloc (*elem,
                                   Mesh.elem_array_len * sizeof(ELEM_INFO));
    if (tmp == NULL) {
      Gen_Error(0, "fatal: insufficient memory");
      *ierr = LB_FATAL;
      return;
    }
    *elem = tmp;

    /* initialize the new spots */
    for (i = Mesh.num_elems; i < Mesh.elem_array_len; i++)
      (*elem)[i].globalID = -1;
  }

  /* now search for the first empty spot */
  for (i = 0; i < Mesh.elem_array_len; i++)
    if ((*elem)[i].globalID == -1) break;

  /* now put the migrated information into the array */
  (*elem)[i].border   = elem_mig->border;
  (*elem)[i].globalID = elem_mig->globalID;
  (*elem)[i].elem_blk = elem_mig->elem_blk;
  (*elem)[i].cpu_wgt  = elem_mig->cpu_wgt;
  (*elem)[i].mem_wgt  = elem_mig->mem_wgt;

  /* and update the Mesh struct */
  Mesh.num_elems++;
  Mesh.eb_cnts[(*elem)[i].elem_blk]++;

  *ierr = LB_OK;
}
