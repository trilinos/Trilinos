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
#include "comm.h"
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void perform_error_checking(LB *lb)
{
/* 
 *  Routine to make sure required functions are defined for the given method.
 *  Num_Objs, comm rtns should be defined for all methods.  
 */

}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void compute_destinations(
  LB *lb,                 /* Load balancing object for the current balance. */
  int num_non_local,      /* Number of non-local objects assigned to the 
                             processor in the new decomposition.            */
  LB_TAG *non_local_objs, /* Array of tags for non-local objects assigned
                             to the processor in the new decomposition.     */
  int *num_export,        /* Returned value:  Number of objs to be exported
                             to other processors to establish the new
                             decomposition.                                 */
  LB_TAG **export_objs    /* Returned value:  Array of tags of objects to 
                             be exported to other processors to establish the
                             new decomposition.                             */
)
{
/*
 *  Routine to compute the inverse map:  Given, for each processor, a list 
 *  of non-local objects assigned to the processor, compute the list of objects
 *  that processor needs to export to other processors to establish the new
 *  decomposition.
 */

int *proc_list = NULL;      /* List of processors from which objs are to be 
                            imported.                                       */
COMM_OBJ *comm_plan;        /* Communication object returned
                            by Bruce and Steve's communication routines     */
LB_TAG *import_objs = NULL; /* Array of import objects used to request the objs
                            from other processors.                          */
int i;

  /*
   *  Build processor's list of requests for non-local objs.
   */

  if (num_non_local > 0) {
    proc_list = (int *) LB_array_alloc(__FILE__, __LINE__, 1,
                                       num_non_local, sizeof(int));
    import_objs = (LB_TAG *) LB_array_alloc(__FILE__, __LINE__, 1,
                                            num_non_local, sizeof(LB_TAG));

    for (i = 0; i < num_non_local; i++) {
      proc_list[i] = non_local_objs[i].Proc;

      import_objs[i].Global_ID = non_local_objs[i].Global_ID;
      import_objs[i].Local_ID  = non_local_objs[i].Local_ID;
      import_objs[i].Proc      = LB_Proc;
    }
  }

  /*
   *  Compute communication map and num_export, the number of objs this
   *  processor has to export to establish the new decomposition.
   */

  comm_plan = comm_create(num_non_local, proc_list, MPI_COMM_WORLD, num_export);

  /*
   *  Allocate space for the object tags that need to be exported.  Communicate
   *  to get the list of objects to be exported.
   */

  if (*num_export > 0)
    *export_objs = (LB_TAG *) LB_array_alloc(__FILE__, __LINE__, 1,
                                             *num_export, sizeof(LB_TAG));
  else
    *export_objs = NULL;
  comm_do(comm_plan, (char *) import_objs, sizeof(LB_TAG), 
          (char *) *export_objs);

  LB_safe_free((void **) &proc_list);
  LB_safe_free((void **) &import_objs);
  
  comm_destroy(comm_plan);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void help_migrate(LB *lb, int num_import, LB_TAG *import_objs, 
                  int num_export, LB_TAG *export_objs)
{
char *yo = "help_migrate";
int size;                /* size (in bytes) of the object data for export.  */
char *export_buf = NULL; /* buffer for packing export data.                 */
char *import_buf = NULL; /* buffer for receiving imported data.             */
char *tmp;               /* temporary pointer into buffers.                 */
int i;                   /* loop counter.                                   */
int tmp_import;          /* number of objects to be imported.               */
int *proc_list = NULL;   /* list of processors to which this proc exports.  */
COMM_OBJ *comm_plan;     /* Communication object returned
                            by Bruce and Steve's communication routines     */

  
  if (LB_Debug > 4)
    printf("DLBLIB %d %s Entering HELP_MIGRATE %d %d\n",
            LB_Proc, yo, num_import, num_export);

  if (lb->Migrate.Pre_Process != NULL)
    lb->Migrate.Pre_Process(num_import, import_objs, num_export, export_objs);

  if (LB_Debug > 5)
    printf("DLBLIB %d %s Done Pre-Process\n", LB_Proc, yo);

  size = lb->Migrate.Get_Obj_Data_Size(lb->Object_Type);

  if (num_export > 0) {
    export_buf = (char *) LB_array_alloc(__FILE__, __LINE__, 1, num_export,
                                         size);

    proc_list = (int *) LB_array_alloc(__FILE__, __LINE__, 1, num_export,
                                       sizeof(int));

    /*
     *  Pack the proc_list (to create the map) and the objects for export.
     */
  
    tmp = export_buf;
    for (i = 0; i < num_export; i++) {
      proc_list[i] = export_objs[i].Proc;
      lb->Migrate.Pack_Obj_Data(&(export_objs[i]), lb->Object_Type, size, tmp);
      tmp += size;
    }
  }

  /*
   *  Compute communication map and tmp_import, the number of objs this
   *  processor has to import to establish the new decomposition.
   */

  comm_plan = comm_create(num_export, proc_list, MPI_COMM_WORLD, &tmp_import);
  if (tmp_import != num_import) {
    fprintf(stderr, "%d  Error in %s:  tmp_import %d != num_import %d\n", 
            LB_Proc, yo, tmp_import, num_import);
  }

  if (num_import > 0)
    import_buf = (char *) LB_array_alloc(__FILE__, __LINE__, 1, num_import,
                                         size);

  /*
   *  Send the export data using the communication plan.
   */

  comm_do(comm_plan, export_buf, size, import_buf);

  /*
   *  Free whatever memory we can.
   */

  comm_destroy(comm_plan);
  LB_safe_free((void **) &proc_list);
  LB_safe_free((void **) &export_buf);

  /*
   *  Unpack the object data.
   */

  tmp = import_buf;
  for (i = 0; i < num_import; i++) {
    lb->Migrate.Unpack_Obj_Data(size, tmp);
    tmp += size;
  }

  LB_safe_free((void **) &import_buf);
  if (LB_Debug > 4)
    printf("DLBLIB %d %s Leaving HELP_MIGRATE %d %d\n",
            LB_Proc, yo, num_import, num_export);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void clean_up(LB *lb)
{
/*
 *  Routine to free the load-balancing object's data structures and 
 *  communication data structures.
 */


}
