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
static char *cvs_lbc_id = "$Id$";
#endif

#include "lb_const.h"
#include "lb_util_const.h"
#include "comm_const.h"
#include "all_allo_const.h"
#include "params_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains all routines implementing the load-balancing interface.
 *  These functions are all callable by the application.  Functions not
 *  intended to be callable by the application should be put in another file,
 *  such as lb_util.c.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Initialize(int argc, char **argv, float *ver)
{
/*
 *  Function to initialize values needed in load balancing tools.
 *  The function should be called after MPI_Init if the application
 *  uses MPI.
 */

int mpi_flag;

  /* 
   *  Test whether MPI is already initialized.  If not, call MPI_Init.
   */

  MPI_Initialized(&mpi_flag);

  if (!mpi_flag) {
    MPI_Init(&argc, &argv);
  }

  /*
   * Now return the version so that the user knows which version of
   * the libarary is being used without having to get the source
   * code.
   */
  *ver = LB_VER;

  return (LB_OK);
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

LB *LB_Create_Object(MPI_Comm communicator)
{
/*
 *  Function to create a load balancing object.  May want more than one
 *  object if using different decompositions with different techniques.
 *  This function allocates and initializes the object.
 *  Output:
 *    LB *               --  Pointer to a LB object.
 *
 */

char *yo = "LB_Create_Object";
LB *lb;
int flag;

  /*
   * Allocate storage for the load-balancing object.
   */

  lb = (LB *) LB_MALLOC(sizeof(LB));
  if (!lb) {
    fprintf(stderr, "Error from %s: Insufficient memory\n", yo);
    return NULL;
  }

  /*
   *  Set MPI values for lb:
   */

  if (communicator == MPI_COMM_NULL) {
    /*
     *  The processor is not in the communicator for the load-balancing
     *  object.  Set lb->Communicator to MPI_COMM_NULL and give dummy 
     *  values to lb->Proc and lb->Num_Proc.
     */
    lb->Communicator = MPI_COMM_NULL;
    lb->Proc = -1;
    lb->Num_Proc = 0;
  }
  else {
    /*
     *  Set Communicator to the communicator passed in.
     */
    MPI_Comm_dup(communicator, &(lb->Communicator));
    MPI_Comm_size(lb->Communicator, &(lb->Num_Proc));
    MPI_Comm_rank(lb->Communicator, &(lb->Proc));
  }

  /*
   *  Set defaults for fields of lb:
   */

  lb->Method = RCB;    
  lb->LB_Fn = LB_rcb;
  lb->Debug = 0;
  lb->Params = NULL;
  lb->Tolerance = 0.9;
  lb->Data_Structure = NULL;

  lb->Get_Num_Edges = NULL;
  lb->Get_Edge_List = NULL;
  lb->Get_Num_Geom = NULL;
  lb->Get_Geom = NULL;
  lb->Get_Num_Obj = NULL;
  lb->Get_Obj_List = NULL;
  lb->Get_First_Obj = NULL;
  lb->Get_Next_Obj = NULL;
  lb->Get_Num_Border_Obj = NULL;
  lb->Get_Border_Obj_List = NULL;
  lb->Get_First_Border_Obj = NULL;
  lb->Get_Next_Border_Obj = NULL;

  lb->Migrate.Help_Migrate = FALSE;
  lb->Migrate.Pre_Process = NULL;
  lb->Migrate.Pack_Obj = NULL;
  lb->Migrate.Unpack_Obj = NULL;
  lb->Migrate.Get_Obj_Size = NULL;
  
  return(lb);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Destroy_Object(LB **lb)
{
/*
 *  Function to free a load balancing object.
 *  Input:
 *    LB *               --  Pointer to a LB object.
 *
 */

char *yo = "LB_Destroy_Object";

  LB_Free_Structure(*lb);

  LB_Free_Params(*lb);

  MPI_Comm_free(&((*lb)->Communicator));

  LB_Free((void **) lb);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Set_Fn(LB *lb, LB_FN_TYPE fn_type, void *fn(), void *data)
{
/*
 *  Function to initialize a given LB interface function.
 *  Input:
 *    lb                --  Pointer to a LB object.
 *    fn_type           --  Enum type indicating the function to be set.
 *    fn                --  Pointer to the function to be used in the
 *                          assignment.
 *    data              --  Pointer to data that the LB library will
 *                          pass as an argument to fn(). May be NULL.
 *  Output:
 *    lb                --  Appropriate field set to value in void *().
 */

char *yo = "LB_Set_Fn";

  switch (fn_type) {
  case LB_NUM_EDGES_FN_TYPE:
    lb->Get_Num_Edges = (LB_NUM_EDGES_FN *) fn;
    lb->Get_Num_Edges_Data = data;
    break;
  case LB_EDGE_LIST_FN_TYPE:
    lb->Get_Edge_List = (LB_EDGE_LIST_FN *) fn;
    lb->Get_Edge_List_Data = data;
    break;
  case LB_NUM_GEOM_FN_TYPE:
    lb->Get_Num_Geom = (LB_NUM_GEOM_FN *) fn;
    lb->Get_Num_Geom_Data = data;
    break;
  case LB_GEOM_FN_TYPE:
    lb->Get_Geom = (LB_GEOM_FN *) fn;
    lb->Get_Geom_Data = data;
    break;
  case LB_NUM_OBJ_FN_TYPE:
    lb->Get_Num_Obj = (LB_NUM_OBJ_FN *) fn;
    lb->Get_Num_Obj_Data = data;
    break;
  case LB_OBJ_LIST_FN_TYPE:
    lb->Get_Obj_List = (LB_OBJ_LIST_FN *) fn;
    lb->Get_Obj_List_Data = data;
    break;
  case LB_FIRST_OBJ_FN_TYPE:
    lb->Get_First_Obj = (LB_FIRST_OBJ_FN *) fn;
    lb->Get_First_Obj_Data = data;
    break;
  case LB_NEXT_OBJ_FN_TYPE:
    lb->Get_Next_Obj = (LB_NEXT_OBJ_FN *) fn;
    lb->Get_Next_Obj_Data = data;
    break;
  case LB_NUM_BORDER_OBJ_FN_TYPE:
    lb->Get_Num_Border_Obj = (LB_NUM_BORDER_OBJ_FN *) fn;
    lb->Get_Num_Border_Obj_Data = data;
    break;
  case LB_BORDER_OBJ_LIST_FN_TYPE:
    lb->Get_Border_Obj_List = (LB_BORDER_OBJ_LIST_FN *) fn;
    lb->Get_Border_Obj_List_Data = data;
    break;
  case LB_FIRST_BORDER_OBJ_FN_TYPE:
    lb->Get_First_Border_Obj = (LB_FIRST_BORDER_OBJ_FN *) fn;
    lb->Get_First_Border_Obj_Data = data;
    break;
  case LB_NEXT_BORDER_OBJ_FN_TYPE:
    lb->Get_Next_Border_Obj = (LB_NEXT_BORDER_OBJ_FN *) fn;
    lb->Get_Next_Border_Obj_Data = data;
    break;
  case LB_PRE_MIGRATE_FN_TYPE:
    lb->Migrate.Pre_Process = (LB_PRE_MIGRATE_FN *) fn;
    lb->Migrate.Pre_Process_Data = data;
    break;
  case LB_OBJ_SIZE_FN_TYPE:
    lb->Migrate.Get_Obj_Size = (LB_OBJ_SIZE_FN *) fn;
    lb->Migrate.Get_Obj_Size_Data = data;
    break;
  case LB_PACK_OBJ_FN_TYPE:
    lb->Migrate.Pack_Obj = (LB_PACK_OBJ_FN *) fn;
    lb->Migrate.Pack_Obj_Data = data;
    break;
  case LB_UNPACK_OBJ_FN_TYPE:
    lb->Migrate.Unpack_Obj = (LB_UNPACK_OBJ_FN *) fn;
    lb->Migrate.Unpack_Obj_Data = data;
    break;
  default:
    fprintf(stderr, "Error from %s:  LB_FN_TYPE %d is invalid.\n", yo, fn_type);
    fprintf(stderr, "Value must be in range 0 to %d\n", LB_MAX_FN_TYPES);
    return (LB_WARN);
  }

  return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Set_Method(LB *lb, char *method_name)
{
/*
 *  Function to set the load balancing method to be used.
 *  Input:
 *    lb                 --  The load balancing object to which this method
 *                           applies.
 *    method_name        --  String specifying the desired method.
 *    params             --  Params needed by desired method.  (This field
 *                           will be better defined later.)
 *  Output:
 *    lbf*               --  Appropriate fields set to designated values.
 */

char *yo = "LB_Set_Method";
int i;

  /*
   *  Compare method_name string with standard strings for methods.
   *  If a match is found, set lb->Method and other pointers.
   *  But first free any left-over data from the previous method.
   */

  LB_Free_Structure(lb);

  if (strcasecmp(method_name, "RCB") == 0) {
    lb->Method = RCB;
    lb->LB_Fn = LB_rcb;
/*
    lb->LB_Comm->Build_Request_Proclist = rcb_build_request_proclist;
    lb->LB_Comm->Build_Send_Request_List = rcb_build_send_request_list;
*/
  }
  else if (strcasecmp(method_name, "OCTPART") == 0) {
    lb->Method = OCTPART;
    lb->LB_Fn = LB_octpart;
  }
  else if (strcasecmp(method_name, "PARMETIS") == 0) {
    lb->Method = PARMETIS;
    lb->LB_Fn = LB_ParMetis;
  }
  else if (strcasecmp(method_name, "NONE") == 0) {
    lb->Method = NONE;
    lb->LB_Fn = NULL;
  }

  /*
   *  SET OTHER METHODS HERE!!
   */

  else {  
    fprintf(stderr, "Error from %s:  Invalid LB method specified:  %s\n", 
            yo, method_name);
    return (LB_FATAL);
  }

  if (lb->Proc == 0) {
    printf("LB:  Load balancing method = %d (%s)\n", lb->Method, method_name);
  }

  return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Set_Migration(struct LB_Struct *lb, int help)
{
/*
 *  Function to set a flag indicating whether the application wants the
 *  load-balancer to help with data migration.   If migration help is
 *  wanted, routines to pack and unpack object data must be provided by
 *  the application (see LB_PACK_OBJ_FN, LB_UNPACK_OBJ_FN).
 *
 *  Input:
 *    struct LB_Struct * --  The load balancing object to which this tolerance
 *                           applies.
 *    int                --  TRUE or FALSE to indicate whether the application
 *                           wants migration help.  Default is FALSE.
 *  Output:
 *    struct LB_Struct * --  Appropriate fields set to designated type.
 */
char *yo = "LB_Set_Migration";

  if (help != TRUE && help != FALSE) {
    fprintf(stderr, "Error from %s:  Invalid value for Help_Migration:  %d\n", 
            yo, help);
    fprintf(stderr, "Value must be between %d or %d.\n", TRUE, FALSE);
    return (LB_FATAL);
  }

  lb->Migrate.Help_Migrate = help;
  if (lb->Proc == 0) {
    printf("LB:  Load balancing Migration flag = %d\n", help);
  }

  return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Balance(
  LB *lb, 
  int *changes,               /* Set to zero or one depending on if the
                                 load balancer determines a new
                                 decomposition or not:
                                 zero - No changes to the decomposition
                                        were made by the load-balancing
                                        algorithm; migration is not needed.
                                 one  - A new decomposition is suggested
                                        by the load-balancer; migration is
                                        needed to establish the new
                                        decomposition.                     */
  int *num_import_objs,       /* The number of non-local objects in the
                                 processor's new decomposition.            */
  LB_GID **import_global_ids, /* Array of global IDs for non-local objects
                                 (i.e., objs to be imported) in
                                 the processor's new decomposition.        */
  LB_LID **import_local_ids,  /* Array of local IDs for non-local objects
                                 (i.e., objs to be imported) in
                                 the processor's new decomposition.        */
  int **import_procs,         /* Array of processor IDs for processors 
                                 currently owning non-local objects (i.e.,
                                 objs to be imported) in this processor's
                                 new decomposition.                        */
  int *num_export_objs,       /* The number of local objects that need to
                                 be exported from the processor to establish
                                 the new decomposition.                    */
  LB_GID **export_global_ids, /* Array of global IDs for objects that need
                                 to be exported (assigned and sent to other
                                 processors) to establish the new 
                                 decomposition.                            */
  LB_LID **export_local_ids,  /* Array of local IDs for objects that need
                                 to be exported (assigned and sent to other
                                 processors) to establish the new 
                                 decomposition.                            */
  int **export_procs          /* Array of destination processor IDs for
                                 objects that need to be exported 
                                 to establish the new decomposition.       */
)
{
/*
 * Main user-call for load-balancing.
 * Input:  a load-balancing object with appropriate function pointers set.
 * Output: 
 *   num_import_objs
 *   import_global_ids
 *   import_local_ids
 *   import_procs
 *   num_export_objs
 *   export_global_ids
 *   export_local_ids
 *   export_procs
 * Return values:
 *   Return zero if the decomposition is unchanged by load-balancing;
 *   Return one otherwise.
 */

char *yo = "LB_Balance";
int gmax_imports;              /* Maximum number of imported objects over 
                                  all processors.                           */
double LB_start_time, LB_end_time;
double LB_time[2] = {0.0,0.0}, LB_max_time[2] = {0.0,0.0};

  LB_start_time = MPI_Wtime();

  /* assume no changes */
  *changes = 0;

  *num_import_objs = *num_export_objs = 0;
  *import_global_ids = NULL;
  *import_local_ids = NULL;
  *import_procs = NULL;
  *export_global_ids = NULL;
  *export_local_ids = NULL;
  *export_procs = NULL;

  /*
   *  Return if this processor is not in the load-balancing object's
   *  communicator.
   */

  if (LB_PROC_NOT_IN_COMMUNICATOR(lb))
    return (LB_OK);

  if (lb->Method == NONE) {
    if (lb->Proc == 0)
      printf("%s Balancing method selected == NONE; no balancing performed\n",
              yo);

    return (LB_WARN);
  }

  /*
   *  Initializations and error checking.
   */

  LB_perform_error_checking(lb);

  lb->LB_Fn(lb, num_import_objs, import_global_ids, import_local_ids,
            import_procs);

  MPI_Allreduce(num_import_objs, &gmax_imports, 1, MPI_INT, MPI_MAX, 
                lb->Communicator);

  if (gmax_imports == 0) {

    /*
     *  Decomposition was not changed by the load balancing; no migration
     *  is needed.
     */

    if (lb->Proc == 0)
      printf("%s No changes to the decomposition due to load-balancing; "
             "no migration is needed.\n", yo);

    return (LB_OK);
  }

  /*
   *  We now know what the new decomposition should look like on the
   *  processor, but we don't know which of our local objects we have
   *  to export to other processors to establish the new decomposition.
   *  Compute the destination map.
   */

  LB_Compute_Destinations(lb, *num_import_objs, *import_global_ids, 
                          *import_local_ids, *import_procs,
                          num_export_objs, export_global_ids,
                          export_local_ids, export_procs);

  LB_end_time = MPI_Wtime();
  LB_time[0] = LB_end_time - LB_start_time;

  if (lb->Debug > 6) {
    int i;
    LB_print_sync_start(lb, TRUE);
    printf("LBLB: Objects to be exported from Proc %d\n", lb->Proc);
    for (i = 0; i < *num_export_objs; i++) {
      printf("    Obj: %10d  Destination: %4d\n", 
             (*export_global_ids)[i], (*export_procs)[i]);
    }
    LB_print_sync_end(lb, TRUE);
  }

  /*
   *  If the Help_Migrate flag is set, perform migration for the application.
   */

  if (lb->Migrate.Help_Migrate) {
    LB_start_time = MPI_Wtime();
    LB_Help_Migrate(lb, *num_import_objs, *import_global_ids,
                    *import_local_ids, *import_procs,
                    *num_export_objs, *export_global_ids,
                    *export_local_ids, *export_procs);
    LB_end_time = MPI_Wtime();
    LB_time[1] = LB_end_time - LB_start_time;
  }
  
  MPI_Allreduce(LB_time, LB_max_time, 2, MPI_DOUBLE, MPI_MAX, lb->Communicator);
  if (lb->Proc == 0) {
    printf("LBLIB LB  Times:  \n");
    printf("LBLIB     Balance:        %f\n", LB_max_time[0]);
    printf("LBLIB     HelpMigrate:    %f\n", LB_max_time[1]);
  }

  *changes = 1;
  return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Compute_Destinations(
  LB *lb,                    /* Load balancing object for current balance.   */
  int num_import,            /* Number of non-local objects assigned to the 
                                processor in the new decomposition.          */
  LB_GID *import_global_ids, /* Array of global IDs for non-local objects 
                                assigned to this processor in the new
                                decomposition.                               */
  LB_LID *import_local_ids,  /* Array of local IDs for non-local objects
                                assigned to the processor in the new
                                decomposition.                               */
  int *import_procs,         /* Array of processor IDs of processors owning
                                the non-local objects that are assigned to
                                this processor in the new decomposition.     */
  int *num_export,           /* Returned value:  Number of objs to be exported
                                to other processors to establish the new
                                decomposition.                               */
  LB_GID **export_global_ids,/* Returned value:  Array of global IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  LB_LID **export_local_ids, /* Returned value:  Array of local IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  int **export_procs         /* Returned value:  Array of processor IDs
                                to which objects will be exported 
                                to establish the new decomposition.          */
)
{
/*
 *  Routine to compute the inverse map:  Given, for each processor, a list 
 *  of non-local objects assigned to the processor, compute the list of objects
 *  that processor needs to export to other processors to establish the new
 *  decomposition.
 */

char *yo = "LB_Compute_Destinations";
int *proc_list = NULL;      /* List of processors from which objs are to be 
                               imported.                                    */
COMM_OBJ *comm_plan;        /* Communication object returned
                               by Bruce and Steve's communication routines  */
LB_TAG *import_objs = NULL; /* Array of import objects used to request objs
                               from other processors.                       */
LB_TAG *export_objs = NULL; /* Array of export objects describing which objs
                               must be sent to other processors.            */
int i;

  /*
   *  Return if this processor is not in the load-balancing object's
   *  communicator.
   */

  if (LB_PROC_NOT_IN_COMMUNICATOR(lb))
    return (LB_OK);

  /*
   *  Build processor's list of requests for non-local objs.
   */

  if (num_import > 0) {
    proc_list = (int *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                       num_import, sizeof(int));
    if (!proc_list) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      return (LB_MEMERR);
    }
    import_objs = (LB_TAG *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                            num_import, sizeof(LB_TAG));
    if (!import_objs) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_Free((void **) &proc_list);
      return (LB_MEMERR);
    }

    for (i = 0; i < num_import; i++) {
      proc_list[i] = import_procs[i];

      LB_SET_GID(import_objs[i].Global_ID, import_global_ids[i]);
      LB_SET_LID(import_objs[i].Local_ID, import_local_ids[i]);
      import_objs[i].Proc = lb->Proc;
    }
  }

  /*
   *  Compute communication map and num_export, the number of objs this
   *  processor has to export to establish the new decomposition.
   */

  comm_plan = LB_comm_create(num_import, proc_list, lb->Communicator,
                             num_export);

  /*
   *  Allocate space for the object tags that need to be exported.  Communicate
   *  to get the list of objects to be exported.
   */

  if (*num_export > 0) {
    export_objs         = (LB_TAG *) LB_Array_Alloc(__FILE__, __LINE__, 1, 
                                                   *num_export, sizeof(LB_TAG));
    if (!export_objs) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_Free((void **) &proc_list);
      LB_Free((void **) &import_objs);
      return (LB_MEMERR);
    }
    *export_global_ids  = (LB_GID *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                                   *num_export, sizeof(LB_GID));
    if (!(*export_global_ids)) { 
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_Free((void **) &proc_list);
      LB_Free((void **) &import_objs);
      LB_Free((void **) &export_objs);
      return (LB_MEMERR);
    }
    *export_local_ids   = (LB_LID *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                                   *num_export, sizeof(LB_LID));
    if (!(*export_local_ids)) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_Free((void **) &proc_list);
      LB_Free((void **) &import_objs);
      LB_Free((void **) &export_objs);
      LB_Free((void **) export_local_ids);
      return (LB_MEMERR);
    }
    *export_procs       = (int *)    LB_Array_Alloc(__FILE__, __LINE__, 1,
                                                   *num_export, sizeof(int));
    if (!(*export_procs)) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_Free((void **) &proc_list);
      LB_Free((void **) &import_objs);
      LB_Free((void **) &export_objs);
      LB_Free((void **) export_local_ids);
      LB_Free((void **) export_procs);
      return (LB_MEMERR);
    }

  }
  else {
    export_objs = NULL;
    *export_global_ids = NULL;
    *export_local_ids = NULL;
    *export_procs = NULL;
  }

  LB_comm_do(comm_plan, (char *) import_objs, sizeof(LB_TAG), 
          (char *) export_objs);

  /*
   *  Put the exported LB_TAGs into the output format.
   */

  for (i = 0; i < *num_export; i++) {
    LB_SET_GID((*export_global_ids)[i], export_objs[i].Global_ID);
    LB_SET_LID((*export_local_ids)[i], export_objs[i].Local_ID);
    (*export_procs)[i]      = export_objs[i].Proc;
  }

  LB_Free((void **) &proc_list);
  LB_Free((void **) &import_objs);
  LB_Free((void **) &export_objs);
  
  LB_comm_destroy(&comm_plan);

  return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Help_Migrate(
  LB *lb,                    /* Load balancing object for current balance.   */
  int num_import,            /* Number of non-local objects assigned to the 
                                processor in the new decomposition.          */
  LB_GID *import_global_ids, /* Array of global IDs for non-local objects 
                                assigned to this processor in the new
                                decomposition; this field can be NULL if 
                                the application does not provide import IDs. */
  LB_LID *import_local_ids,  /* Array of local IDs for non-local objects
                                assigned to the processor in the new
                                decomposition; this field can be NULL if the 
                                application does not provide import IDs.     */
  int *import_procs,         /* Array of processor IDs of processors owning
                                the non-local objects that are assigned to
                                this processor in the new decomposition; this
                                field can be NULL if the application does
                                not provide import IDs.                      */
  int num_export,            /* Number of objs to be exported
                                to other processors to establish the new
                                decomposition.                               */
  LB_GID *export_global_ids, /* Array of global IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  LB_LID *export_local_ids,  /* Array of local IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  int *export_procs          /* Array of processor IDs
                                to which objects will be exported 
                                to establish the new decomposition.          */
)
{
/*
 *  Routine to help perform migration.  If migration pre-processing routine
 *  (LB_PRE_MIGRATE_FN) is specified, this routine first calls that function.
 *  It then calls a function to obtain the size of the migrating objects
 *  (LB_OBJ_SIZE_FN).  The routine next calls an application-specified
 *  object packing routine (LB_PACK_OBJ_FN) for each object
 *  to be exported.  It develops the needed communication map to move the
 *  objects to other processors.  It performs the communication according
 *  to the map, and then calls an application-specified object unpacking
 *  routine (LB_UNPACK_OBJ_FN) for each object imported.
 */

char *yo = "LB_Help_Migrate";
int size;                /* size (in bytes) of the object data for export.  */
char *export_buf = NULL; /* buffer for packing export data.                 */
char *import_buf = NULL; /* buffer for receiving imported data.             */
char *tmp;               /* temporary pointer into buffers.                 */
int i;                   /* loop counter.                                   */
int tmp_import;          /* number of objects to be imported.               */
int *proc_list = NULL;   /* list of processors to which this proc exports.  */
LB_GID global_id;        /* tmp global ID for unpacking objects.            */
COMM_OBJ *comm_plan;     /* Communication object returned
                            by Bruce and Steve's communication routines     */
int ierr = 0;

  /*
   *  Return if this processor is not in the load-balancing object's
   *  communicator.
   */

  if (LB_PROC_NOT_IN_COMMUNICATOR(lb))
    return (LB_OK);

  if (lb->Debug > 4)
    printf("LBLIB %d %s Entering HELP_MIGRATE %d %d\n",
            lb->Proc, yo, num_import, num_export);

  if (lb->Migrate.Get_Obj_Size == NULL) {
    fprintf(stderr, "LBLIB %d %s Error:  Must register an "
           "LB_OBJ_SIZE_FN_TYPE function to use the migration-help tools.\n",
           lb->Proc, yo);
    return (LB_FATAL);
  }

  if (lb->Migrate.Pack_Obj == NULL) {
    fprintf(stderr, "LBLIB %d %s Error:  Must register an "
           "LB_PACK_OBJ_FN_TYPE function to use the migration-help tools.\n",
           lb->Proc, yo);
    return (LB_FATAL);
  }

  if (lb->Migrate.Unpack_Obj == NULL) {
    fprintf(stderr, "LBLIB %d %s Error:  Must register an "
         "LB_UNPACK_OBJ_FN_TYPE function to use the migration-help tools.\n",
         lb->Proc, yo);
    return (LB_FATAL);
  }

  if (lb->Migrate.Pre_Process != NULL) {
    lb->Migrate.Pre_Process(lb->Migrate.Pre_Process_Data,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, &ierr);
    if (ierr) {
      fprintf(stderr, "[%d] %s: Error returned from user defined "
                      "Migrate.Pre_Process function.\n", lb->Proc, yo);
      return (LB_FATAL);
    }

    if (lb->Debug > 5)
      printf("LBLIB %d %s Done Pre-Process\n", lb->Proc, yo);
  }

  size = lb->Migrate.Get_Obj_Size(lb->Migrate.Get_Obj_Size_Data, &ierr);
  if (ierr) {
    fprintf(stderr, "[%d] %s: Error returned from user defined "
                    "Migrate.Get_Obj_Size function.\n", lb->Proc, yo);
    return (LB_FATAL);
  }


  if (num_export > 0) {
    export_buf = (char *) LB_Array_Alloc(__FILE__, __LINE__, 1, num_export,
                                         size);
    if (!export_buf) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      return (LB_FATAL);
    }

    proc_list = (int *) LB_Array_Alloc(__FILE__, __LINE__, 1, num_export,
                                       sizeof(int));
    if (!proc_list) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_Free ((void **) &export_buf);
      return (LB_FATAL);
    }


    /*
     *  Pack the proc_list (to create the map) and the objects for export.
     */
  
    tmp = export_buf;
    for (i = 0; i < num_export; i++) {
      proc_list[i] = export_procs[i];
      lb->Migrate.Pack_Obj(lb->Migrate.Pack_Obj_Data, export_global_ids[i],
                           export_local_ids[i], export_procs[i], size, tmp,
                           &ierr);
      if (ierr) {
        fprintf(stderr, "[%d] %s: Error returned from user defined "
                        "Migrate.Pack_Obj function.\n", lb->Proc, yo);
        return (LB_FATAL);
      }
      tmp += size;
    }
  }

  /*
   *  Compute communication map and tmp_import, the number of objs this
   *  processor has to import to establish the new decomposition.
   */

  comm_plan = LB_comm_create(num_export, proc_list, lb->Communicator,
                             &tmp_import);
  if (tmp_import != num_import) {
    fprintf(stderr, "%d  Error in %s:  tmp_import %d != num_import %d\n", 
            lb->Proc, yo, tmp_import, num_import);
  }

  if (num_import > 0) {
    import_buf = (char *) LB_Array_Alloc(__FILE__, __LINE__, 1, num_import,
                                         size);
    if (!import_buf) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_Free ((void **) &export_buf);
      LB_Free ((void **) &proc_list);
      return (LB_FATAL);
    }

  }

  /*
   *  Send the export data using the communication plan.
   */

  LB_comm_do(comm_plan, export_buf, size, import_buf);

  /*
   *  Free whatever memory we can.
   */

  LB_comm_destroy(&comm_plan);
  LB_Free((void **) &proc_list);
  LB_Free((void **) &export_buf);

  /*
   *  Unpack the object data.
   */

  tmp = import_buf;
  for (i = 0; i < num_import; i++) {
    if (import_global_ids != NULL) 
      LB_SET_GID(global_id, import_global_ids[i]);
    lb->Migrate.Unpack_Obj(lb->Migrate.Unpack_Obj_Data, global_id, size,
                           tmp, &ierr);
    if (ierr) {
      fprintf(stderr, "[%d] %s: Error returned from user defined "
                      "Migrate.Unpack_Obj function.\n", lb->Proc, yo);
      return (LB_FATAL);
    }
    tmp += size;
  }

  LB_Free((void **) &import_buf);
  if (lb->Debug > 4)
    printf("LBLIB %d %s Leaving HELP_MIGRATE %d %d\n",
            lb->Proc, yo, num_import, num_export);

  return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Free_Data(
  LB_GID **import_global_ids, /* Array of global IDs for non-local objects 
                                assigned to this processor in the new
                                decomposition.                               */
  LB_LID **import_local_ids,  /* Array of local IDs for non-local objects
                                assigned to the processor in the new
                                decomposition.                               */
  int **import_procs,         /* Array of processor IDs of processors owning
                                the non-local objects that are assigned to
                                this processor in the new decomposition.     */
  LB_GID **export_global_ids, /* Array of global IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  LB_LID **export_local_ids,  /* Array of local IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  int **export_procs          /* Array of processor IDs
                                to which objects will be exported 
                                to establish the new decomposition.          */
)
{
/*
 *  Routine to free the arrays returning the results of the load balancing.
 */

  LB_Free((void **) import_global_ids);
  LB_Free((void **) import_local_ids);
  LB_Free((void **) import_procs);
  LB_Free((void **) export_global_ids);
  LB_Free((void **) export_local_ids);
  LB_Free((void **) export_procs);

  return (LB_OK);

}
