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
#include "lb.h"
#include "lb_util_const.h"
#include "comm.h"
#include "all_allo_const.h"
#include "par.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Initialize(int argc, char **argv)
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
   *  Set global processor values for the load balacing tools.
   */

  MPI_Comm_rank(MPI_COMM_WORLD, &LB_Proc);
  MPI_Comm_size(MPI_COMM_WORLD, &LB_Num_Proc);

}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

LB *LB_Create_LB_Object()
{
/*
 *  Function to create a load balancing object.  May want more than one
 *  object if using different decompositions with different techniques.
 *  This function allocates and initializes the object.
 *  Output:
 *    LB *               --  Pointer to a LB object.
 *
 */

char *yo = "LB_Create_LB_Object";
LB *lb;

  /*
   * Allocate storage for the load-balancing object.
   */

  lb = (LB *) LB_smalloc(sizeof(LB));

  /*
   *  Set defaults for fields of lb:
   */

  lb->Method = RCB;    
  lb->LB_Fn = lb_rcb;
  lb->Params = NULL;
  lb->Tolerance = 0.9;
  lb->Data_Structure = NULL;
  lb->Object_Type = 0;

  lb->Get_Obj_Weight = NULL;
  lb->Get_Num_Edges = NULL;
  lb->Get_Edge_List = NULL;
  lb->Get_Num_Geom = NULL;
  lb->Get_Obj_Geom = NULL;
  lb->Get_Num_Local_Obj = NULL;
  lb->Get_All_Local_Objs = NULL;
  lb->Get_Next_Local_Obj = NULL;
  lb->Get_Num_Border_Obj = NULL;
  lb->Get_All_Border_Objs = NULL;
  lb->Get_Next_Border_Obj = NULL;

  lb->Migrate.Help_Migrate = FALSE;
  lb->Migrate.Pack_Obj_Data = NULL;
  lb->Migrate.Unpack_Obj_Data = NULL;
  lb->Migrate.Get_Obj_Data_Size = NULL;
  
  return(lb);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Set_LB_Fn(LB *lb, LB_FN_TYPE fn_type, void *fn())
{
/*
 *  Function to initialize a given LB interface function.
 *  Input:
 *    lb                --  Pointer to a LB object.
 *    fn_type           --  Enum type indicating the function to be set.
 *    fn                --  Pointer to the function to be used in the
 *                          assignment.
 *  Output:
 *    lb                --  Appropriate field set to value in void *().
 */

char *yo = "LB_Set_LB_Fn";

  switch (fn_type) {
  case LB_OBJECT_WEIGHT_FN_TYPE:
    lb->Get_Obj_Weight = (LB_OBJECT_WEIGHT_FN *) fn; 
    break;
  case LB_NUM_EDGES_FN_TYPE:
    lb->Get_Num_Edges = (LB_NUM_EDGES_FN *) fn;
    break;
  case LB_EDGE_LIST_FN_TYPE:
    lb->Get_Edge_List = (LB_EDGE_LIST_FN *) fn;
    break;
  case LB_NUM_GEOM_FN_TYPE:
    lb->Get_Num_Geom = (LB_NUM_GEOM_FN *) fn;
    break;
  case LB_GEOM_FN_TYPE:
    lb->Get_Obj_Geom = (LB_GEOM_FN *) fn;
    break;
  case LB_NUM_OBJ_FN_TYPE:
    lb->Get_Num_Local_Obj = (LB_NUM_OBJ_FN *) fn;
    break;
  case LB_GET_LOCAL_OBJECTS_FN_TYPE:
    lb->Get_All_Local_Objs = (LB_GET_LOCAL_OBJECTS_FN *) fn;
    break;
  case LB_NEXT_OBJ_FN_TYPE:
    lb->Get_Next_Local_Obj = (LB_NEXT_OBJ_FN *) fn;
    break;
  case LB_NUM_BORDER_OBJ_FN_TYPE:
    lb->Get_Num_Border_Obj = (LB_NUM_BORDER_OBJ_FN *) fn;
    break;
  case LB_BORDER_OBJ_FN_TYPE:
    lb->Get_All_Border_Objs = (LB_BORDER_OBJ_FN *) fn;
    break;
  case LB_NEXT_BORDER_OBJ_FN_TYPE:
    lb->Get_Next_Border_Obj = (LB_NEXT_BORDER_OBJ_FN *) fn;
    break;
  case LB_PRE_MIGRATE_FN_TYPE:
    lb->Migrate.Pre_Process = (LB_PRE_MIGRATE_FN *) fn;
    break;
  case LB_OBJECT_SIZE_FN_TYPE:
    lb->Migrate.Get_Obj_Data_Size = (LB_OBJECT_SIZE_FN *) fn;
    break;
  case LB_PACK_OBJECT_FN_TYPE:
    lb->Migrate.Pack_Obj_Data = (LB_PACK_OBJECT_FN *) fn;
    break;
  case LB_UNPACK_OBJECT_FN_TYPE:
    lb->Migrate.Unpack_Obj_Data = (LB_UNPACK_OBJECT_FN *) fn;
    break;
  default:
    fprintf(stderr, "Error from %s:  LB_FN_TYPE %d is invalid.\n", yo, fn_type);
    fprintf(stderr, "Value must be in range 0 to %d\n", LB_MAX_FN_TYPES);
    exit(-1);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Set_LB_Method(LB *lb, char *method_name, double *params)
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

char *yo = "LB_Set_LB_Method";
int i;

  /*
   *  Compare method_name string with standard strings for methods.
   *  If a match is found, set lb->Method and other pointers.
   */

  if (strcasecmp(method_name, "RCB") == 0) {
    lb->Method = RCB;
    lb->LB_Fn = lb_rcb;
/*
    lb->LB_Comm->Build_Request_Proclist = rcb_build_request_proclist;
    lb->LB_Comm->Build_Send_Request_List = rcb_build_send_request_list;
*/
  }
  else if (strcasecmp(method_name, "OCTPART") == 0) {
    lb->Method = OCTPART;
    lb->LB_Fn = oct_init;
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
    exit(-1);
  }

  if (LB_Proc == 0) {
    printf("LB:  Load balancing method = %d (%s)\n", lb->Method, method_name);
  }

  /*
   *  Set the parameters pointer if the application specifies parameters.
   */

  if (params != NULL) {
    lb->Params = (double *) LB_array_alloc(__FILE__, __LINE__, 1,
                                           LB_PARAMS_MAX_SIZE, sizeof(double));
    for (i = 0; i < LB_PARAMS_MAX_SIZE; i++) 
      lb->Params[i] = params[i];
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Set_LB_Tolerance(LB *lb, double tolerance)
{
/*
 *  Function to set the tolerance to which the system must be load balanced.
 *  For example, if the tolerance is set to 0.9, 10% load imbalance between
 *  the most heavily loaded processor and the average load will be accepted
 *  as balanced.
 *  Input:
 *    lb                 --  The load balancing object to which this tolerance
 *                           applies.
 *    tolerance          --  The tolerance desired.
 *  Output:
 *    lb                 --  Tolerance field set to appropriate value.
 */

char *yo = "LB_Set_LB_Tolerance";

  /*
   *  Check tolerance for errors.
   */

  if (tolerance < 0.0 || tolerance > 1.0) {
    fprintf(stderr, "Error from %s:  LB Tolerance is invalid:  %f\n", 
            yo, tolerance);
    fprintf(stderr, "Tolerance must be between 0 and 1.\n");
    exit(-1);
  }

  /*
   *  No error; set the tolerance value.
   */

  lb->Tolerance = tolerance;

  if (LB_Proc == 0) {
    printf("LB:  Load balancing tolerance = %f\n", tolerance);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Set_LB_Object_Type(LB *lb, int object_type)
{
/*
 *  Function to set the object type for objects to be balanced.
 *  The object type is an integer value.  It can be used to help distinguish
 *  between the IDs for, say, elements and surfaces.
 *  This value is used only by the application; it is optional as far
 *  as the load-balancer is concerned.
 *  Input:
 *    LB *               --  The load balancing object for this object type.
 *    int                --  An integer representing the object type.
 *  Output:
 *    LB *               --  Appropriate fields set to designated type.
 */


  lb->Object_Type = object_type;
  if (LB_Proc == 0) {
    printf("LB:  Load balancing object type = %d\n", object_type);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Set_Help_Migrate(struct LB_Struct *lb, int help)
{
/*
 *  Function to set a flag indicating whether the application wants the
 *  load-balancer to help with data migration.   If migration help is
 *  wanted, routines to pack and unpack object data must be provided by
 *  the application (see LB_PACK_OBJECT_FN, LB_UNPACK_OBJECT_FN).
 *
 *  Input:
 *    struct LB_Struct * --  The load balancing object to which this tolerance
 *                           applies.
 *    int                --  TRUE or FALSE to indicate whether the application
 *                           wants migration help.  Default is FALSE.
 *  Output:
 *    struct LB_Struct * --  Appropriate fields set to designated type.
 */
char *yo = "LB_Set_Help_Migrate";

  if (help != TRUE && help != FALSE) {
    fprintf(stderr, "Error from %s:  Invalid value for Help_Migration:  %d\n", 
            yo, help);
    fprintf(stderr, "Value must be between %d or %d.\n", TRUE, FALSE);
    exit(-1);
  }

  lb->Migrate.Help_Migrate = help;
  if (LB_Proc == 0) {
    printf("LB:  Load balancing Help_Migrate flag = %d\n", help);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Balance(
  LB *lb, 
  int *num_import_objs,       /* The number of non-local objects in the
                                 processor's new decomposition.            */
  LB_TAG **import_objs,       /* Array of tags for non-local objects in
                                 the processor's new decomposition.        */
  int *num_export_objs,       /* The number of local objects that need to
                                 be exported from the processor to establish
                                 the new decomposition.                    */
  LB_TAG **export_objs        /* Array of tags for objects that need to be
                                 exported to establish the new 
                                 decomposition.                            */
)
{
/*
 * Main user-call for load-balancing.
 * Input:  a load-balancing object with appropriate function pointers set.
 * Output: ???
 *   Return zero if the decomposition is unchanged by load-balancing;
 *   Return one otherwise.
 */

char *yo = "LB_Balance";
int num_objs;                  /* Set to the new number of objects on 
                                  the processor.                            */
int num_keep;                  /* Set to the number of objects the processor
                                  keeps from the old decomposition.         */
int gmax_imports;              /* Maximum number of imported objects over 
                                  all processors.                           */
double LB_start_time, LB_end_time;
double LB_time[2], LB_max_time[2];

  LB_start_time = MPI_Wtime();

  if (lb->Method == NONE) {
    if (LB_Proc == 0)
      printf("%s Balancing method selected == NONE; no balancing performed\n",
              yo);

    return 0;
  }

  perform_error_checking(lb);

  lb->LB_Fn(lb, &num_objs, &num_keep, num_import_objs, import_objs);

  MPI_Allreduce(num_import_objs, &gmax_imports, 1, MPI_INT, MPI_MAX, 
                MPI_COMM_WORLD);

  if (gmax_imports == 0) {

    /*
     *  Decomposition was not changed by the load balancing; no migration
     *  is needed.
     */

    if (LB_Proc == 0)
      printf("%s No changes to the decomposition due to load-balancing; "
             "no migration is needed.\n", yo);

    return 0;
  }

  compute_destinations(lb, *num_import_objs, *import_objs, 
                       num_export_objs, export_objs);

  LB_end_time = MPI_Wtime();
  LB_time[0] = LB_end_time - LB_start_time;

  if (LB_Debug > 6) {
    int i;
    LB_print_sync_start(TRUE);
    printf("LBLB: Objects to be exported from Proc %d\n", LB_Proc);
    for (i = 0; i < *num_export_objs; i++) {
      printf("    Obj: %10d  Destination: %4d\n", 
             (*export_objs)[i].Global_ID, (*export_objs)[i].Proc);
    }
    LB_print_sync_end(TRUE);
  }

  if (lb->Migrate.Help_Migrate) {
    LB_start_time = MPI_Wtime();
    help_migrate(lb, *num_import_objs, *import_objs, 
                 *num_export_objs, *export_objs);
    LB_end_time = MPI_Wtime();
    LB_time[1] = LB_end_time - LB_start_time;
  }
  
  MPI_Allreduce(LB_time, LB_max_time, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (LB_Proc == 0) {
    printf("DLBLIB LB  Times:  \n");
    printf("DLBLIB     Balance:        %f\n", LB_max_time[0]);
    printf("DLBLIB     HelpMigrate:    %f\n", LB_max_time[1]);
  }

  clean_up(lb);

  return 1;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

