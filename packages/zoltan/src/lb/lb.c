#include "lb_const.h"
#include "lb.h"
#include "comm.h"
#include "all_allo_const.h"
#include "par.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void perform_error_checking(LB *);
static void help_migrate(LB *);
static void clean_up(LB *);
static void compute_destinations(LB *, int, int, LB_TAG *, int *, LB_TAG **);

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

#ifdef LB_MPI

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

#endif  /* LB_MPI */
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
  lb->Help_Migrate = FALSE;

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

  /*
   *  SET OTHER METHODS HERE!!
   */

  else {  
    fprintf(stderr, "Error from %s:  Invalid LB method specified:  %s\n", 
            yo, method_name);
    exit(-1);
  }

  if (LB_Proc == 0) {
    printf("LB:  Load balancing method = %d (%s)\n", i, method_name);
  }

  /*
   *  Set the parameters pointer if the application specifies parameters.
   */

  if (params != NULL) {
    lb->Params = params;
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
 *    LB *               --  The load balancing object to which this tolerance
 *                           applies.
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

void LB_Balance(LB *lb)
{
/*
 * Main user-call for load-balancing.
 * Input:  a load-balancing object with appropriate function pointers set.
 * Output: ???
 */

int num_objs;                  /* Set to the new number of objects on 
                                  the processor.                            */
int num_keep;                  /* Set to the number of objects the processor
                                  keeps from the old decomposition.         */
int num_non_local;             /* The number of non-local objects in the
                                  processor's new decomposition.            */
LB_TAG *non_local_objs;        /* Array of tags for non-local objects in
                                  the processor's new decomposition.        */
int num_export_objs;           /* The number of local objects that need to
                                  be exported from the processor to establish
                                  the new decomposition.                    */
LB_TAG *export_objs;           /* Array of tags for objects that need to be
                                  exported to establish the new 
                                  decomposition.  */

  perform_error_checking(lb);

  lb->LB_Fn(lb, &num_objs, &num_keep, &num_non_local, &non_local_objs);

  compute_destinations(lb, num_objs, num_non_local, non_local_objs, 
                       &num_export_objs, &export_objs);

  LB_print_sync_start(TRUE);
  {
    int i;
    printf("LBLB: Objects to be exported from Proc %d\n", LB_Proc);
    for (i = 0; i < num_export_objs; i++) {
      printf("    Obj: %10d  Destination: %4d\n", 
             export_objs[i].Global_ID, export_objs[i].Proc);
    }
  }
  LB_print_sync_end(TRUE);

  if (lb->Help_Migrate) {
    help_migrate(lb);
  }
  clean_up(lb);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void perform_error_checking(LB *lb)
{
/* 
 *  Routine to make sure required functions are defined for the given method.
 *  Num_Objs, comm rtns should be defined for all methods.  
 */

}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void compute_destinations(
  LB *lb,                 /* Load balancing object for the current balance. */
  int num_objs,           /* Number of objects originally assigned to the
                             processor.                                     */
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

int *proc_list;          /* List of processors from which objs are to be 
                            imported.                                       */
COMM_OBJ *comm_plan;     /* Communication object returned
                            by Bruce and Steve's communication routines     */
LB_TAG *import_objs;     /* Array of import objects used to request the objs
                            from other processors.                          */
int i;

  /*
   *  Build processor's list of requests for non-local objs.
   */

  proc_list = (int *) array_alloc(1, num_non_local, sizeof(int));
  import_objs = (LB_TAG *) array_alloc(1, num_non_local, sizeof(LB_TAG));

  for (i = 0; i < num_non_local; i++) {
    proc_list[i] = non_local_objs[i].Proc;

    import_objs[i].Global_ID = non_local_objs[i].Global_ID;
    import_objs[i].Local_ID  = non_local_objs[i].Local_ID;
    import_objs[i].Proc      = LB_Proc;
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

  *export_objs = (LB_TAG *) array_alloc(1, *num_export, sizeof(LB_TAG));
  comm_do(comm_plan, (char *) import_objs, sizeof(LB_TAG), 
          (char *) *export_objs);

  LB_safe_free((void **) &proc_list);
  LB_safe_free((void **) &import_objs);
  
  comm_destroy(comm_plan);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void help_migrate(LB *lb)
{

}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void clean_up(LB *lb)
{
/*
 *  Routine to free the load-balancing object's data structures and 
 *  communication data structures.
 */


}
