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

#include "lb_const.h"
#include "lb_util_const.h"
#include "comm_const.h"
#include "all_allo_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "ha_const.h"
#include <ctype.h>

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

LB *LB_Create(MPI_Comm communicator)
{
/*
 *  Function to create a load balancing structure.  May want more than one
 *  structure if using different decompositions with different techniques.
 *  This function allocates and initializes the structure.
 *  Output:
 *    LB *               --  Pointer to a LB structure.
 *
 */

char *yo = "LB_Create_Object";
LB *lb;

  /*
   * Allocate storage for the load-balancing structure.
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
     *  structure.  Set lb->Communicator to MPI_COMM_NULL and give dummy 
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
  lb->Debug_Level = 0;
  lb->Fortran = 0;
  lb->Machine_Desc = NULL;
  lb->Params = NULL;
  lb->Imbalance_Tol = 1.1;
  lb->Deterministic = TRUE;
  lb->Obj_Weight_Dim = 0;
  lb->Comm_Weight_Dim = 0;
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
  lb->Get_Num_Coarse_Obj = NULL;
  lb->Get_Coarse_Obj_List = NULL;
  lb->Get_First_Coarse_Obj = NULL;
  lb->Get_Next_Coarse_Obj = NULL;
  lb->Get_Num_Child = NULL;
  lb->Get_Child_List = NULL;
  lb->Get_Child_Weight = NULL;

  lb->Get_Num_Edges_Fort = NULL;
  lb->Get_Edge_List_Fort = NULL;
  lb->Get_Num_Geom_Fort = NULL;
  lb->Get_Geom_Fort = NULL;
  lb->Get_Num_Obj_Fort = NULL;
  lb->Get_Obj_List_Fort = NULL;
  lb->Get_First_Obj_Fort = NULL;
  lb->Get_Next_Obj_Fort = NULL;
  lb->Get_Num_Border_Obj_Fort = NULL;
  lb->Get_Border_Obj_List_Fort = NULL;
  lb->Get_First_Border_Obj_Fort = NULL;
  lb->Get_Next_Border_Obj_Fort = NULL;
  lb->Get_Num_Coarse_Obj_Fort = NULL;
  lb->Get_Coarse_Obj_List_Fort = NULL;
  lb->Get_First_Coarse_Obj_Fort = NULL;
  lb->Get_Next_Coarse_Obj_Fort = NULL;
  lb->Get_Num_Child_Fort = NULL;
  lb->Get_Child_List_Fort = NULL;
  lb->Get_Child_Weight_Fort = NULL;

  lb->Migrate.Auto_Migrate = FALSE;
  lb->Migrate.Pre_Process = NULL;
  lb->Migrate.Post_Process = NULL;
  lb->Migrate.Pack_Obj = NULL;
  lb->Migrate.Unpack_Obj = NULL;
  lb->Migrate.Get_Obj_Size = NULL;
  
  lb->Migrate.Pre_Process_Fort = NULL;
  lb->Migrate.Pack_Obj_Fort = NULL;
  lb->Migrate.Unpack_Obj_Fort = NULL;
  lb->Migrate.Get_Obj_Size_Fort = NULL;

  return(lb);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Destroy(LB **lb)
{
/*
 *  Function to free a load balancing structure.
 *  Input:
 *    LB *               --  Pointer to a LB structure.
 *
 */

  LB_Free_Structure(*lb);

  LB_Free_Params(*lb);

  MPI_Comm_free(&((*lb)->Communicator));

  LB_FREE(lb);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Set_Fn(LB *lb, LB_FN_TYPE fn_type, void *fn, void *data)
{
/*
 *  Function to initialize a given LB interface function.
 *  Input:
 *    lb                --  Pointer to a LB structure.
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
  case LB_POST_MIGRATE_FN_TYPE:
    lb->Migrate.Post_Process = (LB_POST_MIGRATE_FN *) fn;
    lb->Migrate.Post_Process_Data = data;
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
  case LB_NUM_COARSE_OBJ_FN_TYPE:
    lb->Get_Num_Coarse_Obj = (LB_NUM_COARSE_OBJ_FN *) fn;
    lb->Get_Num_Coarse_Obj_Data = data;
    break;
  case LB_COARSE_OBJ_LIST_FN_TYPE:
    lb->Get_Coarse_Obj_List = (LB_COARSE_OBJ_LIST_FN *) fn;
    lb->Get_Coarse_Obj_List_Data = data;
    break;
  case LB_FIRST_COARSE_OBJ_FN_TYPE:
    lb->Get_First_Coarse_Obj = (LB_FIRST_COARSE_OBJ_FN *) fn;
    lb->Get_First_Coarse_Obj_Data = data;
    break;
  case LB_NEXT_COARSE_OBJ_FN_TYPE:
    lb->Get_Next_Coarse_Obj = (LB_NEXT_COARSE_OBJ_FN *) fn;
    lb->Get_Next_Coarse_Obj_Data = data;
    break;
  case LB_NUM_CHILD_FN_TYPE:
    lb->Get_Num_Child = (LB_NUM_CHILD_FN *) fn;
    lb->Get_Num_Child_Data = data;
    break;
  case LB_CHILD_LIST_FN_TYPE:
    lb->Get_Child_List = (LB_CHILD_LIST_FN *) fn;
    lb->Get_Child_List_Data = data;
    break;
  case LB_CHILD_WEIGHT_FN_TYPE:
    lb->Get_Child_Weight = (LB_CHILD_WEIGHT_FN *) fn;
    lb->Get_Child_Weight_Data = data;
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
 *    lb                 --  The load balancing structure to which this method
 *                           applies.
 *    method_name        --  String specifying the desired method.
 *
 *  Output:
 *    lbf*               --  Appropriate fields set to designated values.
 */

  char *yo = "LB_Set_Method";
  char *method_upper;
  int i;

  /*
   *  Compare method_name string with standard strings for methods.
   *  If a match is found, set lb->Method and other pointers.
   *  But first free any left-over data from the previous method.
   */

  LB_Free_Structure(lb);

  /*
   *  Convert method_name to all upper case.
   *  Do not change the original string.
   */
  method_upper = (char *)LB_MALLOC((strlen(method_name)+1)*sizeof(char));
  for (i=strlen(method_name); i>=0; i--){
      method_upper[i] = toupper(method_name[i]);
  }

  if (strcmp(method_upper, "RCB") == 0) {
    lb->Method = RCB;
    lb->LB_Fn = LB_rcb;
  }
  else if (strcmp(method_upper, "OCTPART") == 0) {
    lb->Method = OCTPART;
    lb->LB_Fn = LB_octpart;
  }
  else if (strcmp(method_upper, "PARMETIS") == 0) {
    lb->Method = PARMETIS;
    lb->LB_Fn = LB_ParMetis;
  }
  else if (strcmp(method_upper, "JOSTLE") == 0) {
    lb->Method = JOSTLE;
    lb->LB_Fn = LB_Jostle;
  }
  else if (strcmp(method_upper, "REFTREE") == 0) {
    lb->Method = REFTREE;
    lb->LB_Fn = LB_Reftree_Part;
  }
  else if (strcmp(method_upper, "NONE") == 0) {
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
    printf("ZOLTAN Load balancing method = %d (%s)\n", lb->Method, method_name);
  }

  LB_FREE(&method_upper);

  return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Balance(
  LB *lb, 
  int *changes,               /* Set to zero or one depending on if 
                                 Zoltan determines a new
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
 * Input:  a load-balancing structure with appropriate function pointers set.
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
int gmax;    /* Maximum number of imported/exported objects 
                over all processors.                       */
int error;    /* Error code */
double start_time, end_time;
double lb_time[2] = {0.0,0.0};

  if (lb->Debug_Level > 0 && lb->Proc == 0)
    LB_Print_Key_Params(lb);

  start_time = LB_Time();

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
   *  Return if this processor is not in the load-balancing structure's
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


  /*
   *  Construct the heterogenous machine description.
   */

  error = LB_Build_Machine_Desc(lb);

  if (error == LB_FATAL){
    return (error);
  }

  /*
   * Call the actual load-balancing function.
   */

  error = lb->LB_Fn(lb, num_import_objs, import_global_ids, import_local_ids,
          import_procs, num_export_objs, export_global_ids, 
          export_local_ids, export_procs);

  if (error == LB_FATAL){
    printf("[%1d] FATAL ERROR: Zoltan returned error code %d\n", 
      lb->Proc, error);
    return (error);
  }
  else if (error){
    printf("[%1d] WARNING: Zoltan returned error code %d\n", 
      lb->Proc, error);
  }

  if (*num_import_objs >= 0)
    MPI_Allreduce(num_import_objs, &gmax, 1, MPI_INT, MPI_MAX, 
                lb->Communicator);
  else /* use export data */
    MPI_Allreduce(num_export_objs, &gmax, 1, MPI_INT, MPI_MAX, 
                lb->Communicator);

  if (gmax == 0) {

    /*
     *  Decomposition was not changed by the load balancing; no migration
     *  is needed.
     */

    if (lb->Proc == 0)
      printf("%s No changes to the decomposition due to load-balancing; "
             "no migration is needed.\n", yo);

    /*
     *  Reset num_import_objs and num_export_objs; don't want to return
     *  -1 for the arrays that weren't returned by LB_FN.
     */

    *num_import_objs = *num_export_objs = 0;

    return (LB_OK);
  }

  /*
   *  Check whether we know the import data, export data, or both.
   *
   *  If we were given the import data,
   *  we know what the new decomposition should look like on the
   *  processor, but we don't know which of our local objects we have
   *  to export to other processors to establish the new decomposition.
   *  Reverse the argument if we were given the export data.
   *
   *  Unless we were given both maps, compute the inverse map.
   */

  if (*num_import_objs >= 0){
    if (*num_export_objs >= 0)
      /* Both maps already available; nothing to do. */;
    else
      /* Compute export map */
      LB_Compute_Destinations(lb, *num_import_objs, *import_global_ids, 
                              *import_local_ids, *import_procs,
                              num_export_objs, export_global_ids,
                              export_local_ids, export_procs);
  }
  else { /* if (*num_import_objs < 0) */
    if (*num_export_objs >= 0)
      /* Compute export map */
      LB_Compute_Destinations(lb, *num_export_objs, *export_global_ids, 
                              *export_local_ids, *export_procs,
                              num_import_objs, import_global_ids,
                              import_local_ids, import_procs);
    else{
      /* No map at all available */
      printf("%s Error: Load-balancing function returned neither import nor "
             "export data.\n", yo);
      return LB_WARN;
    }
  }

  end_time = LB_Time();
  lb_time[0] = end_time - start_time;

  if (lb->Debug_Level > 6) {
    int i;
    LB_Print_Sync_Start(lb, TRUE);
    printf("ZOLTAN: Objects to be exported from Proc %d\n", lb->Proc);
    for (i = 0; i < *num_export_objs; i++) {
      printf("    Obj: %10d  Destination: %4d\n", 
             (*export_global_ids)[i], (*export_procs)[i]);
    }
    LB_Print_Sync_End(lb, TRUE);
  }

  /*
   *  If the Help_Migrate flag is set, perform migration for the application.
   */

  if (lb->Migrate.Auto_Migrate) {
    start_time = LB_Time();
    LB_Help_Migrate(lb, *num_import_objs, *import_global_ids,
                    *import_local_ids, *import_procs,
                    *num_export_objs, *export_global_ids,
                    *export_local_ids, *export_procs);
    end_time = LB_Time();
    lb_time[1] = end_time - start_time;
  }
  
  /* Print timing info */
  if (lb->Proc == 0) {
    printf("ZOLTAN Times:  \n");
  }
  LB_Print_Stats (lb, lb_time[0], "ZOLTAN     Balance:     ");
  LB_Print_Stats (lb, lb_time[1], "ZOLTAN     HelpMigrate: ");

  *changes = 1;
  if (error)
    return (error);
  else
    return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Compute_Destinations(
  LB *lb,                    /* Load balancing structure for current balance.*/
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
COMM_OBJ *comm_plan;        /* Object returned communication routines  */
LB_TAG *import_objs = NULL; /* Array of import objects used to request objs
                               from other processors.                       */
LB_TAG *export_objs = NULL; /* Array of export objects describing which objs
                               must be sent to other processors.            */
int msgtag, msgtag2;        /* Message tags for communication routines */
int i;


  /*
   *  Return if this processor is not in the load-balancing structure's
   *  communicator.
   */

  if (LB_PROC_NOT_IN_COMMUNICATOR(lb))
    return (LB_OK);

  /*
   *  Build processor's list of requests for non-local objs.
   */

  if (num_import > 0) {
    proc_list = (int *) LB_MALLOC(num_import*sizeof(int));
    if (!proc_list) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      return (LB_MEMERR);
    }
    import_objs = (LB_TAG *) LB_MALLOC(num_import*sizeof(LB_TAG));
    if (!import_objs) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&proc_list);
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

   msgtag = 32767;
   LB_Comm_Create(&comm_plan, num_import, proc_list, lb->Communicator, msgtag,
                  lb->Deterministic, num_export);

  /*
   *  Allocate space for the object tags that need to be exported.  Communicate
   *  to get the list of objects to be exported.
   */

  if (*num_export > 0) {
    export_objs = (LB_TAG *) LB_MALLOC((*num_export)*sizeof(LB_TAG));
    if (!export_objs) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&proc_list);
      LB_FREE(&import_objs);
      return (LB_MEMERR);
    }
    if (!LB_Special_Malloc(lb,(void **)export_global_ids,*num_export,
                           LB_SPECIAL_MALLOC_GID)) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&proc_list);
      LB_FREE(&import_objs);
      LB_FREE(&export_objs);
      return (LB_MEMERR);
    }
    if (!LB_Special_Malloc(lb,(void **)export_local_ids,*num_export,
                           LB_SPECIAL_MALLOC_LID)) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&proc_list);
      LB_FREE(&import_objs);
      LB_FREE(&export_objs);
      LB_Special_Free(lb,(void **)export_global_ids,LB_SPECIAL_MALLOC_GID);
      return (LB_MEMERR);
    }
    if (!LB_Special_Malloc(lb,(void **)export_procs,*num_export,
                           LB_SPECIAL_MALLOC_INT)) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&proc_list);
      LB_FREE(&import_objs);
      LB_FREE(&export_objs);
      LB_Special_Free(lb,(void **)export_global_ids,LB_SPECIAL_MALLOC_GID);
      LB_Special_Free(lb,(void **)export_local_ids,LB_SPECIAL_MALLOC_LID);
      return (LB_MEMERR);
    }

  }
  else {
    export_objs = NULL;
    *export_global_ids = NULL;
    *export_local_ids = NULL;
    *export_procs = NULL;
  }

  msgtag2 = 32766;
  LB_Comm_Do(comm_plan, msgtag2, (char *) import_objs, (int) sizeof(LB_TAG), 
          (char *) export_objs);

  /*
   *  Put the exported LB_TAGs into the output format.
   */

  for (i = 0; i < *num_export; i++) {
    LB_SET_GID((*export_global_ids)[i], export_objs[i].Global_ID);
    LB_SET_LID((*export_local_ids)[i], export_objs[i].Local_ID);
    (*export_procs)[i]      = export_objs[i].Proc;
  }

  LB_FREE(&proc_list);
  LB_FREE(&import_objs);
  LB_FREE(&export_objs);
  
  LB_Comm_Destroy(&comm_plan);

  return (LB_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Help_Migrate(
  LB *lb,                    /* Load balancing structure for current balance.*/
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
int id_size;             /* size (in bytes) of LB_GID + padding for 
                            alignment                                       */
char *export_buf = NULL; /* buffer for packing export data.                 */
char *import_buf = NULL; /* buffer for receiving imported data.             */
char *tmp;               /* temporary pointer into buffers.                 */
int i;                   /* loop counter.                                   */
int tmp_import;          /* number of objects to be imported.               */
int *proc_list = NULL;   /* list of processors to which this proc exports.  */
LB_GID global_id;        /* tmp global ID for unpacking objects.            */
LB_GID *tmp_id;          /* pointer to storage for an LB_GID in comm buf    */
COMM_OBJ *comm_plan;     /* Object returned by communication routines       */
int msgtag, msgtag2;     /* Tags for communication routines                 */
int ierr = 0;

  /*
   *  Return if this processor is not in the load-balancing structure's
   *  communicator.
   */

  if (LB_PROC_NOT_IN_COMMUNICATOR(lb))
    return (LB_OK);

  if (lb->Debug_Level > 4)
    printf("ZOLTAN %d %s Entering HELP_MIGRATE %d %d\n",
            lb->Proc, yo, num_import, num_export);

  if (lb->Migrate.Get_Obj_Size == NULL) {
    fprintf(stderr, "ZOLTAN %d %s Error:  Must register an "
           "LB_OBJ_SIZE_FN_TYPE function to use the migration-help tools.\n",
           lb->Proc, yo);
    return (LB_FATAL);
  }

  if (lb->Migrate.Pack_Obj == NULL) {
    fprintf(stderr, "ZOLTAN %d %s Error:  Must register an "
           "LB_PACK_OBJ_FN_TYPE function to use the migration-help tools.\n",
           lb->Proc, yo);
    return (LB_FATAL);
  }

  if (lb->Migrate.Unpack_Obj == NULL) {
    fprintf(stderr, "ZOLTAN %d %s Error:  Must register an "
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

    if (lb->Debug_Level > 5)
      printf("ZOLTAN %d %s Done Pre-Process\n", lb->Proc, yo);
  }

  /*
   * For each object, allow space for its LB_GID and its data.
   * Zoltan will pack the LB_GIDs; the application must pack the data
   * through the pack routine.  Zoltan needs the LB_GIDs for unpacking,
   * as the order of the data received during communication is not 
   * necessarily the same order as import_global_ids[].
   */
  id_size = sizeof(LB_GID) + LB_pad_for_alignment(sizeof(LB_GID));
  size = lb->Migrate.Get_Obj_Size(lb->Migrate.Get_Obj_Size_Data, &ierr)
       + id_size;
  if (ierr) {
    fprintf(stderr, "[%d] %s: Error returned from user defined "
                    "Migrate.Get_Obj_Size function.\n", lb->Proc, yo);
    return (LB_FATAL);
  }


  if (num_export > 0) {
    export_buf = (char *) LB_MALLOC(num_export*size);
    if (!export_buf) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      return (LB_FATAL);
    }

    proc_list = (int *) LB_MALLOC(num_export*sizeof(int));
    if (!proc_list) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&export_buf);
      return (LB_FATAL);
    }


    /*
     *  Pack the proc_list (to create the map) and the objects for export.
     */
  
    tmp = export_buf;
    for (i = 0; i < num_export; i++) {
      proc_list[i] = export_procs[i];

      /* Pack the object's global ID */
      tmp_id = (LB_GID *) tmp;
      LB_SET_GID(*tmp_id, export_global_ids[i]);
    
      /* Pack the object's data */
      lb->Migrate.Pack_Obj(lb->Migrate.Pack_Obj_Data, export_global_ids[i],
                           export_local_ids[i], export_procs[i], size,
                           tmp+id_size, &ierr);
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

  msgtag = 32767;
  LB_Comm_Create(&comm_plan, num_export, proc_list, lb->Communicator, msgtag,
                 lb->Deterministic, &tmp_import);
  if (tmp_import != num_import) {
    fprintf(stderr, "%d  Error in %s:  tmp_import %d != num_import %d\n", 
            lb->Proc, yo, tmp_import, num_import);
  }

  if (num_import > 0) {
    import_buf = (char *) LB_MALLOC(num_import*size);
    if (!import_buf) {
      fprintf(stderr, "[%d] Error from %s: Insufficient memory\n",
              lb->Proc, yo);
      LB_FREE(&export_buf);
      LB_FREE(&proc_list);
      return (LB_FATAL);
    }

  }

  /*
   *  Send the export data using the communication plan.
   */

  msgtag2 = 32766;
  LB_Comm_Do(comm_plan, msgtag2, export_buf, size, import_buf);

  /*
   *  Free whatever memory we can.
   */

  LB_Comm_Destroy(&comm_plan);
  LB_FREE(&proc_list);
  LB_FREE(&export_buf);

  /*
   *  Unpack the object data.
   */

  tmp = import_buf;
  for (i = 0; i < num_import; i++) {

    /* Unpack the object's global ID */
    tmp_id = (LB_GID *) tmp;
    LB_SET_GID(global_id, *tmp_id);

    /* Unpack the object's data */
    lb->Migrate.Unpack_Obj(lb->Migrate.Unpack_Obj_Data, global_id, size,
                           tmp+id_size, &ierr);
    if (ierr) {
      fprintf(stderr, "[%d] %s: Error returned from user defined "
                      "Migrate.Unpack_Obj function.\n", lb->Proc, yo);
      return (LB_FATAL);
    }
    tmp += size;
  }

  LB_FREE(&import_buf);

  if (lb->Migrate.Post_Process != NULL) {
    lb->Migrate.Post_Process(lb->Migrate.Post_Process_Data,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, &ierr);
    if (ierr) {
      fprintf(stderr, "[%d] %s: Error returned from user defined "
                      "Migrate.Post_Process function.\n", lb->Proc, yo);
      return (LB_FATAL);
    }

    if (lb->Debug_Level > 5)
      printf("ZOLTAN %d %s Done Post-Process\n", lb->Proc, yo);
  }

  if (lb->Debug_Level > 4)
    printf("ZOLTAN %d %s Leaving HELP_MIGRATE %d %d\n",
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

  LB_FREE(import_global_ids);
  LB_FREE(import_local_ids);
  LB_FREE(import_procs);
  LB_FREE(export_global_ids);
  LB_FREE(export_local_ids);
  LB_FREE(export_procs);

  return (LB_OK);

}

/************************************************************************/
/* LB_Eval evaluates the quality of the current partitioning/balance.   */
/************************************************************************/

#define NUM_GSTATS 4 /* Number of graph statistics */

void LB_Eval (LB *lb, int print_stats, 
     int *nobj, float *obj_wgt, int *cut_wgt, int *nboundary,
     int *nadj, int *ierr)
/* 
 * Input:
 *   lb          - pointer to lb structure
 *   print_stats - if > 0, compute and print max and sum of the metrics
 *
 * Output:
 *   nobj      - number of objects (for each proc)
 *   obj_wgt   - obj_wgt[0:lb->Obj_Weight_Dim-1] are the object weights (on each proc)
 *   cut_wgt   - cut size/weight (for each proc)
 *   nboundary - number of boundary objects (for each proc)
 *   nadj      - the number of adjacent procs (for each proc)
 *   ierr      - error code
 *
 * Ouput parameters will only be returned if they are 
 * not NULL on entry (except for the error code ierr).
 */

{
  char *yo = "LB_Eval";
  char *yo2 = "ZOLTAN LB_Eval";
  int i, j, num_obj, max_edges, flag, nedges;
  int num_adj, num_boundary, cut_weight;
  int stats[3*NUM_GSTATS], *ewgts;
  int *proc, *nbors_proc;
  float *tmp_wgt, *vwgts, nproc;
  LB_LID * local_ids;
  LB_GID *global_ids, *nbors_global;
  
  /* Set default error code */
  *ierr = LB_OK;

  /* Set all pointers to NULL */
  global_ids = NULL;
  local_ids = NULL;
  tmp_wgt = NULL;
  vwgts = NULL;
  ewgts = NULL;
  nbors_global = NULL;
  nbors_proc = NULL;
  proc = NULL;

  /* First compute number of objs and object weight on each proc */
  num_obj = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, ierr);

  if (num_obj>0){

    /* Allocate space for object data */
    global_ids = (LB_GID *) LB_MALLOC(num_obj * sizeof(LB_GID));
    local_ids  = (LB_LID *) LB_MALLOC(num_obj * sizeof(LB_LID));
      
    if ((!global_ids) || (!local_ids)){
      *ierr = LB_MEMERR;
      LB_FREE(&global_ids);
      LB_FREE(&local_ids);
      return;
    }
  }
  

  /* Allocate space for weights if needed */
  if (lb->Obj_Weight_Dim>0){
    vwgts   = (float  *) LB_MALLOC(lb->Obj_Weight_Dim*num_obj * sizeof(float));
    tmp_wgt = (float *) LB_MALLOC(3*lb->Obj_Weight_Dim * sizeof(float));
    if ((num_obj && !vwgts) || (!tmp_wgt)){
      *ierr = LB_MEMERR;
      LB_FREE(&global_ids);
      LB_FREE(&local_ids);
      LB_FREE(&vwgts);
      LB_FREE(&tmp_wgt);
      return;
    }
  } 
  
  LB_Get_Obj_List(lb, global_ids, local_ids, lb->Obj_Weight_Dim, vwgts, ierr);
  if (*ierr == LB_FATAL){
    LB_FREE(&global_ids);
    LB_FREE(&local_ids);
    return;
  }
  

  /* Compute object weight sums */
  if (lb->Obj_Weight_Dim>0){
    for (j=0; j<lb->Obj_Weight_Dim; j++)
      tmp_wgt[j] = 0;
    for (i=0; i<num_obj; i++){
      for (j=0; j<lb->Obj_Weight_Dim; j++){
        tmp_wgt[j] += vwgts[i*lb->Obj_Weight_Dim+j];
      }
    }
  }


  /* Compute (weighted) edge cuts, #boundary vertices,
     and # adjacent procs if possible */

  cut_weight = 0;
  num_boundary = 0;
  num_adj = 0;

  if (lb->Get_Num_Edges != NULL) {
    /* Use the basic graph query functions */

    /* First compute max no. of edges so we can allocate the right
       amount of space */
    max_edges = 0;
    for (i=0; i< num_obj; i++){
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], ierr);
      if (*ierr){
        printf("Error in %s: Get_Num_Edges returned error code %d\n", 
          yo, *ierr);
        LB_FREE(&global_ids);
        LB_FREE(&local_ids);
        LB_FREE(&vwgts);
        LB_FREE(&tmp_wgt);
        return;
      }
      if (nedges>max_edges) max_edges = nedges;
    }

    /* Allocate edge list space */
    nbors_global = (LB_GID *)LB_MALLOC(max_edges * sizeof(LB_GID));
    nbors_proc = (int *)LB_MALLOC(max_edges * sizeof(int));
    ewgts = (int *)LB_MALLOC(lb->Comm_Weight_Dim*max_edges * sizeof(int));
    /* Allocate a proc list for computing nadjacent */
    proc = (int *)LB_MALLOC((lb->Num_Proc)* sizeof(int));

    if (max_edges && ((!nbors_global) || (!nbors_proc) || 
        (lb->Comm_Weight_Dim && (!ewgts)) || (!proc))){
      *ierr = LB_MEMERR;
      LB_FREE(&global_ids);
      LB_FREE(&local_ids);
      LB_FREE(&vwgts);
      LB_FREE(&tmp_wgt);
      LB_FREE(&nbors_global);
      LB_FREE(&nbors_proc);
      LB_FREE(&ewgts);
      LB_FREE(&proc);
      return;
    }

    for (i=0; i<lb->Num_Proc; i++)
      proc[i] = 0;

    for (i=0; i<num_obj; i++){
      flag = 0;
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], ierr);
      if (*ierr == LB_FATAL){
        LB_FREE(&global_ids);
        LB_FREE(&local_ids);
        LB_FREE(&vwgts);
        LB_FREE(&tmp_wgt);
        LB_FREE(&nbors_global);
        LB_FREE(&nbors_proc);
        LB_FREE(&ewgts);
        LB_FREE(&proc);
        return;
      }
      lb->Get_Edge_List(lb->Get_Edge_List_Data, global_ids[i], local_ids[i],
          nbors_global, nbors_proc, lb->Comm_Weight_Dim, ewgts, ierr);
      if (*ierr == LB_FATAL){
        LB_FREE(&global_ids);
        LB_FREE(&local_ids);
        LB_FREE(&vwgts);
        LB_FREE(&tmp_wgt);
        LB_FREE(&nbors_global);
        LB_FREE(&nbors_proc);
        LB_FREE(&ewgts);
        LB_FREE(&proc);
        return;
      }
      /* Check for cut edges */
      for (j=0; j<nedges; j++){
        if (nbors_proc[j] != lb->Proc){
          if (lb->Comm_Weight_Dim == 0)
            cut_weight++;
          else if (lb->Comm_Weight_Dim == 1)
            cut_weight += ewgts[j];
          else{
            fprintf(stderr, "Error in %s: Comm_Weight_Dim=%d not supported\n", 
                    yo, lb->Comm_Weight_Dim);
            *ierr = LB_WARN;
          }
          if (flag==0){
            num_boundary++;
            flag = 1;
          }
          proc[nbors_proc[j]]++;
        }
      }
    }
    /* Compute the number of adjacent procs */
    for (j=0; j<lb->Num_Proc; j++)
      if (proc[j]>0) num_adj++;
  }
  else{
    /* No graph query functions available */
  }
  
  if (print_stats){
    /* Global reduction for object weights. */
    if (lb->Obj_Weight_Dim>0){
      MPI_Allreduce(tmp_wgt, &tmp_wgt[lb->Obj_Weight_Dim], lb->Obj_Weight_Dim, MPI_FLOAT, MPI_MAX, 
                    lb->Communicator);
      MPI_Allreduce(tmp_wgt, &tmp_wgt[2*lb->Obj_Weight_Dim], lb->Obj_Weight_Dim, MPI_FLOAT, 
                    MPI_SUM, lb->Communicator);
    }
    stats[0] = num_obj;
    stats[1] = cut_weight;
    stats[2] = num_boundary;
    stats[3] = num_adj;

    /* Compute max and sum in the upper portions of the stats array. */
    MPI_Allreduce(stats, &stats[NUM_GSTATS], NUM_GSTATS, MPI_INT, MPI_MAX, 
                  lb->Communicator);
    MPI_Allreduce(stats, &stats[2*NUM_GSTATS], NUM_GSTATS, MPI_INT, MPI_SUM, 
                  lb->Communicator);

    /* Print max-sum of results */
    nproc = lb->Num_Proc; /* convert to float */
    if (lb->Proc == 0){
      printf("\n%s  Statistics for current partitioning/balance:\n", yo2);
      for (i=0; i<lb->Obj_Weight_Dim; i++)
        printf("%s  Object weight #%1d :  Max = %6.1f, Sum = %7.1f, "
          "Imbal. = %5.3f\n",
          yo2, i+1, tmp_wgt[lb->Obj_Weight_Dim+i], tmp_wgt[2*lb->Obj_Weight_Dim+i], 
          (tmp_wgt[2*lb->Obj_Weight_Dim+i] > 0
            ? tmp_wgt[lb->Obj_Weight_Dim+i]*nproc/tmp_wgt[2*lb->Obj_Weight_Dim+i]
            : 1.));
      printf("%s  No. of objects   :  Max = %6d, Sum = %7d, Imbal. = %5.3f\n",
        yo2, stats[NUM_GSTATS], stats[2*NUM_GSTATS], 
        (stats[2*NUM_GSTATS] > 0
              ? stats[NUM_GSTATS]*nproc/stats[2*NUM_GSTATS]
              : 1.));
      if (lb->Get_Num_Edges != NULL){
        printf("%s  Cut weight       :  Max = %6d, Sum = %7d, Imbal. = %5.3f\n",
          yo2, stats[NUM_GSTATS+1], stats[2*NUM_GSTATS+1], 
          (stats[2*NUM_GSTATS+1] > 0 
                ? stats[NUM_GSTATS+1]*nproc/stats[2*NUM_GSTATS+1]
                : 1.));
        printf("%s  Boundary objects :  Max = %6d, Sum = %7d, Imbal. = %5.3f\n",
          yo2, stats[NUM_GSTATS+2], stats[2*NUM_GSTATS+2], 
          (stats[2*NUM_GSTATS+2] > 0
                ? stats[NUM_GSTATS+2]*nproc/stats[2*NUM_GSTATS+2]
                : 1.));
        printf("%s  Adjacent procs   :  Max = %6d, Sum = %7d, Imbal. = %5.3f\n",
          yo2, stats[NUM_GSTATS+3], stats[2*NUM_GSTATS+3], 
          (stats[2*NUM_GSTATS+3] > 0
                ? stats[NUM_GSTATS+3]*nproc/stats[2*NUM_GSTATS+3]
                : 1.));
      }
      printf("\n");
    }
  }

  /* Copy results to output parameters if desired */
  if (nobj) *nobj = num_obj;
  if (nadj) *nadj = num_adj;
  if (nboundary) *nboundary = num_boundary;
  if (cut_wgt) *cut_wgt = cut_weight;
  if (obj_wgt){
    for (i=0; i<lb->Obj_Weight_Dim; i++) 
      obj_wgt[i] = tmp_wgt[i];
  }

  /* Free data */
  LB_FREE(&global_ids);
  LB_FREE(&local_ids);
  LB_FREE(&tmp_wgt);
  LB_FREE(&vwgts);
  LB_FREE(&ewgts);
  LB_FREE(&nbors_global);
  LB_FREE(&nbors_proc);
  LB_FREE(&proc);
}
