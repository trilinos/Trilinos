#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "all_allo_const.h"
#include "par_const.h"
#include "par.h"
#include "lb_const.h"
#include "lb.h"

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
    MPI_Init(argc, argv);
  }

  /*
   *  Set global processor values for the load balacing tools.
   */

  MPI_Comm_rank(MPI_COMM_WORLD, &Proc);
  MPI_Comm_size(MPI_COMM_WORLD, &Num_Proc);

#endif  /* LB_MPI */
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void *LB_Create_LB_Object()
{
/*
 *  Function to create a load balancing object.  May want more than one
 *  object if using different decompositions with different techniques.
 *  This function allocates and initializes the object.
 *  Output:
 *    void *             --  Pointer to a LB object.
 *
 */

char *yo = "LB_Create_LB_Object";
LB *lb;

  /*
   * Allocate storage for the load-balancing object.
   */

  lb = (LB *) smalloc(sizeof(LB));

  /*
   *  Set defaults for fields of lb:
   */

  lb->Method = RCB;    
  lb->Tolerance = 0.9;
  lb->Data_Structure = NULL;

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
  
  /*
   *  Cast return value.
   */

  return((void *) lb);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Set_LB_Fn(void *lbv, LB_FN_TYPE fn_type, void *fn())
{
/*
 *  Function to initialize a given LB interface function.
 *  Input:
 *    lbv               --  Pointer to a LB object.
 *    fn_type           --  Enum type indicating the function to be set.
 *    fn                --  Pointer to the function to be used in the
 *                          assignment.
 *  Output:
 *    lbv               --  Appropriate field set to value in void *().
 */

char *yo = "LB_Set_LB_Fn";
LB *lb = (LB *) lbv;

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

void LB_Set_LB_Method(void *lbv, char *method_name, int *params)
{
/*
 *  Function to set the load balancing method to be used.
 *  Input:
 *    lbv                --  The load balancing object to which this method
 *                           applies.
 *    method_name        --  String specifying the desired method.
 *    params             --  Params needed by desired method.  (This field
 *                           will be better defined later.)
 *  Output:
 *    lbf*               --  Appropriate fields set to designated values.
 */

char *yo = "LB_Set_LB_Method";
LB *lb = (LB *) lbv;
int i;

  /*
   *  Compare method_name string with standard strings for methods.
   *  If a match is found, set lb->Method.
   */

  for (i = 0; i < LB_MAX_METHODS; i++) {
    if (strcasecmp(method_name, LB_Method_Strings[i]) == 0) {
      lb->Method = i;
      if (Proc == 0) {
        printf("LB:  Load balancing method = %d (%s)\n", i, method_name);
      }
      break;
    }
  }

  if (i == LB_MAX_METHODS) {
    fprintf(stderr, "Error from %s:  Invalid LB method specified:  %s\n", 
            yo, method_name);
    exit(-1);
  }

  /*
   *  For now, don't know how params should be defined.  It will probably
   *  depend on the method.  We'll ignore it for now.
   */

  if (params != NULL) {
    fprintf(stderr, "Warning from %s:  Params array ignored.\n", yo);
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void LB_Set_LB_Tolerance(void *lbv, double tolerance)
{
/*
 *  Function to set the tolerance to which the system must be load balanced.
 *  For example, if the tolerance is set to 0.9, 10% load imbalance between
 *  the most heavily loaded processor and the average load will be accepted
 *  as balanced.
 *  Input:
 *    lbv                --  The load balancing object to which this tolerance
 *                           applies.
 *    tolerance          --  The tolerance desired.
 *  Output:
 *    lbv                --  Tolerance field set to appropriate value.
 */

char *yo = "LB_Set_LB_Tolerance";
LB *lb = (LB *) lbv;

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

  if (Proc == 0) {
    printf("LB:  Load balancing tolerance = %f\n", tolerance);
  }
}
