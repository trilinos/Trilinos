
#include "lbi_const.h"

/*
 *  Define the possible load balancing methods allowed.
 *  The order of this type MUST be the same as the LB_Method_Strings array
 *  defined in lb.h.
 */

typedef enum LB_Method {
  RCB = 0,
  WHEAT,
  LB_MAX_METHODS                  /*  This entry should always be last.      */
} LB_METHOD;

extern char *LB_Method_Strings[];

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Define a load balancing object.  It will contain pointers to the
 *  appropriate functions for interfacing with applications and 
 *  pointers to the data structure used for load balancing.
 */

typedef struct LB_Struct {
  LB_METHOD Method;               /*  Method to be used for load balancing.  */
  double Tolerance;               /*  Tolerance to which to load balance;
                                      tolerance = 0.9 implies 10% imbalance
                                      is acceptable.                         */
  void *Data_Structure;           /*  Data structure used by the load 
                                      balancer; cast by the method routines
                                      to the appropriate data type.          */
  LB_OBJECT_WEIGHT_FN *Get_Obj_Weight;         /* Fn ptr to get an object's
                                                  weight.                    */
  LB_NUM_EDGES_FN *Get_Num_Edges;              /* Fn ptr to get an object's
                                                  number of edges.           */
  LB_EDGE_LIST_FN *Get_Edge_List;              /* Fn ptr to get an object's
                                                  edge list.                 */
  LB_NUM_GEOM_FN *Get_Num_Geom;                /* Fn ptr to get an object's
                                                  number of geometry values. */
  LB_GEOM_FN *Get_Obj_Geom;                    /* Fn ptr to get an object's
                                                  geometry values.           */
  LB_NUM_OBJ_FN *Get_Num_Local_Obj;            /* Fn ptr to get a proc's  
                                                  number of local objects.   */
  LB_GET_LOCAL_OBJECTS_FN *Get_All_Local_Objs; /* Fn ptr to get all local
                                                  objects on a proc.         */
  LB_NEXT_OBJ_FN *Get_Next_Local_Obj;          /* Fn ptr to get the next   
                                                  local obj on a proc.       */
  LB_NUM_BORDER_OBJ_FN *Get_Num_Border_Obj;    /* Fn ptr to get a proc's 
                                                  number of border objs wrt
                                                  a given processor.         */
  LB_BORDER_OBJ_FN *Get_All_Border_Objs;       /* Fn ptr to get all objects
                                                  sharing a border with a
                                                  given processor.           */
  LB_NEXT_BORDER_OBJ_FN *Get_Next_Border_Obj;  /* Fn ptr to get the next 
                                                  object sharing a border 
                                                  with a given processor.    */
} LB;
