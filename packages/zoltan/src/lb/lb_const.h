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
static char *cvs_lbconsth_id = "$Id$";
#endif


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "par_const.h"
#include "lbi_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Type definitions.
 */

typedef struct LB_Struct LB;
typedef struct LB_Comm_Struct LB_COMM;
typedef struct LB_Tag_Struct LB_TAG;

typedef void LB_FN(LB *, int*, int*, int *, LB_TAG **);
typedef void LB_COMM_BUILD_REQUEST_PROCLIST_FN_TYPE(LB *, int n_cells_orig, int *n_requests);
typedef void LB_COMM_BUILD_SEND_REQUEST_LIST_FN_TYPE(LB *, int n_cells_orig);
typedef int  LB_COMM_OBJ_DATA_SIZE_FN_TYPE(int object_type);
typedef void LB_COMM_MIGRATE_OBJ_DATA_FN_TYPE(LB_ID object, char *start_pos_in_buffer);

/*
 *  Define the possible load balancing methods allowed.
 */

typedef enum LB_Method {
  RCB = 0,
  WHEAT,
  LB_MAX_METHODS                  /*  This entry should always be last.      */
} LB_METHOD;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/



/* 
 *  Define a communication structure for load balancing results.  This
 *  structure will be used by Steve Plimpton's and Bruce Hendrickson's
 *  communication library.
 */

struct LB_Comm_Struct {
  /*  
   *  Pointers to routines that depend on the LB Data_Structure field, and,
   *  thus, on the balancing method used.
   */
  LB_COMM_BUILD_REQUEST_PROCLIST_FN_TYPE *Build_Request_Proclist;
                                       /* Routine that build a list of procs
                                          from which object data are needed.
                                          There is one entry in the list for
                                          each remote object.  The value of
                                          the entry is the processor number
                                          of the processor owning the object
                                          when the load balancer was invoked.*/
  LB_COMM_BUILD_SEND_REQUEST_LIST_FN_TYPE *Build_Send_Request_List;
                                       /* Routine that build a list of objects
                                          needed from other processors. 
                                          Each entry contains the tracking 
                                          data built before the load balancing
                                          was invoked.                       */

  /*
   *  Pointers to routines that depend on the application.
   */

  LB_COMM_OBJ_DATA_SIZE_FN_TYPE *Get_Obj_Data_Size;
                                       /* Function that returns the size of
                                          contiguous memory needed to store
                                          the data for a single object for
                                          migration.                         */
  LB_COMM_MIGRATE_OBJ_DATA_FN_TYPE *Pack_Obj_Data;
                                       /* Routine that packs object data for
                                          a given object into contiguous 
                                          memory for migration.              */
  LB_COMM_MIGRATE_OBJ_DATA_FN_TYPE *Unpack_Obj_Data;
                                       /* Routine that unpacks object data for
                                          a given object from contiguous 
                                          memory after migration.            */
                                        
  /*
   *  Pointers to temporary storage for communication.
   */

  int *Proc_List;                      /* Array of processor numbers for
                                          requesting and sending objects.
                                          There is one entry per requested
                                          or sent object; its value is the
                                          processor number to which the 
                                          request or object is sent.         */
  struct Request_Struct *Send_Request; /* Array of requests for objects on
                                          other processors.
                                          There is one entry per requested 
                                          object; each entry contains the
                                          tracking data giving the object's
                                          original location.                 */
  struct Request_Struct *Recv_Request; /* Array of requests for objects
                                          needed by other processors.  
                                          There is one entry per requested 
                                          object; each entry contains the
                                          tracking data giving the object's
                                          original location.                 */
  char *Send_Data;                     /* Buffer containing object data to 
                                          be sent to other processors.       */
  char *Recv_Data;                     /* Buffer containing object data to 
                                          be received from other processors. */

};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


/*
 *  Define a load balancing object.  It will contain pointers to the
 *  appropriate functions for interfacing with applications and 
 *  pointers to the data structure used for load balancing.
 */

struct LB_Struct {
  LB_METHOD Method;               /*  Method to be used for load balancing.  */
  LB_FN *LB_Fn;                   /*  Pointer to the function that performs
                                      the load balancing; this ptr is set
                                      based on the method used.              */
  double *Params;                 /*  Array of parameters passed to the 
                                      load balancing function.               */
  double Tolerance;               /*  Tolerance to which to load balance;
                                      tolerance = 0.9 implies 10% imbalance
                                      is acceptable.                         */
  int Object_Type;                /*  The application-specified object type
                                      for objects being balanced.  The
                                      application can use this value to 
                                      distinguish which objects (e.g.,
                                      elements or surfaces) are being used
                                      in this load-balancing object.  This
                                      value is not used specifically by the
                                      load balancer.                         */
  int Help_Migrate;               /*  Flag indicating whether the load
                                      balancer should help the application
                                      migrate data.  Some applications may
                                      prefer to do it themselves.            */
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
  LB_COMM LB_Comm;                             /* Communication struct for
                                                  load balancing results.    */
};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Structure used to describe the results of the new decomposition.  The
 *  load-balancing routines (lb_rcb, etc.) return an array of LB_Tag_Structs
 *  with one entry for each non-local (i.e., imported) object in the new
 *  new decomposition for a given processor.
 *  This structure is the minimum structure required in the load-balancing
 *  data structures.
 */

struct LB_Tag_Struct {
  LB_ID Global_ID;         /* The global ID of the related object.           */
  LB_ID Local_ID;          /* The local ID of the related object.            */
  int Proc;                /* The original processor of the object.  Also 
                              used for target processor in inverse comm. map */
};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */

extern LB_FN lb_rcb;
extern LB_FN lb_wheat;
