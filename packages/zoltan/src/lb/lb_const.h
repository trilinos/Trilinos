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
#ifndef __LB_CONST_H
#define __LB_CONST_H

#ifndef lint
static char *cvs_lbconsth_id = "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

#include "par_const.h"
#include "lbi_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Type definitions.
 */

/*
 *  Value used to initialize the parameters when the
 *  load-balancer allocates the parameters array.
 */

#define LB_PARAMS_INIT_VALUE -1.

#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */

typedef struct LB_Struct LB;
typedef struct LB_Comm_Struct LB_COMM;

typedef void LB_FN(LB *, int*, int*, int *, LB_TAG **);

/*
 *  Define the possible load balancing methods allowed.
 */

typedef enum LB_Method {
  NONE = -1,
  RCB,
  OCTPART,
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

struct LB_Migrate_Struct {
  int Help_Migrate;                   /*  Flag indicating whether the load
                                          balancer should help the application
                                          migrate data.  Some applications may
                                          prefer to do it themselves.        */
  /*
   *  Pointers to routines that depend on the application.
   */

  LB_PRE_MIGRATE_FN *Pre_Process;      /* Function that performs application
                                          specific pre-processing.  Optional
                                          for help-migration.                */
  LB_OBJECT_SIZE_FN *Get_Obj_Data_Size;/* Function that returns the size of
                                          contiguous memory needed to store
                                          the data for a single object for
                                          migration.                         */
  LB_PACK_OBJECT_FN *Pack_Obj_Data;    /* Routine that packs object data for
                                          a given object into contiguous 
                                          memory for migration.              */
  LB_UNPACK_OBJECT_FN *Unpack_Obj_Data;
                                       /* Routine that unpacks object data for
                                          a given object from contiguous 
                                          memory after migration.            */
};

typedef struct LB_Migrate_Struct LB_MIGRATE;

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
  LB_MIGRATE Migrate;                          /* Struct with info for helping
                                                  with migration.            */
};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */

extern LB_FN lb_rcb;
extern LB_FN lb_wheat;
extern LB_FN lb_oct_init;

#endif
