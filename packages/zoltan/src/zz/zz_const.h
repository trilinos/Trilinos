// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __ZOLTAN_CONST_H
#define __ZOLTAN_CONST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>

#include "zoltan_util.h"
#include "zoltan.h"
#include "lb_const.h"
#include "order_const.h"
#include "zz_id_const.h"
#include "par_const.h"
#include "third_library_const.h"
#include "zoltan_timer.h"

#ifdef _MSC_VER
#define __func__ __FUNCTION__
#endif

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Type definitions and forward declarations.
 */

/*
 * Strings to define the library name, and the version number
 * so that it is easier to keep track of what code each user
 * has.
 */
#define UTIL_NAME "zoltan"

/*
 * Type used to store linked list of new values for parameters.
 */
struct Param_List;

#define MIN(A,B)                (((A) < (B)) ? (A) : (B))
#define MAX(A,B)                (((A) > (B)) ? (A) : (B))

#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */

/*
 *  Define the debug levels allowed.
 *    ZOLTAN_DEBUG_NONE = 0           - quiet mode; no debugging information.
 *    ZOLTAN_DEBUG_PARAMS = 1         - print values of all parameters used 
 *                                           plus backtrace on Zoltan error.
 *    ZOLTAN_DEBUG_ZTIME = 2          - print Zoltan timing information.
 *    ZOLTAN_DEBUG_ATIME = 3          - print algorithm's timing info, if the
 *                                  algorithm supports this level.
 *    ZOLTAN_DEBUG_TRACE_ZERO = 5     - print trace info on processor 0 only.
 *    ZOLTAN_DEBUG_TRACE_ALL = 6      - print trace info on all processors.
 *    ZOLTAN_DEBUG_TRACE_DETAIL = 7   - print detailed trace info on all processors.
 *    ZOLTAN_DEBUG_LIST = 8           - print lists of objects to be imported 
 *                                  and exported.
 *    ZOLTAN_DEBUG_ALL = 10           - print all debug information available.
 */
#define ZOLTAN_DEBUG_NONE 0     
#define ZOLTAN_DEBUG_PARAMS 1
#define ZOLTAN_DEBUG_ZTIME 2
#define ZOLTAN_DEBUG_ATIME 3
#define ZOLTAN_DEBUG_TRACE_SINGLE 5
#define ZOLTAN_DEBUG_TRACE_ALL 6
#define ZOLTAN_DEBUG_TRACE_DETAIL 7 
#define ZOLTAN_DEBUG_LIST 8
#define ZOLTAN_DEBUG_ALL 10

/*
 ******************************************************
 * Define default values for key parameters.
 * These are used in both lb.c and key_params.c.
 ******************************************************
 */
#define ZOLTAN_DEBUG_LEVEL_DEF    ZOLTAN_DEBUG_PARAMS
#define ZOLTAN_DEBUG_PROC_DEF     0
#define ZOLTAN_OBJ_WEIGHT_DEF     0
#define ZOLTAN_EDGE_WEIGHT_DEF    0
#define ZOLTAN_DETERMINISTIC_DEF  TRUE
#define ZOLTAN_NUM_ID_ENTRIES_DEF 1
#define ZOLTAN_TIMER_DEF          ZOLTAN_TIME_WALL
#define ZOLTAN_TFLOPS_SPECIAL_DEF FALSE

/*****************************************************************************/
/*****************************************************************************/
/* Symbols used by various Zoltan routines. */
#define ZOLTAN_LOCAL  1
#define ZOLTAN_GLOBAL 2

/*****************************************************************************/
/*****************************************************************************/
/* Structs */
struct Zoltan_Migrate_Struct;
struct Zoltan_LB_Struct;
struct Zoltan_Order_Struct;
struct Zoltan_TPL_Order_Struct;

/*****************************************************************************/
/*****************************************************************************/
 
/**************************************************************************/
/* The data structure below for hetero. machines is not being used yet!   */
/* The structure will almost certainly change in the next release.        */
/**************************************************************************/

typedef struct {
  int nnodes;           /* the number of subnodes */
  int ntypes;           /* the number of different types of subnodes */
  int *type;            /* type[i] is the `node type pointer' of subnode i */
                        /* if (ntypes == 1)        */
                           /* specify only type[0] */
                        /* else */
                           /* specify type[0] ... type[nnodes-1] */

  int top_id;          /* See `topology types' defined below */

  /************************************************/
  /* specify if (nnodes == 1)                     */
  /************************************************/
  int power;             /* if (nnodes == 1) specify power of the processor */
  int memory;            /* if (nnodes == 1) specify memory of the processor */
  /************************************************/

  /*****************************/
  /* specify if (top_id 0 and 1) */
  /*****************************/
  int bandwidth;         /* specify the bandwidth of the topology */
  int latency;           /* specify the latency of the topology */
  /*****************************/

  /*****************************/
  /*  specify if (top_id == 1) */
  /*****************************/
  int ndims;           /* Number of dimensions of the mesh */
  int cart_dim[3];     /* Number of nodes in each dimension */
  int wrap_around[3];  /* Wrap around in each dimension? (y/n) */
  /******************************/

  /*****************************/
  /*  specify if (top_id == 2) */
  /*****************************/
  int *xadj;           /* pointer to link list */
  int *adjncy;         /* link list */
  int *adj_band;       /* bandwidth of each link */
  int *adj_lat;        /* latency of each link */
  /******************************/
} MachineType;

  /*******************************/
  /* `topology types'
  top_id 0   =>   flat network or SMP
  top_id 1   =>   mesh
  top_id 2   =>   hypercube
  top_id 3   =>   user-specified */
  /*******************************/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Define a Zoltan structure.  It will contain pointers to the
 *  appropriate functions for interfacing with applications and 
 *  pointers to the data structure used for load balancing.
 */

struct Zoltan_Struct {
  MPI_Comm Communicator;          /*  The MPI Communicator.                  */
  int Proc;                       /*  The processor's ID within the MPI
                                      Communicator.                          */
  int Num_Proc;                   /*  The number of processors in the MPI
                                      Communicator.                          */
  int Num_GID;                    /*  The number of entries in Global IDs.   */
  int Num_LID;                    /*  The number of entries in Local IDs.    */
  int Debug_Level;                /*  Debug level for this instance of
                                      load balancing.                        */
  int Debug_Proc;                 /*  Print from this processor any debugging 
                                      info that is printed from only one 
                                      processor.                             */
  int Fortran;                    /*  1 if created from Fortran, 0 otherwise */
  int Tflops_Special;             /*  Flag to indicate if we should use some
                                      MPI constructs (0) or not (1) on tflops*/
  unsigned int Seed;              /*  Zoltan_Rand seed: default or provided 
                                      by user.     */
  struct Param_List *Params;      /*  List of parameter names & new vals     */
  int Deterministic;              /*  Flag indicating whether algorithms used
                                      should be forced to be deterministic.
                                      Default = TRUE.                        */
  int Obj_Weight_Dim;             /*  Dimension of the object weights, 
                                      usually 0 (no weights) or 1            */
  int Edge_Weight_Dim;            /*  Dimension of the edge weights, 
                                      usually 0 (no weights) or 1            */
  int Timer;                      /*  Timer type that is currently active */
  struct Zoltan_Timer *ZTime;     /*  Timer structure for persistent timing. */
  /***************************************************************************/
  ZOLTAN_PART_MULTI_FN *Get_Part_Multi;/* Fn ptr to get objects'
                                          part assignments.     */
  ZOLTAN_PART_MULTI_FORT_FN *Get_Part_Multi_Fort;
                                       /* Fortran version            */
  void *Get_Part_Multi_Data;           /* Ptr to user defined data to be 
                                          passed to Get_Part_Multi()    */
  ZOLTAN_PART_FN *Get_Part;            /* Fn ptr to get an object's
                                          part assignment.      */
  ZOLTAN_PART_FORT_FN *Get_Part_Fort;  /* Fortran version            */
  void *Get_Part_Data;                 /* Ptr to user defined data
                                          to be passed to Get_Part()    */
  /***************************************************************************/
  ZOLTAN_NUM_EDGES_FN *Get_Num_Edges;  /* Fn ptr to get an object's
                                          number of edges.           */
  ZOLTAN_NUM_EDGES_FORT_FN *Get_Num_Edges_Fort;
                                       /* Fortran version            */
  void *Get_Num_Edges_Data;            /* Ptr to user defined data
                                          to be passed to Get_Num_Edges()    */
  ZOLTAN_NUM_EDGES_MULTI_FN *Get_Num_Edges_Multi;  
                                       /* Fn ptr to get multiple objects'
                                          number of edges.           */
  ZOLTAN_NUM_EDGES_MULTI_FORT_FN *Get_Num_Edges_Multi_Fort;
                                       /* Fortran version            */
  void *Get_Num_Edges_Multi_Data;      /* Ptr to user defined data
                                          to be passed to Get_Num_Edges_Multi */
  /***************************************************************************/
  ZOLTAN_EDGE_LIST_FN *Get_Edge_List;  /* Fn ptr to get an object's edge list.*/
  ZOLTAN_EDGE_LIST_FORT_FN *Get_Edge_List_Fort;
                                       /* Fortran version            */
  void *Get_Edge_List_Data;            /* Ptr to user defined data
                                          to be passed to Get_Edge_List()    */
  ZOLTAN_EDGE_LIST_MULTI_FN *Get_Edge_List_Multi;  
                                       /* Fn ptr to get an object's edge list.*/
  ZOLTAN_EDGE_LIST_MULTI_FORT_FN *Get_Edge_List_Multi_Fort;
                                       /* Fortran version            */
  void *Get_Edge_List_Multi_Data;            /* Ptr to user defined data
                                          to be passed to Get_Edge_List()    */
  /***************************************************************************/
  ZOLTAN_NUM_GEOM_FN *Get_Num_Geom;    /* Fn ptr to get an object's
                                          number of geometry values. */
  ZOLTAN_NUM_GEOM_FORT_FN *Get_Num_Geom_Fort;  
                                       /* Fortran version            */
  void *Get_Num_Geom_Data;             /* Ptr to user defined data
                                          to be passed to Get_Num_Geom()     */
  /***************************************************************************/
  ZOLTAN_GEOM_MULTI_FN *Get_Geom_Multi;        
                                       /* Fn ptr to get all objects'
                                          geometry values.           */
  ZOLTAN_GEOM_MULTI_FORT_FN *Get_Geom_Multi_Fort; 
                                       /* Fortran version         */
  void *Get_Geom_Multi_Data;           /* Ptr to user defined data
                                          to be passed to Get_Geom_Multi()   */
  ZOLTAN_GEOM_FN *Get_Geom;            /* Fn ptr to get an object's
                                          geometry values.           */
  ZOLTAN_GEOM_FORT_FN *Get_Geom_Fort;  /* Fortran version            */
  void *Get_Geom_Data;                 /* Ptr to user defined data
                                          to be passed to Get_Geom()         */
  /***************************************************************************/
  ZOLTAN_NUM_OBJ_FN *Get_Num_Obj;      /* Fn ptr to get a proc's  
                                          number of local objects.   */
  ZOLTAN_NUM_OBJ_FORT_FN *Get_Num_Obj_Fort;    
                                       /* Fortran version            */
  void *Get_Num_Obj_Data;              /* Ptr to user defined data
                                          to be passed to Get_Num_Obj()      */
  /***************************************************************************/
  ZOLTAN_OBJ_LIST_FN *Get_Obj_List;    /* Fn ptr to get all local
                                          objects on a proc.         */
  ZOLTAN_OBJ_LIST_FORT_FN *Get_Obj_List_Fort;  
                                       /* Fortran version            */
  void *Get_Obj_List_Data;             /* Ptr to user defined data
                                          to be passed to Get_Obj_List()     */
  ZOLTAN_FIRST_OBJ_FN *Get_First_Obj;  /* Fn ptr to get the first   
                                          local obj on a proc.       */
  ZOLTAN_FIRST_OBJ_FORT_FN *Get_First_Obj_Fort;
                                       /* Fortran version            */
  void *Get_First_Obj_Data;            /* Ptr to user defined data
                                          to be passed to Get_First_Obj()    */
  ZOLTAN_NEXT_OBJ_FN *Get_Next_Obj;    /* Fn ptr to get the next   
                                          local obj on a proc.       */
  ZOLTAN_NEXT_OBJ_FORT_FN *Get_Next_Obj_Fort;  
                                       /* Fortran version            */
  void *Get_Next_Obj_Data;             /* Ptr to user defined data
                                          to be passed to Get_Next_Obj()     */
  /***************************************************************************/
  ZOLTAN_NUM_BORDER_OBJ_FN *Get_Num_Border_Obj;
                                       /* Fn ptr to get a proc's 
                                          number of border objs wrt
                                          a given processor.         */
  ZOLTAN_NUM_BORDER_OBJ_FORT_FN *Get_Num_Border_Obj_Fort; 
                                       /* Fortran version     */ 
  void *Get_Num_Border_Obj_Data;       /* Ptr to user defined data
                                          to be passed to
                                          Get_Num_Border_Obj()       */
  ZOLTAN_BORDER_OBJ_LIST_FN *Get_Border_Obj_List;  
                                       /* Fn ptr to get all objects
                                          sharing a border with a
                                          given processor.           */
  ZOLTAN_BORDER_OBJ_LIST_FORT_FN *Get_Border_Obj_List_Fort; 
                                       /* Fortran version   */
  void *Get_Border_Obj_List_Data;      /* Ptr to user defined data
                                          to be passed to
                                          Get_Border_Obj_List()      */
  ZOLTAN_FIRST_BORDER_OBJ_FN *Get_First_Border_Obj;
                                       /* Fn ptr to get the first 
                                          object sharing a border 
                                          with a given processor.    */
  ZOLTAN_FIRST_BORDER_OBJ_FORT_FN *Get_First_Border_Obj_Fort; 
                                       /* Fortran version */
  void *Get_First_Border_Obj_Data;     /* Ptr to user defined data
                                          to be passed to
                                          Get_First_Border_Obj()     */
  ZOLTAN_NEXT_BORDER_OBJ_FN *Get_Next_Border_Obj;  
                                       /* Fn ptr to get the next 
                                          object sharing a border 
                                          with a given processor.    */
  ZOLTAN_NEXT_BORDER_OBJ_FORT_FN *Get_Next_Border_Obj_Fort; 
                                       /* Fortran version   */
  void *Get_Next_Border_Obj_Data;      /* Ptr to user defined data
                                          to be passed to
                                          Get_Next_Border_Obj()      */
  /***************************************************************************/
  ZOLTAN_NUM_COARSE_OBJ_FN *Get_Num_Coarse_Obj;
                                       /* Fn ptr to get the number of
                                          elements in the coarse grid*/
  ZOLTAN_NUM_COARSE_OBJ_FORT_FN *Get_Num_Coarse_Obj_Fort; 
                                       /* Fortran version     */
  void *Get_Num_Coarse_Obj_Data;       /* Ptr to user defined data
                                          to be passed to
                                          Get_Num_Coarse_Obj()       */
  /***************************************************************************/
  ZOLTAN_COARSE_OBJ_LIST_FN *Get_Coarse_Obj_List;  
                                       /* Fn ptr to get all
                                          elements in the coarse grid*/
  ZOLTAN_COARSE_OBJ_LIST_FORT_FN *Get_Coarse_Obj_List_Fort; 
                                       /* Fortran version   */
  void *Get_Coarse_Obj_List_Data;      /* Ptr to user defined data
                                          to be passed to
                                          Get_Coarse_Obj_List()      */
  ZOLTAN_FIRST_COARSE_OBJ_FN *Get_First_Coarse_Obj;
                                       /* Fn ptr to get the first coarse
                                          obj on a proc.             */
  ZOLTAN_FIRST_COARSE_OBJ_FORT_FN *Get_First_Coarse_Obj_Fort; 
                                       /* Fortran version */
  void *Get_First_Coarse_Obj_Data;     /* Ptr to user defined data
                                          to be passed to
                                          Get_First_Coarse_Obj()     */
  ZOLTAN_NEXT_COARSE_OBJ_FN *Get_Next_Coarse_Obj;  
                                       /* Fn ptr to get the next coarse
                                          obj on a proc.             */
  ZOLTAN_NEXT_COARSE_OBJ_FORT_FN *Get_Next_Coarse_Obj_Fort; 
                                       /* Fortran version   */
  void *Get_Next_Coarse_Obj_Data;      /* Ptr to user defined data
                                          to be passed to
                                          Get_Next_Coarse_Obj()      */
  /***************************************************************************/
  ZOLTAN_NUM_CHILD_FN *Get_Num_Child;  /* Fn ptr to get the number of
                                          children of an element     */
  ZOLTAN_NUM_CHILD_FORT_FN *Get_Num_Child_Fort;
                                       /* Fortran version            */
  void *Get_Num_Child_Data;            /* Ptr to user defined data
                                          to be passed to
                                          Get_Num_Child()            */
  /***************************************************************************/
  ZOLTAN_CHILD_LIST_FN *Get_Child_List;        
                                       /* Fn ptr to get all
                                          children of an element     */
  ZOLTAN_CHILD_LIST_FORT_FN *Get_Child_List_Fort;  
                                       /* Fortran version            */
  void *Get_Child_List_Data;           /* Ptr to user defined data
                                          to be passed to
                                          Get_Child_List()           */
  /***************************************************************************/
  ZOLTAN_CHILD_WEIGHT_FN *Get_Child_Weight;    
                                       /* Fn ptr to get the weight
                                          of an element              */
  ZOLTAN_CHILD_WEIGHT_FORT_FN *Get_Child_Weight_Fort; 
                                       /* Fortran version         */
  void *Get_Child_Weight_Data;         /* Ptr to user defined data
                                          to be passed to
                                          Get_Child_Weight()         */
  /***************************************************************************/
  ZOLTAN_HG_SIZE_CS_FN *Get_HG_Size_CS;    
                                       /* Fn ptr to get size and format of
                                          hypergraph compressed storage.  */
  ZOLTAN_HG_SIZE_CS_FORT_FN *Get_HG_Size_CS_Fort;
                                       /* Fortran version            */
  void *Get_HG_Size_CS_Data;         /* Ptr to user defined data
                                        to be passed to Get_HG_Size_CS() */
  /***************************************************************************/
  ZOLTAN_HG_CS_FN *Get_HG_CS;    
                                       /* Fn ptr to get hypergraph pins
                                          in a compressed storage format.  */
  ZOLTAN_HG_CS_FORT_FN *Get_HG_CS_Fort;
                                       /* Fortran version            */
  void *Get_HG_CS_Data;                /* Ptr to user defined data
                                        to be passed to Get_HG_CS() */
  /***************************************************************************/
  ZOLTAN_HG_SIZE_EDGE_WTS_FN *Get_HG_Size_Edge_Wts;    
                                       /* Fn ptr to get size of hypergraph
                                          edge weights to be returned.  */
  ZOLTAN_HG_SIZE_EDGE_WTS_FORT_FN *Get_HG_Size_Edge_Wts_Fort;
                                       /* Fortran version            */
  void *Get_HG_Size_Edge_Wts_Data;     /* Ptr to user defined data to be
                                          passed to Get_HG_Size_Edge_Wts() */
  /***************************************************************************/
  ZOLTAN_HG_EDGE_WTS_FN *Get_HG_Edge_Wts;    
                                       /* Fn ptr to get hyperedge weights */
  ZOLTAN_HG_EDGE_WTS_FORT_FN *Get_HG_Edge_Wts_Fort;
                                       /* Fortran version            */
  void *Get_HG_Edge_Wts_Data;                /* Ptr to user defined data
                                        to be passed to Get_HG_Edge_Wts() */
  /***************************************************************************/
  ZOLTAN_NUM_FIXED_OBJ_FN *Get_Num_Fixed_Obj;
                                            /* Fn ptr to get a processor's
                                               number of fixed objects.      */
  ZOLTAN_NUM_FIXED_OBJ_FORT_FN *Get_Num_Fixed_Obj_Fort;  
                                            /* Fortran version               */
  void *Get_Num_Fixed_Obj_Data;             /* Ptr to user defined data to be
                                               passed to Get_Num_Fixed_Obj() */
  /***************************************************************************/
  ZOLTAN_FIXED_OBJ_LIST_FN *Get_Fixed_Obj_List;
                                            /* Fn ptr to get a processor's
                                               number of fixed objects.      */
  ZOLTAN_FIXED_OBJ_LIST_FORT_FN *Get_Fixed_Obj_List_Fort;  
                                            /* Fortran version               */
  void *Get_Fixed_Obj_List_Data;            /* Ptr to user defined data to be
                                               passed to Get_Fixed_Obj_List()*/
  /***************************************************************************/
  ZOLTAN_HIER_NUM_LEVELS_FN *Get_Hier_Num_Levels;
                                       /* Function that returns the number
                                          of levels for hierarchical 
                                          partitioning */
  ZOLTAN_HIER_NUM_LEVELS_FORT_FN *Get_Hier_Num_Levels_Fort;
                                       /* Fortran version             */
  void *Get_Hier_Num_Levels_Data;      /* Ptr to user defined data to be passed
                                          to Get_Hier_Num_Levels() */
  /***************************************************************************/
  ZOLTAN_HIER_PART_FN *Get_Hier_Part;  /* Function that returns the part
                                          for process at a given level in
                                          hierarchical partitioning */
  ZOLTAN_HIER_PART_FORT_FN *Get_Hier_Part_Fort;
                                       /* Fortran version             */
  void *Get_Hier_Part_Data;            /* Ptr to user defined data to be passed
                                          to Get_Hier_Part() */
  /***************************************************************************/
  ZOLTAN_HIER_METHOD_FN *Get_Hier_Method;
                                       /* Function that allows app to set the
                                          LB method and params for process 
                                          at a given level in
                                          hierarchical partitioning */
  ZOLTAN_HIER_METHOD_FORT_FN *Get_Hier_Method_Fort;
                                       /* Fortran version             */
  void *Get_Hier_Method_Data;          /* Ptr to user defined data to be passed
                                          to Get_Hier_Method() */
  /***************************************************************************/
  ZOLTAN_OBJ_SIZE_FN *Get_Obj_Size;    /* Function that returns the size of
                                          contiguous memory needed to store
                                          the data for a single object for
                                          migration.                         */
  ZOLTAN_OBJ_SIZE_FORT_FN *Get_Obj_Size_Fort;
                                       /* Fortran version                    */
  void *Get_Obj_Size_Data;             /* Ptr to user defined data to be
                                          passed to Get_Obj_Size()           */
  ZOLTAN_OBJ_SIZE_MULTI_FN *Get_Obj_Size_Multi;
                                       /* Function that returns the size of
                                          contiguous memory needed to store
                                          the data for multiple objects for
                                          migration.                         */
  ZOLTAN_OBJ_SIZE_MULTI_FORT_FN *Get_Obj_Size_Multi_Fort;
                                       /* Fortran version                    */
  void *Get_Obj_Size_Multi_Data;       /* Ptr to user defined data to be
                                          passed to Get_Obj_Size_Multi()     */
  /***************************************************************************/
  ZOLTAN_PACK_OBJ_FN *Pack_Obj;        /* Routine that packs object data for
                                          a given object into contiguous
                                          memory for migration.              */
  ZOLTAN_PACK_OBJ_FORT_FN *Pack_Obj_Fort;
                                       /* Fortran version                    */
  void *Pack_Obj_Data;                 /* Ptr to user defined data to be
                                          passed to Pack_Obj()               */
  ZOLTAN_PACK_OBJ_MULTI_FN *Pack_Obj_Multi;
                                       /* Routine that packes object data for
                                          multiple objects into contiguous
                                          memory for migration.              */
  ZOLTAN_PACK_OBJ_MULTI_FORT_FN *Pack_Obj_Multi_Fort;
                                       /* Fortran version                    */
  void *Pack_Obj_Multi_Data;                 /* Ptr to user defined data to be
                                          passed to Pack_Obj_Multi()         */


  /***************************************************************************/
  ZOLTAN_UNPACK_OBJ_FN *Unpack_Obj;    /* Routine that unpacks object data for
                                          a given object from contiguous
                                          memory after migration.            */
  ZOLTAN_UNPACK_OBJ_FORT_FN *Unpack_Obj_Fort;
                                       /* Fortran version                    */
  void *Unpack_Obj_Data;               /* Ptr to user defined data to be
                                          passed to Unpack_Obj()             */
  ZOLTAN_UNPACK_OBJ_MULTI_FN *Unpack_Obj_Multi;
                                       /* Routine that unpacks object data for
                                          multiple objects from contiguous
                                          memory after migration.            */
  ZOLTAN_UNPACK_OBJ_MULTI_FORT_FN *Unpack_Obj_Multi_Fort;
                                       /* Fortran version                    */
  void *Unpack_Obj_Multi_Data;         /* Ptr to user defined data to be
                                          passed to Unpack_Obj_Multi()       */
  /***************************************************************************/
  ZOLTAN_PROC_NAME_FN *Get_Processor_Name; 
                                       /* Fn ptr to get proc name   */
  void *Get_Processor_Name_Data;       /* Ptr to user defined data   */

  /***************************************************************************/
  struct Zoltan_LB_Struct LB;          /* Struct with info for load balancing */
  struct Zoltan_Order_Struct  Order;   /* Struct with info for ordering       */
  struct Zoltan_TPL_Order_Struct  TPL_Order; /* Struct with info for ordering */
  struct Zoltan_Migrate_Struct Migrate;/* Struct with info for migration.     */
};

typedef struct Zoltan_Struct ZZ;

/*
 *  A structure to hold coordinate transformations for degenerate
 *  geometries.
 */

struct Zoltan_Transform_Struct{
  int Target_Dim; /* Number of dimensions in transformed geometry:
                    -1 - not computed yet
                     0 - geometry is not degenerate, don't transform it
                     1 or 2 - geometry will be projected to line or plane   */
  double Transformation[3][3]; /* transforms degenerate geometry to 2D or 1D */
  int Permutation[3]; /* if trans. is simple coord switch, use this instead */
  double CM[3];       /* saved center of mass */
  double Evecs[3][3]; /* saved eigenvectors */
  int Axis_Order[3];  /* saved axis alignment information */

};

typedef struct Zoltan_Transform_Struct ZZ_Transform;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* MACROS  */

/*
 *  Test whether the processor is in the given Zoltan structure's
 *  communicator.  Used to exit from balancing routines for processors
 *  that are not included in the load-balancing communicator.
 */

#define ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz) ((zz)->Proc == -1) 

/*  
 *  Print trace information.
 */

#define ZOLTAN_TRACE_ENTER(zz,yo) do { \
  if ((zz)->Debug_Level >= ZOLTAN_DEBUG_PARAMS)  \
    Zoltan_add_back_trace((yo));           \
  if ((zz)->Debug_Level >= ZOLTAN_DEBUG_TRACE_ALL || \
     ((zz)->Proc == (zz)->Debug_Proc && \
      (zz)->Debug_Level == ZOLTAN_DEBUG_TRACE_SINGLE)) \
    ZOLTAN_TRACE_IN((zz)->Proc, (yo), NULL); } while (0)

#define ZOLTAN_TRACE_EXIT(zz,yo) do { \
  if ((zz)->Debug_Level >= ZOLTAN_DEBUG_PARAMS)  \
    Zoltan_remove_back_trace();           \
  if ((zz)->Debug_Level >= ZOLTAN_DEBUG_TRACE_ALL || \
     ((zz)->Proc == (zz)->Debug_Proc && \
      (zz)->Debug_Level == ZOLTAN_DEBUG_TRACE_SINGLE)) \
    ZOLTAN_TRACE_OUT((zz)->Proc, (yo), NULL); } while (0)

#define ZOLTAN_TRACE_DETAIL(zz,yo,string) do { \
  if ((zz)->Debug_Level >= ZOLTAN_DEBUG_TRACE_DETAIL) \
    ZOLTAN_PRINT_INFO((zz)->Proc, (yo), (string)); } while (0)


  /* Error Handling macro, used in PHG, coloring, matrix, graph ... */

#define MEMORY_ERROR do { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
} while (0)

#define FATAL_ERROR(s) do { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, s); \
  ierr = ZOLTAN_FATAL; \
  goto End; \
} while (0)

#define CHECK_FOR_MPI_ERROR(rc) do { \
  if (rc != MPI_SUCCESS){ \
    char _mpi_err_str[MPI_MAX_ERROR_STRING]; \
    int _mpi_err_len; \
    MPI_Error_string(rc, _mpi_err_str, &_mpi_err_len);  \
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, _mpi_err_str); \
    ierr = ZOLTAN_FATAL; \
    goto End; \
  } } while (0)


#define CHECK_IERR do {   if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) \
    goto End;  } while (0)

/*
 *  Debugging macro for Tflop architecture.
 *  ZOLTAN_HEAP_INFO(proc_number, string) prints information about the heap,
 *  such as number of fragments, total free memory, largest free chunk 
 *  of memory, and total used memory.  The processor number and provided
 *  string are printed to help instrument the code.
 *  On architectures other than Tflop, ZOLTAN_HEAP_INFO compiles 
 *  but has no effect.
 */
#ifdef TFLOP
#define ZOLTAN_HEAP_INFO(Proc,a) \
 {int frag, tfree, lfree, tused; \
  heap_info(&frag,&tfree,&lfree,&tused); \
  printf("HI%d %s frags = %d  tot free = %d  lar free = %d  tot used = %d\n", \
         Proc, a, frag, tfree, lfree, tused); \
 }
#else
#define ZOLTAN_HEAP_INFO(Proc,a) ;
#endif


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */
extern int Zoltan_Get_Obj_List(ZZ *, int *, 
              ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int, float**, int **);
extern int Zoltan_Get_Obj_List_Special_Malloc(ZZ *, int *, 
              ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int, float**, int **);

extern int Zoltan_Print_Obj_List( ZZ *zz, ZOLTAN_ID_PTR Gids, ZOLTAN_ID_PTR Lids,
  int wdim, float *Weights, int *Parts, int howMany);

extern int Zoltan_Get_Coordinates(ZZ *, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
  int *, double **);

extern void Zoltan_Initialize_Transformation(ZZ_Transform *tr);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
