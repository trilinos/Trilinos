/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __LB_CONST_H
#define __LB_CONST_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __STDC__
#include <string.h>
#else
#include <strings.h>
#endif  /* __STDC__ */

#include "lbi_const.h"
#include "lb_id_const.h"
#include "zoltan_util.h"
#include "mem_const.h"
#include "par_const.h"
#include "timer_const.h"

/*
 *  See bottom for other included files.
 */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Type definitions.
 */

/*
 * Strings to define the library name, and the version number
 * so that it is easier to keep track of what code each user
 * has.
 */
#define UTIL_NAME "zoltan"
#define LB_VER   1.241


/*
 * Type used to store linked list of new values for parameters.
 */
   
typedef struct LB_Param {
  char *name;
  char *new_val;
  struct LB_Param *next;
} LB_PARAM;
	  


#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */

typedef struct LB_Struct LB;

typedef int LB_FN(LB *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **,
                        int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **);

/*
 *  Define the possible load balancing methods allowed.
 */

typedef enum LB_Method {
  NONE = -1,
  RCB,
  OCTPART,
  PARMETIS,
  JOSTLE,
  REFTREE,
  RIB,
  SFC,
  LB_MAX_METHODS                  /*  This entry should always be last.      */
} LB_METHOD;

/*
 *  Define the debug levels allowed.
 *    LB_DEBUG_NONE = 0           - quiet mode; no debugging information.
 *    LB_DEBUG_PARAMS = 1         - print values of all parameters used.
 *    LB_DEBUG_ZTIME = 2          - print Zoltan timing information.
 *    LB_DEBUG_ATIME = 3          - print algorithm's timing info, if the
 *                                  algorithm supports this level.
 *    LB_DEBUG_TRACE_ZERO = 5     - print trace info on processor 0 only.
 *    LB_DEBUG_TRACE_ALL = 6      - print trace info on all processors.
 *    LB_DEBUG_TRACE_DETAIL = 7   - print detailed trace info on all processors.
 *    LB_DEBUG_LIST = 8           - print lists of objects to be imported 
 *                                  and exported.
 *    LB_DEBUG_ALL = 10           - print all debug information available.
 */
#define LB_DEBUG_NONE 0     
#define LB_DEBUG_PARAMS 1
#define LB_DEBUG_ZTIME 2
#define LB_DEBUG_ATIME 3
#define LB_DEBUG_TRACE_SINGLE 5
#define LB_DEBUG_TRACE_ALL 6
#define LB_DEBUG_TRACE_DETAIL 7 
#define LB_DEBUG_LIST 8
#define LB_DEBUG_ALL 10

/*
 * Values indicating which lists (import, export, both, or none) should
 * be returned by LB_Balance.  LB_NO_LISTS must always be zero; other
 * values should always be greater than zero.
 */
#define LB_NO_LISTS 0
#define LB_IMPORT_LISTS 1
#define LB_EXPORT_LISTS 2
#define LB_ALL_LISTS 3

/*
 ******************************************************
 * Define default values for key parameters.
 * These are used in both lb.c and key_params.c.
 ******************************************************
 */
#define LB_IMBALANCE_TOL_DEF  1.1
#define LB_DEBUG_LEVEL_DEF    LB_DEBUG_PARAMS
#define LB_DEBUG_PROC_DEF     0
#define LB_OBJ_WEIGHT_DEF     0
#define LB_COMM_WEIGHT_DEF    0
#define LB_AUTO_MIGRATE_DEF   FALSE
#define LB_DETERMINISTIC_DEF  TRUE
#define LB_NUM_ID_ENTRIES_DEF 1
#define LB_RETURN_LISTS_DEF   LB_ALL_LISTS
#define LB_TIMER_DEF          LB_TIME_WALL
#define LB_TFLOPS_SPECIAL_DEF FALSE

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* 
 *  Define a communication structure for load balancing results.  This
 *  structure will be used by Steve Plimpton's and Bruce Hendrickson's
 *  communication library.
 */

struct LB_Migrate_Struct {
  int Auto_Migrate;                   /*  Flag indicating whether the load
                                          balancer should automatically
                                          help the application
                                          migrate data.  Some applications may
                                          prefer to do it themselves.        */
  /*
   *  Pointers to routines that depend on the application.
   */

  LB_PRE_MIGRATE_FN *Pre_Migrate;      /* Function that performs application
                                          specific pre-processing.  Optional
                                          for help-migration.                */
  LB_PRE_MIGRATE_FORT_FN *Pre_Migrate_Fort; /* Fortran version               */
  void *Pre_Migrate_Data;              /* Ptr to user defined data to be
                                          passed to Pre_Migrate()            */
  LB_MID_MIGRATE_FN *Mid_Migrate;      /* Function that performs application
                                          specific processing between packing
                                          and unpacking.  Optional
                                          for help-migration.                */
  LB_MID_MIGRATE_FORT_FN *Mid_Migrate_Fort; /* Fortran version               */
  void *Mid_Migrate_Data;              /* Ptr to user defined data to be
                                          passed to Mid_Migrate()            */
  LB_POST_MIGRATE_FN *Post_Migrate;    /* Function that performs application
                                          specific post-processing.  Optional
                                          for help-migration.                */
  LB_POST_MIGRATE_FORT_FN *Post_Migrate_Fort; /* Fortran version             */
  void *Post_Migrate_Data;             /* Ptr to user defined data to be
                                          passed to Post_Migrate()           */
  LB_OBJ_SIZE_FN *Get_Obj_Size;        /* Function that returns the size of
                                          contiguous memory needed to store
                                          the data for a single object for
                                          migration.                         */
  LB_OBJ_SIZE_FORT_FN *Get_Obj_Size_Fort; /* Fortran version                 */
  void *Get_Obj_Size_Data;             /* Ptr to user defined data to be
                                          passed to Get_Obj_Size()           */
  LB_PACK_OBJ_FN *Pack_Obj;            /* Routine that packs object data for
                                          a given object into contiguous 
                                          memory for migration.              */
  LB_PACK_OBJ_FORT_FN *Pack_Obj_Fort;  /* Fortran version                    */
  void *Pack_Obj_Data;                 /* Ptr to user defined data to be
                                          passed to Pack_Obj()               */
  LB_UNPACK_OBJ_FN *Unpack_Obj;        /* Routine that unpacks object data for
                                          a given object from contiguous 
                                          memory after migration.            */
  LB_UNPACK_OBJ_FORT_FN *Unpack_Obj_Fort; /* Fortran version                 */
  void *Unpack_Obj_Data;               /* Ptr to user defined data to be
                                          passed to Unpack_Obj()             */
};

typedef struct LB_Migrate_Struct LB_MIGRATE;


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
 *  Define a load balancing structure.  It will contain pointers to the
 *  appropriate functions for interfacing with applications and 
 *  pointers to the data structure used for load balancing.
 */

struct LB_Struct {
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
  int Return_Lists;               /*  Flag indicating which lists (if any)
                                      should be returned by LB_Balance.      */
  int Tflops_Special;             /*  Flag to indicate if we should use some
                                      MPI constructs (0) or not (1) on tflops*/
  MachineType *Machine_Desc;      /*  Machine description for hetero. arch.  */
  LB_METHOD Method;               /*  Method to be used for load balancing.  */
  LB_FN *LB_Fn;                   /*  Pointer to the function that performs
                                      the load balancing; this ptr is set
                                      based on the method used.              */
  LB_PARAM *Params;               /*  List of parameter names & new vals     */
  double Imbalance_Tol;           /*  Tolerance to which to load balance;
                                      Imbalance_Tol = 1.1 implies 10% imbalance
                                      is acceptable, i.e. max/avg = 1.1.     */
  int Deterministic;              /*  Flag indicating whether algorithms used
                                      should be forced to be deterministic.
                                      Default = TRUE.                        */
  int Obj_Weight_Dim;             /*  Dimension of the object weights, 
                                      usually 0 (no weights) or 1            */
  int Comm_Weight_Dim;            /*  Dimension of the communication weights, 
                                      usually 0 (no weights) or 1            */
  int Timer;                      /*  Timer type that is currently active */
  void *Data_Structure;           /*  Data structure used by the load 
                                      balancer; cast by the method routines
                                      to the appropriate data type.          */
  LB_NUM_EDGES_FN *Get_Num_Edges;              /* Fn ptr to get an object's
                                                  number of edges.           */
  LB_NUM_EDGES_FORT_FN *Get_Num_Edges_Fort;    /* Fortran version            */
  void *Get_Num_Edges_Data;                    /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Num_Edges()            */
  LB_EDGE_LIST_FN *Get_Edge_List;              /* Fn ptr to get an object's
                                                  edge list.                 */
  LB_EDGE_LIST_FORT_FN *Get_Edge_List_Fort;    /* Fortran version            */
  void *Get_Edge_List_Data;                    /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Edge_List()            */
  LB_NUM_GEOM_FN *Get_Num_Geom;                /* Fn ptr to get an object's
                                                  number of geometry values. */
  LB_NUM_GEOM_FORT_FN *Get_Num_Geom_Fort;      /* Fortran version            */
  void *Get_Num_Geom_Data;                     /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Num_Geom()             */
  LB_GEOM_FN *Get_Geom;                        /* Fn ptr to get an object's
                                                  geometry values.           */
  LB_GEOM_FORT_FN *Get_Geom_Fort;              /* Fortran version            */
  void *Get_Geom_Data;                         /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Geom()                 */
  LB_NUM_OBJ_FN *Get_Num_Obj;                  /* Fn ptr to get a proc's  
                                                  number of local objects.   */
  LB_NUM_OBJ_FORT_FN *Get_Num_Obj_Fort;        /* Fortran version            */
  void *Get_Num_Obj_Data;                      /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Num_Obj()              */
  LB_OBJ_LIST_FN *Get_Obj_List;                /* Fn ptr to get all local
                                                  objects on a proc.         */
  LB_OBJ_LIST_FORT_FN *Get_Obj_List_Fort;      /* Fortran version            */
  void *Get_Obj_List_Data;                     /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Obj_List()             */
  LB_FIRST_OBJ_FN *Get_First_Obj;              /* Fn ptr to get the first   
                                                  local obj on a proc.       */
  LB_FIRST_OBJ_FORT_FN *Get_First_Obj_Fort;    /* Fortran version            */
  void *Get_First_Obj_Data;                    /* Ptr to user defined data
                                                  to be passed to
                                                  Get_First_Obj()            */
  LB_NEXT_OBJ_FN *Get_Next_Obj;                /* Fn ptr to get the next   
                                                  local obj on a proc.       */
  LB_NEXT_OBJ_FORT_FN *Get_Next_Obj_Fort;      /* Fortran version            */
  void *Get_Next_Obj_Data;                     /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Next_Obj()             */
  LB_NUM_BORDER_OBJ_FN *Get_Num_Border_Obj;    /* Fn ptr to get a proc's 
                                                  number of border objs wrt
                                                  a given processor.         */
  LB_NUM_BORDER_OBJ_FORT_FN *Get_Num_Border_Obj_Fort; /* Fortran version     */ 
  void *Get_Num_Border_Obj_Data;               /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Num_Border_Obj()       */
  LB_BORDER_OBJ_LIST_FN *Get_Border_Obj_List;  /* Fn ptr to get all objects
                                                  sharing a border with a
                                                  given processor.           */
  LB_BORDER_OBJ_LIST_FORT_FN *Get_Border_Obj_List_Fort; /* Fortran version   */
  void *Get_Border_Obj_List_Data;              /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Border_Obj_List()      */
  LB_FIRST_BORDER_OBJ_FN *Get_First_Border_Obj;/* Fn ptr to get the first 
                                                  object sharing a border 
                                                  with a given processor.    */
  LB_FIRST_BORDER_OBJ_FORT_FN *Get_First_Border_Obj_Fort; /* Fortran version */
  void *Get_First_Border_Obj_Data;             /* Ptr to user defined data
                                                  to be passed to
                                                  Get_First_Border_Obj()     */
  LB_NEXT_BORDER_OBJ_FN *Get_Next_Border_Obj;  /* Fn ptr to get the next 
                                                  object sharing a border 
                                                  with a given processor.    */
  LB_NEXT_BORDER_OBJ_FORT_FN *Get_Next_Border_Obj_Fort; /* Fortran version   */
  void *Get_Next_Border_Obj_Data;              /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Next_Border_Obj()      */
  LB_NUM_COARSE_OBJ_FN *Get_Num_Coarse_Obj;    /* Fn ptr to get the number of
                                                  elements in the coarse grid*/
  LB_NUM_COARSE_OBJ_FORT_FN *Get_Num_Coarse_Obj_Fort; /* Fortran version     */
  void *Get_Num_Coarse_Obj_Data;               /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Num_Coarse_Obj()       */
  LB_COARSE_OBJ_LIST_FN *Get_Coarse_Obj_List;  /* Fn ptr to get all
                                                  elements in the coarse grid*/
  LB_COARSE_OBJ_LIST_FORT_FN *Get_Coarse_Obj_List_Fort; /* Fortran version   */
  void *Get_Coarse_Obj_List_Data;              /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Coarse_Obj_List()      */
  LB_FIRST_COARSE_OBJ_FN *Get_First_Coarse_Obj;/* Fn ptr to get the first coarse
                                                  obj on a proc.             */
  LB_FIRST_COARSE_OBJ_FORT_FN *Get_First_Coarse_Obj_Fort; /* Fortran version */
  void *Get_First_Coarse_Obj_Data;             /* Ptr to user defined data
                                                  to be passed to
                                                  Get_First_Coarse_Obj()     */
  LB_NEXT_COARSE_OBJ_FN *Get_Next_Coarse_Obj;  /* Fn ptr to get the next coarse
                                                  obj on a proc.             */
  LB_NEXT_COARSE_OBJ_FORT_FN *Get_Next_Coarse_Obj_Fort; /* Fortran version   */
  void *Get_Next_Coarse_Obj_Data;              /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Next_Coarse_Obj()      */
  LB_NUM_CHILD_FN *Get_Num_Child;              /* Fn ptr to get the number of
                                                  children of an element     */
  LB_NUM_CHILD_FORT_FN *Get_Num_Child_Fort;    /* Fortran version            */
  void *Get_Num_Child_Data;                    /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Num_Child()            */
  LB_CHILD_LIST_FN *Get_Child_List;            /* Fn ptr to get all
                                                  children of an element     */
  LB_CHILD_LIST_FORT_FN *Get_Child_List_Fort;  /* Fortran version            */
  void *Get_Child_List_Data;                   /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Child_List()           */
  LB_CHILD_WEIGHT_FN *Get_Child_Weight;        /* Fn ptr to get the weight
                                                  of an element              */
  LB_CHILD_WEIGHT_FORT_FN *Get_Child_Weight_Fort; /* Fortran version         */
  void *Get_Child_Weight_Data;                 /* Ptr to user defined data
                                                  to be passed to
                                                  Get_Child_Weight()         */
  LB_MIGRATE Migrate;                          /* Struct with info for helping
                                                  with migration.            */
  void *Migrate_Data;                          /* Ptr to user defined data
                                                  to be passed to
                                                  Migrate()                  */
  LB_GET_PROCESSOR_NAME_FN *Get_Processor_Name; /* Fn ptr to get proc name   */
  void *Get_Processor_Name_Data;               /* Ptr to user defined data   */
};
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* MACROS  */

/*
 *  Test whether the processor is in the given load-balancing structure's
 *  communicator.  Used to exit from balancing routines for processors
 *  that are not included in the load-balancing communicator.
 */

#define LB_PROC_NOT_IN_COMMUNICATOR(lb) ((lb)->Proc == -1) 

/*
 *  Macros for consistently printing error and warning messages.
 */

#define LB_PRINT_ERROR(proc,yo,str) \
  ZOLTAN_PRINT_ERROR(proc, yo, str)

#define LB_PRINT_WARN(proc,yo,str) \
  ZOLTAN_PRINT_WARN(proc, yo, str)

/*  
 *  Print trace information.
 */
#define LB_TRACE_ENTER(lb,yo) \
  if ((lb)->Debug_Level >= LB_DEBUG_TRACE_ALL || \
     ((lb)->Proc == (lb)->Debug_Proc && \
      (lb)->Debug_Level == LB_DEBUG_TRACE_SINGLE)) \
    ZOLTAN_TRACE_ENTER((lb)->Proc, (yo), NULL);

#define LB_TRACE_EXIT(lb,yo) \
  if ((lb)->Debug_Level >= LB_DEBUG_TRACE_ALL || \
     ((lb)->Proc == (lb)->Debug_Proc && \
      (lb)->Debug_Level == LB_DEBUG_TRACE_SINGLE)) \
    ZOLTAN_TRACE_EXIT((lb)->Proc, (yo), NULL);

#define LB_TRACE_DETAIL(lb,yo,string) \
  if ((lb)->Debug_Level >= LB_DEBUG_TRACE_DETAIL) \
    ZOLTAN_PRINT_INFO((lb)->Proc, (yo), (string));

/*
 *  Debugging macro for Tflop architecture.
 *  LB_HEAP_INFO(proc_number, string) prints information about the heap,
 *  such as number of fragments, total free memory, largest free chunk 
 *  of memory, and total used memory.  The processor number and provided
 *  string are printed to help instrument the code.
 *  On architectures other than Tflop, LB_HEAP_INFO compiles but has no effect.
 */
#ifdef TFLOP
#define LB_HEAP_INFO(Proc,a) \
 {int frag, tfree, lfree, tused; \
  heap_info(&frag,&tfree,&lfree,&tused); \
  printf("HI%d %s frags = %d  tot free = %d  lar free = %d  tot used = %d\n", \
         Proc, a, frag, tfree, lfree, tused); \
 }
#else
#define LB_HEAP_INFO(Proc,a) ;
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */

extern LB_FN LB_rcb;
extern LB_FN LB_octpart;
extern LB_FN LB_ParMetis;
extern LB_FN LB_Jostle;
extern LB_FN LB_Reftree_Part;
extern LB_FN LB_rib;
extern LB_FN LB_sfc;


#endif
