/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifndef __ZOLTAN_CONST_H
#define __ZOLTAN_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __STDC__
#include <string.h>
#else
#include <strings.h>
#endif  /* __STDC__ */

#include "zoltan.h"
#include "lb_const.h"
#include "zz_id_const.h"
#include "zoltan_util.h"
#include "par_const.h"
#include "timer_const.h"

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
#define ZOLTAN_VER   1.53

/*
 * Type used to store linked list of new values for parameters.
 */
struct Param_List;


#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */


/*
 *  Define the debug levels allowed.
 *    ZOLTAN_DEBUG_NONE = 0           - quiet mode; no debugging information.
 *    ZOLTAN_DEBUG_PARAMS = 1         - print values of all parameters used.
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
  MachineType *Machine_Desc;      /*  Machine description for hetero. arch.  */
  struct Param_List *Params;      /*  List of parameter names & new vals     */
  int Deterministic;              /*  Flag indicating whether algorithms used
                                      should be forced to be deterministic.
                                      Default = TRUE.                        */
  int Obj_Weight_Dim;             /*  Dimension of the object weights, 
                                      usually 0 (no weights) or 1            */
  int Edge_Weight_Dim;            /*  Dimension of the edge weights, 
                                      usually 0 (no weights) or 1            */
  int Timer;                      /*  Timer type that is currently active */
  /***************************************************************************/
  ZOLTAN_PARTITION_MULTI_FN *Get_Partition_Multi;
                                       /* Fn ptr to get objects'
                                          partition assignments.     */
  ZOLTAN_PARTITION_MULTI_FORT_FN *Get_Partition_Multi_Fort;
                                       /* Fortran version            */
  void *Get_Partition_Multi_Data;      /* Ptr to user defined data to be 
                                          passed to Get_Partition_Multi()    */
  ZOLTAN_PARTITION_FN *Get_Partition;          
                                       /* Fn ptr to get an object's
                                          partition assignment.      */
  ZOLTAN_PARTITION_FORT_FN *Get_Partition_Fort;
                                       /* Fortran version            */
  void *Get_Partition_Data;            /* Ptr to user defined data
                                          to be passed to Get_Partition()    */
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
  ZOLTAN_NUM_HG_EDGES_FN *Get_Num_HG_Edges;    
                                       /* Fn ptr to get an object's
                                          number of hypergraph edges.*/
  ZOLTAN_NUM_HG_EDGES_FORT_FN *Get_Num_HG_Edges_Fort;
                                       /* Fortran version      */
  void *Get_Num_HG_Edges_Data;         /* Ptr to user defined data
                                          to be passed to
                                          Get_Num_HG_Edges()         */
  /***************************************************************************/
  ZOLTAN_HG_EDGE_LIST_FN *Get_HG_Edge_List;    
                                       /* Fn ptr to get an object's
                                          hyper-edge list.           */
  ZOLTAN_HG_EDGE_LIST_FORT_FN *Get_HG_Edge_List_Fort;/* Fortran version      */
  void *Get_HG_Edge_List_Data;         /* Ptr to user defined data
                                          to be passed to Get_HG_Edge_List() */
  /***************************************************************************/
  ZOLTAN_NUM_HG_PINS_FN *Get_Num_HG_Pins;      
                                       /* Fn ptr to get the total
                                          # pins in local hypergraph */
  ZOLTAN_NUM_HG_PINS_FORT_FN *Get_Num_HG_Pins_Fort; 
                                       /* Fortran version             */
  void *Get_Num_HG_Pins_Data;          /* Ptr to user defined data
                                           to be passed to Get_Num_HG_Pins() */
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
  ZOLTAN_GET_PROCESSOR_NAME_FN *Get_Processor_Name; 
                                       /* Fn ptr to get proc name   */
  void *Get_Processor_Name_Data;       /* Ptr to user defined data   */

  /***************************************************************************/
  struct Zoltan_LB_Struct LB;          /* Struct with info for load balancing */
  struct Zoltan_Migrate_Struct Migrate;/* Struct with info for migration.     */
};

typedef struct Zoltan_Struct ZZ;

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
#define ZOLTAN_TRACE_ENTER(zz,yo) \
  if ((zz)->Debug_Level >= ZOLTAN_DEBUG_TRACE_ALL || \
     ((zz)->Proc == (zz)->Debug_Proc && \
      (zz)->Debug_Level == ZOLTAN_DEBUG_TRACE_SINGLE)) \
    ZOLTAN_TRACE_IN((zz)->Proc, (yo), NULL);

#define ZOLTAN_TRACE_EXIT(zz,yo) \
  if ((zz)->Debug_Level >= ZOLTAN_DEBUG_TRACE_ALL || \
     ((zz)->Proc == (zz)->Debug_Proc && \
      (zz)->Debug_Level == ZOLTAN_DEBUG_TRACE_SINGLE)) \
    ZOLTAN_TRACE_OUT((zz)->Proc, (yo), NULL);

#define ZOLTAN_TRACE_DETAIL(zz,yo,string) \
  if ((zz)->Debug_Level >= ZOLTAN_DEBUG_TRACE_DETAIL) \
    ZOLTAN_PRINT_INFO((zz)->Proc, (yo), (string));

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
extern int Zoltan_Get_Obj_List(ZZ *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, 
                               int, float**, int **);

extern int Zoltan_Get_Coordinates(ZZ *, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
  int *, double **);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
