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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __STDC__
#include <string.h>
#else
#include <strings.h>
#endif  /* __STDC__ */

#include "lbi_const.h"

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
#define LB_VER   1.02


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
typedef struct LB_Tag_Struct LB_TAG;

typedef int LB_FN(LB *, int *, LB_GID **, LB_LID **, int **,
                        int *, LB_GID **, LB_LID **, int **);

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
  LB_MAX_METHODS                  /*  This entry should always be last.      */
} LB_METHOD;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Define a data structure for object information.
 */

struct LB_Tag_Struct {
  LB_GID Global_ID;               /* Global ID for the object; provided by 
                                     the application.                        */
  LB_LID Local_ID;                /* Local ID for the object; the application
                                     determines what is meant by "local ID";
                                     the load-balancer stores this field only
                                     so the application can take advantage of
                                     local indexing in the query routines.   */
  int   Proc;                     /* A processor ID for the tag.  Could be 
                                     the destination of an object to be 
                                     exported or the source of an object to
                                     be imported.                            */
};

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

  LB_PRE_MIGRATE_FN *Pre_Process;      /* Function that performs application
                                          specific pre-processing.  Optional
                                          for help-migration.                */
  LB_PRE_MIGRATE_FORT_FN *Pre_Process_Fort; /* Fortran version               */
  void *Pre_Process_Data;              /* Ptr to user defined data to be
                                          passed to Pre_Process()            */
  LB_POST_MIGRATE_FN *Post_Process;    /* Function that performs application
                                          specific post-processing.  Optional
                                          for help-migration.                */
  LB_POST_MIGRATE_FORT_FN *Post_Process_Fort; /* Fortran version             */
  void *Post_Process_Data;             /* Ptr to user defined data to be
                                          passed to Post_Process()           */
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
 
/**************************************************************/
/* The data structures below for HA are not being used yet!   */
/* They are merely Erik's suggestions.                        */
/**************************************************************/

/* Data structure for topology */
struct Topology_Type {
  int topo_id;
  int cart_dim[3];
  int wrap[3];
  int *xadj, *adjncy;  /* EB: Need a more general graph structure that allows for edge weights */
};

typedef struct Topology_Type TopologyType;

/* Data structure for machine */
struct Machine_Type {
  int nnodes;
  int cpu_power;
  int memory_size;
  TopologyType *network;
  struct Machine_Type **node_desc; /* Pointers to an array of nnodes descriptors */
};

typedef struct Machine_Type MachineType;

/* Machine type constants */
#define TOPO_FLAT    0
#define TOPO_SMP     0
#define TOPO_1D      1
#define TOPO_2D      2
#define TOPO_3D      3
#define TOPO_HCUBE   4
#define TOPO_GENERIC 5

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Define a load balancing object.  It will contain pointers to the
 *  appropriate functions for interfacing with applications and 
 *  pointers to the data structure used for load balancing.
 */

struct LB_Struct {
  MPI_Comm Communicator;          /*  The MPI Communicator.                  */
  int Proc;                       /*  The processor's ID within the MPI
                                      Communicator.                          */
  int Num_Proc;                   /*  The number of processors in the MPI
                                      Communicator.                          */
  int Debug_Level;                /*  Debug level for this instance of
                                      load balancing.                        */
  int Fortran;                    /*  1 if created from Fortran, 0 otherwise */
  MachineType *Machine_Desc;      /*  Machine description for hetero. arch. */
  LB_METHOD Method;               /*  Method to be used for load balancing.  */
  LB_FN *LB_Fn;                   /*  Pointer to the function that performs
                                      the load balancing; this ptr is set
                                      based on the method used.              */
  LB_PARAM *Params;               /*  List of parameter names & new vals */
  double Imbalance_Tol;           /*  Tolerance to which to load balance;
                                      Imbalance_Tol = 1.1 implies 10% imbalance
                                      is acceptable, i.e. max/avg = 1.1.     */
  int Obj_Weight_Dim;             /*  Dimension of the object weights, 
                                      usually 0 (no weights) or 1 */
  int Comm_Weight_Dim;            /*  Dimension of the communication weights, 
                                      usually 0 (no weights) or 1 */
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
 *  Test whether the processor is in the given load-balancing object's
 *  communicator.  Used to exit from balancing routines for processors
 *  that are not included in the load-balancing communicator.
 */

#define LB_PROC_NOT_IN_COMMUNICATOR(lb) ((lb)->Proc == -1) 

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

#include "par_const.h"

#endif
