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

#ifndef __LBI_CONST_H
#define __LBI_CONST_H

#ifndef lint
static char *cvs_lbiconsth_id = "$Id$";
#endif

#include <mpi.h>

/*
 *  Include user-defined data types and comparison macros for LB_GID and LB_LID.
 */
#include "lb_user_const.h"

/*
 *  Data types and functions describing the interface between the
 *  application and the load balancing tool.
 */

/*
 *  Enumerated type used to indicate which function is to be set by
 *  LB_Set_Fn.
 */

enum LB_Fn_Type {
  LB_NUM_EDGES_FN_TYPE,
  LB_EDGE_LIST_FN_TYPE,
  LB_NUM_GEOM_FN_TYPE,
  LB_GEOM_FN_TYPE,
  LB_NUM_OBJ_FN_TYPE,
  LB_OBJ_LIST_FN_TYPE,
  LB_FIRST_OBJ_FN_TYPE,
  LB_NEXT_OBJ_FN_TYPE,
  LB_NUM_BORDER_OBJ_FN_TYPE,
  LB_BORDER_OBJ_LIST_FN_TYPE,
  LB_FIRST_BORDER_OBJ_FN_TYPE,
  LB_NEXT_BORDER_OBJ_FN_TYPE,
  LB_PRE_MIGRATE_FN_TYPE,
  LB_OBJ_SIZE_FN_TYPE,
  LB_PACK_OBJ_FN_TYPE,
  LB_UNPACK_OBJ_FN_TYPE,
  LB_MAX_FN_TYPES               /*  This entry should always be last.        */
};

typedef enum LB_Fn_Type LB_FN_TYPE;

/*
 *  Other common definitions:
 */

struct LB_Struct;


/*
 *  Maximum number of parameters to be passed to any load-balancing
 *  method.
 */

#define LB_PARAMS_MAX_SIZE 7

/*
 * Error codes for DLB library
 *   DLB_OK     - no errors
 *   DLB_WARN   - some warning occurred in DLB library; application should be
 *                able to continue running
 *   DLB_FATAL  - a fatal error occurred
 *   DLB_MEMERR - memory allocation failed; with this error, it could be
 *                possible to try a different, more memory-friendly, algorithm
 */
#define DLB_OK     0
#define DLB_WARN   1
#define DLB_FATAL  -1
#define DLB_MEMERR -2

/* DLB_ should be replaced with LB_ to be consistent with 
 * the rest of Zoltan. For now, keep both sets of definitions.
 */
#define LB_OK      0
#define LB_WARN    1
#define LB_FATAL  -1
#define LB_MEMERR -2

/*****************************************************************************/
/*****************************************************************************/
/**********************  Functions to query application  *********************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Function to return, for the object with a given ID,
 *  the object's number of edges (i.e., the number of objects with which
 *  the given object must communicate).
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    LB_GID global_id          --  the Global ID for the object
 *    LB_LID local_id           --  the Local ID for the object
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the number of neighbor objects.
 */

typedef int LB_NUM_EDGES_FN(void *data, LB_GID global_id, LB_LID local_id,
                            int *ierr);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID, 
 *  the object's edge list (i.e., objects with which the given object must
 *  communicate.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    LB_GID global_id          --  the Global ID for the object
 *    LB_LID local_id           --  the Local ID for the object
 *    int    wdim               --  dimension of edge weights, or 0 if
 *                                  edge weights are not sought.
 *  Output:
 *    LB_GID *nbor_global_ids   --  Array of Global IDs of neighboring objects.
 *    int    *nbor_procs        --  Array of neighboring procs.
 *    int    *nbor_ewgts        --  Array of edge weights, where 
 *                                  nbor_ewgts[i*wdim:(i+1)*wdim-1]
 *                                  corresponds to the weight of edge i
 *    int *ierr                 --  error code
 */

typedef void LB_EDGE_LIST_FN(void *data, LB_GID global_id, LB_LID local_id,
                             LB_GID *nbor_global_id, int *nbor_procs,
                             int wdim, int *nbor_ewgts, int *ierr);

/*****************************************************************************/
/*
 *  Function to return
 *  the number of geometry fields per object (e.g., the number of values
 *  used to express the coordinates of the object).
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the number of geometry fields.
 */

typedef int LB_NUM_GEOM_FN(void *data, int *ierr);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID,
 *  the geometry information for the object (e.g., coordinates).
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    LB_GID global_id          --  the Global ID for the object
 *    LB_LID local_id           --  the Local ID for the object
 *  Output:
 *    double *geom_vec          --  the geometry info for the object
 *                                  (e.g., coordinates)
 *    int *ierr                 --  error code
 */

typedef void LB_GEOM_FN(void *data, LB_GID global_id, LB_LID local_id,
                        double *geom_vec, int *ierr);

/*****************************************************************************/
/*
 *  Function to return, for the calling processor, the number of objects 
 *  located in that processor's memory.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the number of local objects.
 */

typedef int LB_NUM_OBJ_FN(void *data, int *ierr);

/*****************************************************************************/
/*
 *  Function to return a list of all local objects on a processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int wdim                  --  dimension of object weights, or 0 if
 *                                  object weights are not sought. 
 *  Output:
 *    LB_GID *global_ids        --  array of Global IDs of all objects on the
 *                                  processor.
 *    LB_LID *local_ids         --  array of Local IDs of all objects on the
 *                                  processor.
 *    float *objwgts            --  objwgts[i*wdim:(i+1)*wdim-1] correponds
 *                                  to the weight of object i 
 *    int *ierr                 --  error code
 */

typedef void LB_OBJ_LIST_FN(void *data, LB_GID *global_ids, LB_LID *local_ids,
                            int wdim, float *objwgts, int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for local objects; return the first local object on
 *  the processor.  This function should be used with LB_NEXT_OBJ_FN.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *  Output:
 *    LB_GID *first_global_id   --  Global ID of the first object; NULL if no
 *                                  objects.
 *    LB_LID *first_local_id    --  Local ID of the first object; NULL if no
 *                                  objects.
 *    float *first_obj_wgt      --  weight vector for first object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist on the processor.
 */

typedef int LB_FIRST_OBJ_FN(void *data, LB_GID *first_global_id,
                            LB_LID *first_local_id, 
                            int wdim, float *first_obj_wgt, int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for local objects; return the next local object.
 *  This function should be used with LB_FIRST_OBJ_FN.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    LB_GID global_id          --  Global ID of the previous object.
 *    LB_LID local_id           --  Local ID of the previous object.
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *  Output:
 *    LB_GID *next_global_id    --  Global ID of the next object; NULL if no
 *                                  more objects.
 *    LB_LID *next_local_id     --  Local ID of the next object; NULL if no
 *                                  more objects.
 *    float *next_obj_wgt       --  weight vector for the next object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist (i.e., global_id is
 *                                  the last object).
 */

typedef int LB_NEXT_OBJ_FN(void *data, LB_GID global_id, LB_LID local_id,
                           LB_GID *next_global_id, LB_LID *next_local_id,
                           int wdim, float *next_obj_wgt, int *ierr);

/*****************************************************************************/
/*
 *  Function to return, for the calling processor, the number of objects 
 *  sharing a subdomain border with a given processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int nbor_proc             --  processor ID of the neighboring processor.
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the number of local objects.
 */

typedef int LB_NUM_BORDER_OBJ_FN(void *data, int nbor_proc, int *ierr);

/*****************************************************************************/
/*
 *  Function to return a list of all objects sharing a subdomain border 
 *  with a given processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int nbor_proc             --  processor ID of the neighboring processor.
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weights are not sought.
 *  Output:
 *    LB_GID *global_ids        --  array of Global IDs of all objects on the
 *                                  processor border with the given neighboring
 *                                  processor.
 *    LB_LID *local_ids         --  array of Local IDs of all objects on the 
 *                                  processor border with the given neighboring 
 *                                  processor.
 *    float *objwgts            --  objwgts[i*wdim:(i+1)*wdim-1] correponds
 *                                  to the weight of object i 
 *                                  (objwgts is undefined if wdim=0)
 *    int *ierr                 --  error code
 */

typedef void LB_BORDER_OBJ_LIST_FN(void *data, int nbor_proc,
                                   LB_GID *global_ids, LB_LID *local_ids,
                                   int wdim, float *objwgts, int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for border objects; return the first local object 
 *  along the subdomain boundary with a given processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int nbor_proc             --  processor ID of the neighboring processor.
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *  Output:
 *    LB_GID *first_global_id   --  Global ID of the first object; NULL if no
 *                                  objects.
 *    LB_LID *first_local_id    --  Local ID of the first object; NULL if no 
 *                                  objects.
 *    float *first_obj_wgt      --  weight vector for the first object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist (i.e., global_id is
 *                                  the last object).
 */

typedef int LB_FIRST_BORDER_OBJ_FN(void *data, int nbor_proc,
                                   LB_GID *first_global_id,
                                   LB_LID *first_local_id, 
                                   int wdim, float *first_obj_wgt,
                                   int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for border objects; return the next local object 
 *  along the subdomain boundary with a given processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    LB_GID global_id          --  Global ID of the previous object.
 *    LB_LID local_id           --  Local ID of the previous object.
 *    int nbor_proc             --  processor ID of the neighboring processor.
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *  Output:
 *    LB_GID *next_global_id    --  Global ID of the next object; NULL if no
 *                                  more objects.
 *    LB_LID *next_local_id     --  Local ID of the next object; NULL if no 
 *                                  more objects.
 *    float *next_obj_wgt       --  weight vector for the next object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist (i.e., global_id is
 *                                  the last object).
 */

typedef int LB_NEXT_BORDER_OBJ_FN(void *data, LB_GID global_id,
                                  LB_LID local_id, int nbor_proc,
                                  LB_GID *next_global_id,
                                  LB_LID *next_local_id, 
                                  int wdim, float *next_obj_wgt,
                                  int *ierr);

/*****************************************************************************/
/*
 *  Function to return the size (in bytes) of data to be migrated.
 *  This function is needed only when the application
 *  wants the load-balancer to help migrate the data.  It is used by the
 *  comm.c routines to allocate message buffers.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the size of data of local objects.
 */

typedef int LB_OBJ_SIZE_FN(void *data, int *ierr);

/*****************************************************************************/
/*
 *  Function called as a pre-processor to the migration.  This function is 
 *  optional, and is used only when the application wants the load-balancer 
 *  to help migrate the data.  The application can perform any type of 
 *  pre-processing in this function.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_import            --  Number of objects to be imported.
 *    LB_GID *import_global_ids --  Global IDs of objects to be imported.
 *    LB_LID *import_local_ids  --  Local IDs of objects to be imported.
 *    int *import_procs         --  Processor IDs of importing processors.
 *    int num_export            --  Number of objects to be exported.
 *    LB_GID *export_global_ids --  Global IDs of objects to be exported.
 *    LB_LID *export_local_ids  --  Local IDs of objects to be exported.
 *    int *export_procs         --  Processor IDs of processors to receive
 *                                  the objects.
 *  Output:
 *    int *ierr                 --  error code
 */

typedef void LB_PRE_MIGRATE_FN(void *data, int num_import,
                               LB_GID *import_global_ids,
                               LB_LID *import_local_ids, int *import_procs,
                               int num_export, LB_GID *export_global_ids,
                               LB_LID *export_local_ids, int *export_procs,
                               int *ierr);

/*****************************************************************************/
/*
 *  Function to pack data to be migrated for the given object.
 *  This function is needed only when the application wants the load-balancer 
 *  to help migrate the data.  It packs all data related to the given object
 *  into a communication buffer, the starting address of which is provided
 *  by the load-balancer.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    LB_GID global_id          --  Global ID of the object to be packed.
 *    LB_LID local_id           --  Local ID of the object to be packed.
 *    int dest_proc             --  Processor ID of the destination processor.
 *    int size                  --  number of bytes allowed for the object to
 *                                  be packed.
 *    char *buf                 --  starting address of buffer into which to
 *                                  pack the object.
 *  Output:
 *    char *buf                 --  the buffer is rewritten with the packed
 *                                  data.
 *    int *ierr                 --  error code
 */

typedef void LB_PACK_OBJ_FN(void *data, LB_GID global_id, LB_LID local_id,
                            int dest_proc, int size, char *buf, int *ierr);

/*****************************************************************************/
/*
 *  Function to unpack data for an object migrated to a new processor.
 *  This function is needed only when the application wants the load-balancer 
 *  to help migrate the data.  The data is stored in a buffer (char *); the
 *  size of the data for the object is included.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    LB_GID global_id          --  Global ID of the object to be unpacked.
 *    int size                  --  number of bytes in the buffer for the
 *                                  object.
 *    char *buf                 --  starting address of buffer into which to
 *                                  pack the object.
 *  Output:
 *    int *ierr                 --  error code
 */

typedef void LB_UNPACK_OBJ_FN(void *data, LB_GID global_id, int size,
                              char *buf, int *ierr);

/*****************************************************************************/
/*****************************************************************************/
/**********************  Functions to set-up LB object ***********************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Function to initialize values needed in load balancing tools, and
 *  returns which version of the library this is. If the application
 *  uses MPI, call this function after calling MPI_Init. If the
 *  application does not use MPI, this function calls MPI_Init for
 *  use by the load balancer. This function returns the version of
 *  the DLB library.
 *  Input:
 *    int argc                   --  Argument count from main()
 *    char **argv                --  Argument list from main()
 *  Output:
 *    float *ver                 --  Version of DLB library
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Initialize(int argc, char **argv, float *ver);

/*****************************************************************************/
/*
 *  Function to create a load balancing object.  May want more than one
 *  object if using different decompositions with different techniques.
 *  This function allocates and initializes the object.
 *  Input:
 *    MPI_Comm communicator      --  MPI Communicator to be used for this
 *                                   load-balancing object.
 *    KDD_DLB  --  The communicator is not yet used in the algorithms!
 *    KDD_DLB  --  It will have to be incorporated appropriately.
 *    KDD_DLB  --  But I wanted to get it into the interface now!
 *  Returned value:
 *    struct LB_Struct *         --  Pointer to a LB object.
 *                                   If there is an error, NULL is returned.
 *                                   Any error in this function should be
 *                                   considered fatal.
 */

extern struct LB_Struct *LB_Create_Object(MPI_Comm communicator);

/*****************************************************************************/
/*
 *  Function to initialize a given LB interface function.
 *  Input:
 *    struct LB_Struct *lb       --  Pointer to a LB object.
 *    LB_FN_TYPE fn_type         --  Enum type indicating the function to be
 *                                   set.
 *    void *()fn_ptr             --  Pointer to the function to be used in the 
 *                                   assignment.
 *    void *data_ptr             --  Pointer to data that the DLB library will
 *                                   pass as an argument to fn(). May be NULL.
 *  Output:
 *    struct LB_Struct *lb       --  Appropriate field set to value in fn_ptr.
 *  Returned value:
 *    int                        --  Error code
 */
extern int LB_Set_Fn(struct LB_Struct *lb, LB_FN_TYPE fn_type,
                     void *fn_ptr(), void *data_ptr);

/*****************************************************************************/
/*
 *  Function to set the load balancing method to be used.
 *  Input:
 *    struct LB_Struct *lb       --  The load balancing object to which this
 *                                   method applies.
 *    char *string               --  String specifying the desired method.
 *    double *params             --  Params needed by desired method.
 *                                   (This field depends upon the particular
 *                                   method.)
 *  Output:
 *    struct LB_Struct *lb       --  Appropriate fields set to designated
 *                                   values.
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Set_Method(struct LB_Struct *lb, char *string, double *params);

/*****************************************************************************/
/*
 *  Function to set the tolerance to which the system must be load balanced.
 *  For example, if the tolerance is set to 0.9, 10% load imbalance between
 *  the most heavily loaded processor and the average load will be accepted
 *  as balanced.
 *  Input:
 *    struct LB_Struct *lb       --  The load balancing object to which this 
 *                                   tolerance applies.
 *    double tolerance           --  The tolerance desired.
 *  Output:
 *    struct LB_Struct *lb       --  Appropriate fields set to designated value.
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Set_Tolerance(struct LB_Struct *lb, double tolerance);

/*****************************************************************************/
/*
 *  Function to set a flag indicating whether the application wants the
 *  load-balancer to help with data migration.   If migration help is
 *  wanted, routines to pack and unpack object data must be provided by
 *  the application (see LB_OBJ_SIZE_FN, LB_PACK_OBJ_FN, LB_UNPACK_OBJ_FN).
 *
 *  Input:
 *    struct LB_Struct *lb       --  The load balancing object to which this 
 *                                   flag applies.
 *    int auto_migrate_flag      --  TRUE or FALSE to indicate whether the 
 *                                   application wants migration help.
 *                                   Default is FALSE.
 *  Output:
 *    struct LB_Struct *lb       --  Appropriate fields set to designated type.
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Set_Migration(struct LB_Struct *lb, int auto_migrate_flag);

/*****************************************************************************/
/*
 *  Function to initialize an array to pass parameters to the load-balancing
 *  methods.  This function is provided so that the load-balancer can
 *  look for array entries not set by the application and use default values
 *  for those entries.
 *
 *  Input/Output:
 *    double *params             --  Pointer to the array to be used to pass
 *                                   parameters to the load-balancing methods.
 *                                   Upon return, the values in this array are
 *                                   initialized to an initial value determined
 *                                   by the load-balancer.
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Initialize_Params_Array(double *params);

/*****************************************************************************/
/*
 *  Function to invoke the load-balancer.
 *
 *  Input:
 *    struct LB_Struct *lb       --  The load balancing object containing info 
 *                                   about this load-balancing invocation.
 *  Output:
 *    int *changes               --  This value tells if the load balancer
 *                                   came up with a new decomposition or
 *                                   not. It can be either a one or a zero:
 *                                   zero - No changes to the decomposition
 *                                          were made by the load-balancing
 *                                          algorithm; migration is not needed.
 *                                   one  - A new decomposition is suggested
 *                                          by the load-balancer; migration
 *                                          is needed to establish the new
 *                                          decomposition.
 *    int *num_import            --  The number of non-local objects in the 
 *                                   processor's new decomposition (i.e.,
 *                                   number of objects to be imported).
 *    LB_GID **import_global_ids --  Pointer to array of Global IDs for the
 *                                   objects to be imported.
 *    LB_LID **import_local_ids  --  Pointer to array of Local IDs for the 
 *                                   objects to be imported (local to the
 *                                   exporting processor).
 *    int **import_procs         --  Pointer to array of Processor IDs for the 
 *                                   objects to be imported (processor IDs of
 *                                   source processor).
 *    int *num_export            --  The number of local objects that need to be
 *                                   exported from the processor to establish
 *                                   the new decomposition.
 *    LB_GID **export_global_ids --  Pointer to array of Global IDs for the
 *                                   objects to be exported from the current
 *                                   processor.
 *    LB_LID **export_local_ids  --  Pointer to array of Local IDs for the
 *                                   objects to be exported (local to the
 *                                   current processor).
 *    int **export_procs         --  Pointer to array of Processor IDs for the
 *                                   objects to be exported (processor IDs of
 *                                   destination processors).
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Balance(struct LB_Struct *lb, int *changes,
                      int *num_import, LB_GID **import_global_ids,
                      LB_LID **import_local_ids, int **import_procs,
                      int *num_export, LB_GID **export_global_ids,
                      LB_LID **export_local_ids, int **export_procs);

/*****************************************************************************/
/*
 *  Routine to compute the inverse map:  Given, for each processor, a list
 *  of non-local objects assigned to the processor, compute the list of objects
 *  that processor needs to export to other processors to establish the new
 *  decomposition.
 *
 *  Input:
 *    struct LB_Struct *lb       --  Load balancing object for current balance.
 *    int num_import             --  Number of non-local objects assigned to the
 *                                   processor in the new decomposition.
 *    LB_GID *import_global_ids  --  Array of global IDs for non-local objects
 *                                   assigned to this processor in the new
 *                                   decomposition.
 *    LB_LID *import_local_ids   --  Array of local IDs for non-local objects
 *                                   assigned to the processor in the new
 *                                   decomposition.
 *    int *import_procs          --  Array of processor IDs of processors owning
 *                                   the non-local objects that are assigned to
 *                                   this processor in the new decomposition.
 *  Output:
 *    int *num_export            --  The number of local objects that need to be
 *                                   exported from the processor to establish
 *                                   the new decomposition.
 *    LB_GID **export_global_ids --  Pointer to array of Global IDs for the
 *                                   objects to be exported from the current
 *                                   processor.
 *    LB_LID **export_local_ids  --  Pointer to array of Local IDs for the
 *                                   objects to be exported (local to the
 *                                   current processor).
 *    int **export_procs         --  Pointer to array of Processor IDs for the
 *                                   objects to be exported (processor IDs of
 *                                   destination processors).
 *  Returned value:
 *    int                        --  Error code
 */


extern int LB_Compute_Destinations(struct LB_Struct *lb,
                                   int num_import, LB_GID *import_global_ids,
                                   LB_LID *import_local_ids, int *import_procs, 
                                   int *num_export, LB_GID **export_global_ids,
                                   LB_LID **export_local_ids,
                                   int **export_procs);

/*****************************************************************************/
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
 *
 *  Input:
 *    struct LB_Struct *lb       --  Load balancing object for current balance.
 *    int num_import             --  Number of non-local objects assigned to the
 *                                   processor in the new decomposition.
 *    LB_GID *import_global_ids  --  Array of global IDs for non-local objects
 *                                   assigned to this processor in the new
 *                                   decomposition.
 *    LB_LID *import_local_ids   --  Array of local IDs for non-local objects
 *                                   assigned to the processor in the new
 *                                   decomposition.
 *    int *import_procs          --  Array of processor IDs of processors owning
 *                                   the non-local objects that are assigned to
 *                                   this processor in the new decomposition.
 *    int num_export             --  The number of local objects that need to be
 *                                   exported from the processor to establish
 *                                   the new decomposition.
 *    LB_GID *export_global_ids  --  Array of Global IDs for the objects to be 
 *                                   exported from the current processor.
 *    LB_LID *export_local_ids   --  Array of Local IDs for the objects to be 
 *                                   exported (local to the current processor).
 *    int *export_procs          --  Array of Processor IDs for the objects to
 *                                   be exported (processor IDs of destination
 *                                   processor).
 *  Output:
 *    none                       --  The objects are migrated to their new
 *                                   processor locations.  The input arrays
 *                                   are unchanged.
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Help_Migrate(struct LB_Struct *lb,
                           int num_import, LB_GID *import_global_ids,
                           LB_LID *import_local_ids, int *import_procs,
                           int num_export, LB_GID *export_global_ids,
                           LB_LID *export_local_ids, int *export_procs);

/*****************************************************************************/
/*
 *  Routine to free the data arrays returned by LB_Balance.  The arrays
 *  are freed and the pointers are set to NULL.
 *
 *  Input:
 *    LB_GID **import_global_ids --  Pointer to array of global IDs for imported
 *                                   objects.
 *    LB_LID **import_local_ids  --  Pointer to array of local IDs for imported 
 *                                   objects.
 *    int **import_procs         --  Pointer to array of processor IDs of 
 *                                   imported objects.
 *    LB_GID **export_global_ids --  Pointer to array of global IDs for exported
 *                                   objects.
 *    LB_LID **export_local_ids  --  Pointer to array of local IDs for exported
 *                                   objects.
 *    int **export_procs         --  Pointer to array of destination processor
 *                                   IDs of exported objects.
 *  Returned value:
 *    int                        --  Error code
 */
extern int LB_Free_Data(LB_GID **import_global_ids, LB_LID **import_local_ids,
                        int **import_procs,
                        LB_GID **export_global_ids, LB_LID **export_local_ids,
                        int **export_procs);

/*****************************************************************************/
#endif
