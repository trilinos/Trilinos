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


/*
 *  Data types and functions describing the interface between the
 *  application and the load balancing tool.
 */

/*
 *  Enumerated type used to indicate which function is to be set by
 *  LB_Set_Fn.
 */

enum LB_Fn_Type {
  LB_OBJECT_WEIGHT_FN_TYPE = 0,
  LB_NUM_EDGES_FN_TYPE,
  LB_EDGE_LIST_FN_TYPE,
  LB_NUM_GEOM_FN_TYPE,
  LB_GEOM_FN_TYPE,
  LB_NUM_OBJ_FN_TYPE,
  LB_GET_LOCAL_OBJECTS_FN_TYPE,
  LB_NEXT_OBJ_FN_TYPE,
  LB_NUM_BORDER_OBJ_FN_TYPE,
  LB_BORDER_OBJ_FN_TYPE,
  LB_NEXT_BORDER_OBJ_FN_TYPE,
  LB_PRE_MIGRATE_FN_TYPE,
  LB_OBJECT_SIZE_FN_TYPE,
  LB_PACK_OBJECT_FN_TYPE,
  LB_UNPACK_OBJECT_FN_TYPE,
  LB_MAX_FN_TYPES               /*  This entry should always be last.        */
};

typedef enum LB_Fn_Type LB_FN_TYPE;

/*
 *  Other common definitions:
 */

typedef int LB_ID;

struct LB_Struct;

/*
 *  Maximum number of parameters to be passed to any load-balancing
 *  method.
 */

#define LB_PARAMS_MAX_SIZE 5

/*****************************************************************************/
/*****************************************************************************/
/**********************  Functions to query application  *********************/
/*****************************************************************************/
/*****************************************************************************/

/* 
 *  Function to return the weight of the object with a given ID and
 *  object type.
 *  Input:  
 *    LB_ID             --  the Global ID for the object
 *    LB_ID             --  the Local ID for the object
 *  Returned value:
 *    double            --  the weight for the object.
 */

typedef double LB_OBJECT_WEIGHT_FN(LB_ID, LB_ID);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID and object type,
 *  the object's number of edges (i.e., the number of objects with which
 *  the given object must communicate).
 *  Input:  
 *    LB_ID             --  the Global ID for the object
 *    LB_ID             --  the Local ID for the object
 *  Returned value:
 *    int               --  the number of neighbor objects.
 */

typedef int LB_NUM_EDGES_FN(LB_ID, LB_ID);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID and object type, 
 *  the object's edge list (i.e., objects with which the given object must
 *  communicate.
 *  Input:  
 *    LB_ID             --  the Global ID for the object
 *    LB_ID             --  the Local ID for the object
 *  Output:
 *    LB_ID *           --  Array of Global IDs of neighboring objects.
 *    LB_ID *           --  Array of Local IDs of neighboring objects.
 */

typedef void LB_EDGE_LIST_FN(LB_ID, LB_ID, LB_ID *, LB_ID *);

/*****************************************************************************/
/*
 *  Function to return, for a given object type,
 *  the number of geometry fields per object (e.g., the number of values
 *  used to express the coordinates of the object).
 *  Input:  
 *  Returned value:
 *    int               --  the number of geometry fields.
 */

typedef int LB_NUM_GEOM_FN();

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID and object type,
 *  the geometry information for the object (e.g., coordinates).
 *  Input:  
 *    LB_ID             --  the Global ID for the object
 *    LB_ID             --  the Local ID for the object
 *  Output:
 *    double *          --  the geometry info for the object (e.g., coordinates)
 */

typedef void LB_GEOM_FN(LB_ID, LB_ID, double *);

/*****************************************************************************/
/*
 *  Function to return, for the calling processor, the number of objects 
 *  located in that processor's memory.
 *  Input:  
 *    none
 *  Returned value:
 *    int               --  the number of local objects.
 */

typedef int LB_NUM_OBJ_FN();

/*****************************************************************************/
/*
 *  Function to return a list of all local objects on a processor.
 *  Input:  
 *    none
 *  Output:
 *    LB_ID *           --  array of Global IDs of all objects on the processor.
 *    LB_ID *           --  array of Local IDs of all objects on the processor.
 */

typedef void LB_GET_LOCAL_OBJECTS_FN(LB_ID *, LB_ID *);

/*****************************************************************************/
/*
 *  Iterator function for local objects; return the next local object.
 *  Input:  
 *    LB_ID             --  Global ID of the previous object; NULL if requesting
 *                          first object.
 *    LB_ID             --  Local ID of the previous object; NULL if requesting
 *                          first object.
 *  Output:
 *    LB_ID *           --  Global ID of the next object; NULL if no more
 *                          objects.
 *    LB_ID *           --  Local ID of the next object; NULL if no more
 *                          objects.
 */

typedef void LB_NEXT_OBJ_FN(LB_ID, LB_ID, LB_ID *, LB_ID *);

/*****************************************************************************/
/*
 *  Function to return, for the calling processor, the number of objects 
 *  sharing a subdomain border with a given processor.
 *  Input:  
 *    int               --  processor ID of the neighboring processor.
 *  Returned value:
 *    int               --  the number of local objects.
 */

typedef int LB_NUM_BORDER_OBJ_FN(int);

/*****************************************************************************/
/*
 *  Function to return a list of all objects sharing a subdomain border 
 *  with a given processor.
 *  Input:  
 *    int               --  processor ID of the neighboring processor.
 *  Output:
 *    LB_ID *           --  array of Global IDs of all objects on the processor.
 *    LB_ID *           --  array of Local IDs of all objects on the processor.
 */

typedef void LB_BORDER_OBJ_FN(int, LB_ID *, LB_ID *);

/*****************************************************************************/
/*
 *  Iterator function for border objects; return the next local object 
 *  along the subdomain boundary with a given processor.
 *  Input:  
 *    LB_ID             --  Global ID of the previous object; NULL if requesting
 *                          first object.
 *    LB_ID             --  Local ID of the previous object; NULL if requesting
 *                          first object.
 *    int               --  processor ID of the neighboring processor.
 *  Output:
 *    LB_ID *           --  Global ID of the next object; NULL if no more
 *                          objects.
 *    LB_ID *           --  Local ID of the next object; NULL if no more
 *                          objects.
 */

typedef LB_ID LB_NEXT_BORDER_OBJ_FN(LB_ID, LB_ID, int, LB_ID *, LB_ID *);

/*****************************************************************************/
/*
 *  Function to return the size (in bytes) of data to be migrated for the
 *  given object type.  This function is needed only when the application
 *  wants the load-balancer to help migrate the data.  It is used by the
 *  comm.c routines to allocate message buffers.
 *  Input:  
 *    none
 *  Returned value:
 *    int               --  the size of data of local objects.
 */

typedef int LB_OBJECT_SIZE_FN();

/*****************************************************************************/
/*
 *  Function called as a pre-processor to the migration.  This function is 
 *  optional, and is used only when the application wants the load-balancer 
 *  to help migrate the data.  The application can perform any type of 
 *  pre-processing in this function.
 *  Input:  
 *    int               --  Number of objects to be imported.
 *    LB_ID *           --  Global IDs of objects to be imported.
 *    LB_ID *           --  Local IDs of objects to be imported.
 *    int *             --  Processor IDs of importing processors.
 *    int               --  Number of objects to be exported.
 *    LB_ID *           --  Global IDs of objects to be exported.
 *    LB_ID *           --  Local IDs of objects to be exported.
 *    int *             --  Processor IDs of processors to receive the objects.
 *  Output:
 *    none              --  the application performs pre-processing.
 */

typedef void LB_PRE_MIGRATE_FN(int, LB_ID *, LB_ID *, int *,
                               int, LB_ID *, LB_ID *, int *);

/*****************************************************************************/
/*
 *  Function to pack data to be migrated for the given object and object type.
 *  This function is needed only when the application wants the load-balancer 
 *  to help migrate the data.  It packs all data related to the given object
 *  into a communication buffer, the starting address of which is provided
 *  by the load-balancer.
 *  Input:  
 *    LB_ID             --  Global ID of the object to be packed.
 *    LB_ID             --  Local ID of the object to be packed.
 *    int               --  Processor ID of the destination processor.
 *    int               --  number of bytes allowed for the object to be packed.
 *    char *            --  starting address of buffer into which to pack the
 *                          object.
 *  Output:
 *    char *            --  the buffer is rewritten with the packed data.
 */

typedef void LB_PACK_OBJECT_FN(LB_ID, LB_ID, int, int, char *);

/*****************************************************************************/
/*
 *  Function to unpack data for an object migrated to a new processor.
 *  This function is needed only when the application wants the load-balancer 
 *  to help migrate the data.  The data is stored in a buffer (char *); the
 *  size of the data for the object is included.
 *  Input:  
 *    LB_ID             --  Global ID of the object to be unpacked.
 *    int               --  number of bytes in the buffer for the object.
 *    char *            --  starting address of buffer into which to pack the
 *                          object.
 *  Output:
 *    none              --  the routine processors the data in the buffer.
 */

typedef void LB_UNPACK_OBJECT_FN(LB_ID, int, char *);

/*****************************************************************************/
/*****************************************************************************/
/**********************  Functions to set-up LB object ***********************/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Function to initialize values needed in load balancing tools.
 *  This function must be called with the argc, argv arguments from main.
 *  If the application uses MPI, call this function after calling MPI_Init.
 *  If the application does not use MPI, this function calls MPI_Init for
 *  use by the load balancer.
 */
extern void LB_Initialize(int argc, char **argv);

/*****************************************************************************/
/*
 *  Function to create a load balancing object.  May want more than one
 *  object if using different decompositions with different techniques.
 *  This function allocates and initializes the object.
 *  Returned value:
 *    struct LB_Struct * --  Pointer to a LB object.
 */
extern struct LB_Struct *LB_Create_Object();

/*****************************************************************************/
/*
 *  Function to initialize a given LB interface function.
 *  Input:
 *    struct LB_Struct *--  Pointer to a LB object.
 *    LB_FN_TYPE        --  Enum type indicating the function to be set.
 *    void *()          --  Pointer to the function to be used in the 
 *                          assignment.
 *  Output:
 *    void *            --  Appropriate field set to value in void *().
 */
extern void LB_Set_Fn(struct LB_Struct *, LB_FN_TYPE, void *());

/*****************************************************************************/
/*
 *  Function to set the load balancing method to be used.
 *  Input:
 *    struct LB_Struct * --  The load balancing object to which this method
 *                           applies.
 *    char *             --  String specifying the desired method.
 *    double *           --  Params needed by desired method.  (This field
 *                           will be better defined later.)
 *  Output:
 *    struct LB_Struct * --  Appropriate fields set to designated values.
 */
extern void LB_Set_Method(struct LB_Struct *, char *, double *);

/*****************************************************************************/
/*
 *  Function to set the tolerance to which the system must be load balanced.
 *  For example, if the tolerance is set to 0.9, 10% load imbalance between
 *  the most heavily loaded processor and the average load will be accepted
 *  as balanced.
 *  Input:
 *    struct LB_Struct * --  The load balancing object to which this tolerance
 *                           applies.
 *    double             --  The tolerance desired.
 *  Output:
 *    struct LB_Struct * --  Appropriate fields set to designated value.
 */
extern void LB_Set_Tolerance(struct LB_Struct *, double);

/*****************************************************************************/
/*
 *  Function to set a flag indicating whether the application wants the
 *  load-balancer to help with data migration.   If migration help is
 *  wanted, routines to pack and unpack object data must be provided by
 *  the application (see LB_OBJECT_SIZE_FN, LB_PACK_OBJECT_FN, 
 *  LB_UNPACK_OBJECT_FN).
 *
 *  Input:
 *    struct LB_Struct * --  The load balancing object to which this tolerance
 *                           applies.
 *    int                --  TRUE or FALSE to indicate whether the application
 *                           wants migration help.  Default is FALSE.
 *  Output:
 *    struct LB_Struct * --  Appropriate fields set to designated type.
 */
extern void LB_Set_Migration(struct LB_Struct *, int);

/*****************************************************************************/
/*
 *  Function to initialize an array to pass parameters to the load-balancing
 *  methods.  This function is provided so that the load-balancer 
 *  look for array entries not set by the application and use default values
 *  for those entries.
 *
 *  Input/Output:
 *    double *           --  Pointer to the array to be used to pass
 *                           parameters to the load-balancing methods.
 *                           Upon return, the values in this array are
 *                           initialized to an initial value determined
 *                           by the load-balancer.
 */

extern void LB_Initialize_Params_Array(double *);

/*****************************************************************************/
/*
 *  Function to invoke the load-balancer.
 *
 *  Input:
 *    struct LB_Struct * -- The load balancing object containing info about
 *                           this load-balancing invocation.
 *  Output:
 *    int *              --  The number of non-local objects in the 
 *                           processor's new decomposition (i.e., number of
 *                           objects to be imported).
 *    LB_ID **           --  Pointer to array of Global IDs for the objects to
 *                           be imported.
 *    LB_ID **           --  Pointer to array of Local IDs for the objects to
 *                           be imported (local to the exporting processor).
 *    int **             --  Pointer to array of Processor IDs for the objects
 *                           to be imported (processor IDs of source processor).
 *    int *              --  The number of local objects that need to be 
 *                           exported from the processor to establish the
 *                           new decomposition.
 *    LB_ID **           --  Pointer to array of Global IDs for the objects to
 *                           be exported from the current processor.
 *    LB_ID **           --  Pointer to array of Local IDs for the objects to 
 *                           be exported (local to the current processor).
 *    int **             --  Pointer to array of Processor IDs for the objects
 *                           to be exported (processor IDs of destination 
 *                           processors).
 *  
 *  Returned value: zero --  No changes to the decomposition were made by the
 *                           load-balancing algorithm; migration is not needed.
 *                  one  --  A new decomposition is suggested by the
 *                           load-balancer; migration is needed to establish
 *                           the new decomposition.
 */

extern int LB_Balance(struct LB_Struct *, int *, LB_ID **, LB_ID **, int**,
                      int *, LB_ID **, LB_ID **, int **);

/*****************************************************************************/
/*
 *  Routine to compute the inverse map:  Given, for each processor, a list
 *  of non-local objects assigned to the processor, compute the list of objects
 *  that processor needs to export to other processors to establish the new
 *  decomposition.
 *
 *  Input:
 *    struct LB_Struct * --  Load balancing object for current balance.
 *    int                --  Number of non-local objects assigned to the
 *                           processor in the new decomposition.
 *    LB_ID *            --  Array of global IDs for non-local objects
 *                           assigned to this processor in the new
 *                           decomposition.
 *    LB_ID *            --  Array of local IDs for non-local objects
 *                           assigned to the processor in the new
 *                           decomposition.
 *    int *              --  Array of processor IDs of processors owning
 *                           the non-local objects that are assigned to
 *                           this processor in the new decomposition.
 *  Output:
 *    int *              --  The number of local objects that need to be 
 *                           exported from the processor to establish the
 *                           new decomposition.
 *    LB_ID **           --  Pointer to array of Global IDs for the objects 
 *                           to be exported from the current processor.
 *    LB_ID **           --  Pointer to array of Local IDs for the objects 
 *                           to be exported (local to the current processor).
 *    int **             --  Pointer to array of Processor IDs for the objects
 *                           to be exported (processor IDs of destination
 *                           processors).
 */


extern void LB_Compute_Destinations(struct LB_Struct *,
                                    int, LB_ID *, LB_ID *, int *, 
                                    int *, LB_ID **, LB_ID **, int **);

/*****************************************************************************/
/*
 *  Routine to help perform migration.  If migration pre-processing routine
 *  (LB_PRE_MIGRATE_FN) is specified, this routine first calls that function.
 *  It then calls a function to obtain the size of the migrating objects
 *  (LB_OBJECT_SIZE_FN).  The routine next calls an application-specified
 *  object packing routine (LB_PACK_OBJECT_FN) for each object
 *  to be exported.  It develops the needed communication map to move the
 *  objects to other processors.  It performs the communication according
 *  to the map, and then calls an application-specified object unpacking 
 *  routine (LB_UNPACK_OBJECT_FN) for each object imported.
 *
 *  Input:
 *    struct LB_Struct * --  Load balancing object for current balance.
 *    int                --  Number of non-local objects assigned to the
 *                           processor in the new decomposition.
 *    LB_ID *            --  Array of global IDs for non-local objects
 *                           assigned to this processor in the new
 *                           decomposition.
 *    LB_ID *            --  Array of local IDs for non-local objects
 *                           assigned to the processor in the new
 *                           decomposition.
 *    int *              --  Array of processor IDs of processors owning
 *                           the non-local objects that are assigned to
 *                           this processor in the new decomposition.
 *    int                --  The number of local objects that need to be 
 *                           exported from the processor to establish the
 *                           new decomposition.
 *    LB_ID *            --  Array of Global IDs for the objects to be 
 *                           exported from the current processor.
 *    LB_ID *            --  Array of Local IDs for the objects to be 
 *                           exported (local to the current processor).
 *    int *              --  Array of Processor IDs for the objects to be 
 *                           exported (processor IDs of destination processor).
 *  Output:
 *    none               --  The objects are migrated to their new processor
 *                           locations.  The input arrays are unchanged.
 */

extern void LB_Help_Migrate(struct LB_Struct *,
                            int, LB_ID *, LB_ID *, int *,
                            int, LB_ID *, LB_ID *, int *);

/*****************************************************************************/
/*
 *  Routine to free the data arrays returned by LB_Balance.  The arrays
 *  are freed and the pointers are set to NULL.
 *
 *  Input:
 *    LB_ID **           --  Pointer to array of global IDs for imported
 *                           objects.
 *    LB_ID **           --  Pointer to array of local IDs for imported objects.
 *    int **             --  Pointer to array of processor IDs of imported
 *                           objects.
 *    LB_ID **           --  Pointer to array of global IDs for exported
 *                           objects.
 *    LB_ID **           --  Pointer to array of local IDs for exported objects.
 *    int **             --  Pointer to array of destination processor IDs of
 *                           exported objects.
 */
extern void LB_Free_Data(LB_ID **, LB_ID **, int **,
                         LB_ID **, LB_ID **, int **);

/*****************************************************************************/
#endif
