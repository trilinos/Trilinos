
/*
 *  Data types and functions describing the interface between the
 *  application and the load balancing tool.
 */

/*
 *  Enumerated type used to indicate which function is to be set by
 *  LB_Set_LB_Fn.
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
  LB_MAX_FN_TYPES               /*  This entry should always be last.        */
};

typedef enum LB_Fn_Type LB_FN_TYPE;

/*
 *  Other common definitions:
 */

typedef int LB_OBJECT_TYPE;
typedef int LB_ID;

struct LB_Struct;

/*****************************************************************************/
/*****************************************************************************/
/**********************  Functions to query application  *********************/
/*****************************************************************************/
/*****************************************************************************/

/* 
 *  Function to return the weight of the object with a given ID and
 *  object type.
 *  Input:  
 *    LB_ID             --  the ID for the object
 *    LB_OBJECT_TYPE    --  the object type for the object
 *  Returned value:
 *    double            --  the weight for the object.
 */

typedef double LB_OBJECT_WEIGHT_FN(LB_ID, LB_OBJECT_TYPE);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID and object type,
 *  the object's number of edges (i.e., the number of objects with which
 *  the given object must communicate).
 *  Input:  
 *    LB_ID             --  the ID for the object
 *    LB_OBJECT_TYPE    --  the object type for the object
 *  Returned value:
 *    int               --  the number of neighbor objects.
 */

typedef int LB_NUM_EDGES_FN(LB_ID, LB_OBJECT_TYPE);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID and object type, 
 *  the object's edge list (i.e., objects with which the given object must
 *  communicate.
 *  Input:  
 *    LB_ID             --  the ID for the object
 *    LB_OBJECT_TYPE    --  the object type for the object
 *  Output:
 *    LB_ID *           --  Array of IDs of neighboring objects.
 */

typedef void LB_EDGE_LIST_FN(LB_ID, LB_OBJECT_TYPE, LB_ID *);

/*****************************************************************************/
/*
 *  Function to return, for a given object type,
 *  the number of geometry fields per object (e.g., the number of values
 *  used to express the coordinates of the object).
 *  Input:  
 *    LB_OBJECT_TYPE    --  the object type for the object
 *  Returned value:
 *    int               --  the number of geometry fields.
 */

typedef int LB_NUM_GEOM_FN(LB_OBJECT_TYPE);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID and object type,
 *  the geometry information for the object (e.g., coordinates).
 *  Input:  
 *    LB_ID             --  the ID for the object
 *    LB_OBJECT_TYPE    --  the object type for the object
 *  Output:
 *    double *          --  the geometry info for the object (e.g., coordinates)
 */

typedef void LB_GEOM_FN(LB_ID, LB_OBJECT_TYPE, double *);

/*****************************************************************************/
/*
 *  Function to return, for the calling processor, the number of objects 
 *  of type LB_OBJECT_TYPE located in that processor's memory.
 *  Input:  
 *    LB_OBJECT_TYPE    --  the object type to be counted.
 *  Returned value:
 *    int               --  the number of local objects of type LB_OBJECT_TYPE.
 */

typedef int LB_NUM_OBJ_FN(LB_OBJECT_TYPE);

/*****************************************************************************/
/*
 *  Function to return a list of all local objects of type LB_OBJECT_TYPE
 *  on a processor.
 *  Input:  
 *    LB_OBJECT_TYPE    --  the object type to be counted.
 *  Output:
 *    LB_ID *           --  array of IDs of all objects on the processor.
 */

typedef void LB_GET_LOCAL_OBJECTS_FN(LB_OBJECT_TYPE, LB_ID *);

/*****************************************************************************/
/*
 *  Iterator function for local objects; return the next local object of
 *  type LB_OBJECT_TYPE.
 *  Input:  
 *    LB_ID             --  ID of the previous object; NULL if requesting
 *                          first object.
 *    LB_OBJECT_TYPE    --  the object type of the object returned.
 *  Returned value:
 *    LB_ID             --  ID of the next object; NULL if no more objects.
 */

typedef LB_ID LB_NEXT_OBJ_FN(LB_ID, LB_OBJECT_TYPE);

/*****************************************************************************/
/*
 *  Function to return, for the calling processor, the number of objects 
 *  of type LB_OBJECT_TYPE sharing a subdomain border with a given processor.
 *  Input:  
 *    LB_OBJECT_TYPE    --  the object type to be counted.
 *    int               --  processor ID of the neighboring processor.
 *  Returned value:
 *    int               --  the number of local objects of type LB_OBJECT_TYPE.
 */

typedef int LB_NUM_BORDER_OBJ_FN(LB_OBJECT_TYPE, int);

/*****************************************************************************/
/*
 *  Function to return a list of all objects of type LB_OBJECT_TYPE
 *  sharing a subdomain border with a given processor.
 *  Input:  
 *    LB_OBJECT_TYPE    --  the object type to be counted.
 *    int               --  processor ID of the neighboring processor.
 *  Output:
 *    LB_ID *           --  array of IDs of all objects on the processor.
 */

typedef void LB_BORDER_OBJ_FN(LB_OBJECT_TYPE, int, LB_ID *);

/*****************************************************************************/
/*
 *  Iterator function for border objects; return the next local object of
 *  type LB_OBJECT_TYPE along the subdomain boundary with a given processor.
 *  Input:  
 *    LB_ID             --  ID of the previous object; NULL if requesting
 *                          first object.
 *    LB_OBJECT_TYPE    --  the object type of the object returned.
 *    int               --  processor ID of the neighboring processor.
 *  Returned value:
 *    LB_ID             --  ID of the next object; NULL if no more objects.
 */

typedef LB_ID LB_NEXT_BORDER_OBJ_FN(LB_ID, LB_OBJECT_TYPE, int);


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
extern struct LB_Struct *LB_Create_LB_Object();

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
extern void LB_Set_LB_Fn(struct LB_Struct *, LB_FN_TYPE, void *());

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
extern void LB_Set_LB_Method(struct LB_Struct *, char *, double *);

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
extern void LB_Set_LB_Tolerance(struct LB_Struct *, double);

/*****************************************************************************/
/*
 *  Function to set the object type for objects to be balanced.
 *  The object type is an integer value.  It can be used to help distinguish
 *  between the IDs for, say, elements and surfaces.
 *  This value is used only by the application; it is optional as far
 *  as the load-balancer is concerned.
 *  Input:
 *    struct LB_Struct * --  The load balancing object to which this tolerance
 *                           applies.
 *    int                --  An integer representing the object type.
 *  Output:
 *    struct LB_Struct * --  Appropriate fields set to designated type.
 */
extern void LB_Set_LB_Object_Type(struct LB_Struct *, int);
