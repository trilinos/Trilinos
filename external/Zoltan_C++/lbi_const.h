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
 *    Revision: 1.57 $
 ****************************************************************************/

#ifndef __LBI_CONST_H
#define __LBI_CONST_H


#ifdef __cplusplus
#define MPI_NO_CPPBIND
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <mpi.h>

/*
 *  Data type LB_ID for global and local identifiers used in Zoltan.
 */

typedef unsigned int    LB_ID_TYPE;
typedef LB_ID_TYPE     *LB_ID_PTR;
#define LB_ID_MPI_TYPE  MPI_UNSIGNED

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
  LB_MID_MIGRATE_FN_TYPE,
  LB_POST_MIGRATE_FN_TYPE,
  LB_OBJ_SIZE_FN_TYPE,
  LB_PACK_OBJ_FN_TYPE,
  LB_UNPACK_OBJ_FN_TYPE,
  LB_NUM_COARSE_OBJ_FN_TYPE,
  LB_COARSE_OBJ_LIST_FN_TYPE,
  LB_FIRST_COARSE_OBJ_FN_TYPE,
  LB_NEXT_COARSE_OBJ_FN_TYPE,
  LB_NUM_CHILD_FN_TYPE,
  LB_CHILD_LIST_FN_TYPE,
  LB_CHILD_WEIGHT_FN_TYPE,
  LB_GET_PROCESSOR_NAME_FN_TYPE,
  LB_MAX_FN_TYPES               /*  This entry should always be last.        */
};

typedef enum LB_Fn_Type LB_FN_TYPE;

/*
 * Enumerated type used to indicate what type of refinement was used when
 * building a refinement tree.
 */

enum LB_Ref_Type {
  LB_OTHER_REF,      /* unspecified type of refinement */
  LB_IN_ORDER,       /* user provides the order of the children */
  LB_TRI_BISECT,     /* bisection of triangles */
  LB_QUAD_QUAD,      /* quadrasection of quadralaterals */
  LB_HEX3D_OCT       /* octasection of hexahedra */
};

typedef enum LB_Ref_Type LB_REF_TYPE;

/*
 *  Other common definitions:
 */

struct LB_Struct;

/*
 * Error codes for Zoltan library
 *   LB_OK     - no errors
 *   LB_WARN   - some warning occurred in Zoltan library; application should be
 *               able to continue running
 *   LB_FATAL  - a fatal error occurred
 *   LB_MEMERR - memory allocation failed; with this error, it could be
 *               possible to try a different, more memory-friendly, algorithm
 */
#define LB_OK     0
#define LB_WARN   1
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
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  the Global ID for the object
 *    LB_ID_PTR local_id        --  the Local ID for the object
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the number of neighbor objects.
 */

typedef int LB_NUM_EDGES_FN(void *data, 
                            int num_gid_entries, int num_lid_entries,
                            LB_ID_PTR global_id, LB_ID_PTR local_id,
                            int *ierr);

typedef int LB_NUM_EDGES_FORT_FN(void *data, 
                                 int *num_gid_entries, int *num_lid_entries,
                                 LB_ID_PTR global_id, LB_ID_PTR local_id, 
                                 int *ierr);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID, 
 *  the object's edge list (i.e., objects with which the given object must
 *  communicate.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  the Global ID for the object
 *    LB_ID_PTR local_id        --  the Local ID for the object
 *    int    wdim               --  dimension of edge weights, or 0 if
 *                                  edge weights are not sought.
 *  Output:
 *    LB_ID_PTR nbor_global_ids --  Array of Global IDs of neighboring objects.
 *    int    *nbor_procs        --  Array of neighboring procs.
 *    int    *nbor_ewgts        --  Array of edge weights, where 
 *                                  nbor_ewgts[i*wdim:(i+1)*wdim-1]
 *                                  corresponds to the weight of edge i
 *    int *ierr                 --  error code
 */

typedef void LB_EDGE_LIST_FN(void *data, 
                             int num_gid_entries, int num_lid_entries,
                             LB_ID_PTR global_id, LB_ID_PTR local_id,
                             LB_ID_PTR nbor_global_id, int *nbor_procs,
                             int wdim, float *nbor_ewgts, int *ierr);

typedef void LB_EDGE_LIST_FORT_FN(void *data, 
                                  int *num_gid_entries, int *num_lid_entries,
                                  LB_ID_PTR global_id, LB_ID_PTR local_id, 
                                  LB_ID_PTR nbor_global_id,
                                  int *nbor_procs, int *wdim, 
                                  float *nbor_ewgts, int *ierr);

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

typedef int LB_NUM_GEOM_FORT_FN(void *data, int *ierr);

/*****************************************************************************/
/*
 *  Function to return, for the object with a given ID,
 *  the geometry information for the object (e.g., coordinates).
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  the Global ID for the object
 *    LB_ID_PTR local_id        --  the Local ID for the object
 *  Output:
 *    double *geom_vec          --  the geometry info for the object
 *                                  (e.g., coordinates)
 *    int *ierr                 --  error code
 */

typedef void LB_GEOM_FN(void *data, int num_gid_entries, int num_lid_entries,
                        LB_ID_PTR global_id, LB_ID_PTR local_id,
                        double *geom_vec, int *ierr);

typedef void LB_GEOM_FORT_FN(void *data, 
                             int *num_gid_entries, int *num_lid_entries,
                             LB_ID_PTR global_id, LB_ID_PTR local_id,
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

typedef int LB_NUM_OBJ_FORT_FN(void *data, int *ierr);

/*****************************************************************************/
/*
 *  Function to return a list of all local objects on a processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    int wdim                  --  dimension of object weights, or 0 if
 *                                  object weights are not sought. 
 *  Output:
 *    LB_ID_PTR global_ids      --  array of Global IDs of all objects on the
 *                                  processor.
 *    LB_ID_PTR local_ids       --  array of Local IDs of all objects on the
 *                                  processor.
 *    float *objwgts            --  objwgts[i*wdim:(i+1)*wdim-1] correponds
 *                                  to the weight of object i 
 *    int *ierr                 --  error code
 */

typedef void LB_OBJ_LIST_FN(void *data, 
                            int num_gid_entries, int num_lid_entries,
                            LB_ID_PTR global_ids, LB_ID_PTR local_ids,
                            int wdim, float *objwgts, int *ierr);

typedef void LB_OBJ_LIST_FORT_FN(void *data, 
                                 int *num_gid_entries, int *num_lid_entries,
                                 LB_ID_PTR global_ids, LB_ID_PTR local_ids,
                                 int *wdim, float *objwgts,
                                 int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for local objects; return the first local object on
 *  the processor.  This function should be used with LB_NEXT_OBJ_FN.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *  Output:
 *    LB_ID_PTR first_global_id --  Global ID of the first object; NULL if no
 *                                  objects.
 *    LB_ID_PTR first_local_id  --  Local ID of the first object; NULL if no
 *                                  objects.
 *    float *first_obj_wgt      --  weight vector for first object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist on the processor.
 */

typedef int LB_FIRST_OBJ_FN(void *data, 
                            int num_gid_entries, int num_lid_entries,
                            LB_ID_PTR first_global_id,
                            LB_ID_PTR first_local_id, 
                            int wdim, float *first_obj_wgt, int *ierr);

typedef int LB_FIRST_OBJ_FORT_FN(void *data, 
                                 int *num_gid_entries, int *num_lid_entries,
                                 LB_ID_PTR first_global_id,
                                 LB_ID_PTR first_local_id, int *wdim,
                                 float *first_obj_wgt, int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for local objects; return the next local object.
 *  This function should be used with LB_FIRST_OBJ_FN.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  Global ID of the previous object.
 *    LB_ID_PTR local_id        --  Local ID of the previous object.
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *  Output:
 *    LB_ID_PTR next_global_id  --  Global ID of the next object; NULL if no
 *                                  more objects.
 *    LB_ID_PTR next_local_id   --  Local ID of the next object; NULL if no
 *                                  more objects.
 *    float *next_obj_wgt       --  weight vector for the next object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist (i.e., global_id is
 *                                  the last object).
 */

typedef int LB_NEXT_OBJ_FN(void *data, int num_gid_entries, int num_lid_entries,
                           LB_ID_PTR global_id, LB_ID_PTR local_id,
                           LB_ID_PTR next_global_id, LB_ID_PTR next_local_id,
                           int wdim, float *next_obj_wgt, int *ierr);

typedef int LB_NEXT_OBJ_FORT_FN(void *data, 
                                int *num_gid_entries, int *num_lid_entries,
                                LB_ID_PTR global_id, LB_ID_PTR local_id,
                                LB_ID_PTR next_global_id,
                                LB_ID_PTR next_local_id,
                                int *wdim, float *next_obj_wgt, int *ierr);

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

typedef int LB_NUM_BORDER_OBJ_FORT_FN(void *data, int *nbor_proc, int *ierr);

/*****************************************************************************/
/*
 *  Function to return a list of all objects sharing a subdomain border 
 *  with a given processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    int nbor_proc             --  processor ID of the neighboring processor.
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weights are not sought.
 *  Output:
 *    LB_ID_PTR global_ids      --  array of Global IDs of all objects on the
 *                                  processor border with the given neighboring
 *                                  processor.
 *    LB_ID_PTR local_ids       --  array of Local IDs of all objects on the 
 *                                  processor border with the given neighboring 
 *                                  processor.
 *    float *objwgts            --  objwgts[i*wdim:(i+1)*wdim-1] correponds
 *                                  to the weight of object i 
 *                                  (objwgts is undefined if wdim=0)
 *    int *ierr                 --  error code
 */

typedef void LB_BORDER_OBJ_LIST_FN(void *data, 
                                   int num_gid_entries, int num_lid_entries,
                                   int nbor_proc,
                                   LB_ID_PTR global_ids, LB_ID_PTR local_ids,
                                   int wdim, float *objwgts, int *ierr);

typedef void LB_BORDER_OBJ_LIST_FORT_FN(void *data,
                                     int *num_gid_entries, int *num_lid_entries,
                                     int *nbor_proc,
                                     LB_ID_PTR global_ids, LB_ID_PTR local_ids,
                                     int *wdim, float *objwgts, int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for border objects; return the first local object 
 *  along the subdomain boundary with a given processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    int nbor_proc             --  processor ID of the neighboring processor.
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *  Output:
 *    LB_ID_PTR first_global_id --  Global ID of the first object; NULL if no
 *                                  objects.
 *    LB_ID_PTR first_local_id  --  Local ID of the first object; NULL if no 
 *                                  objects.
 *    float *first_obj_wgt      --  weight vector for the first object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist (i.e., global_id is
 *                                  the last object).
 */

typedef int LB_FIRST_BORDER_OBJ_FN(void *data, 
                                   int num_gid_entries, int num_lid_entries, 
                                   int nbor_proc,
                                   LB_ID_PTR first_global_id,
                                   LB_ID_PTR first_local_id, 
                                   int wdim, float *first_obj_wgt,
                                   int *ierr);

typedef int LB_FIRST_BORDER_OBJ_FORT_FN(void *data, 
                                    int *num_gid_entries, int *num_lid_entries, 
                                    int *nbor_proc,
                                    LB_ID_PTR first_global_id,
                                    LB_ID_PTR first_local_id, 
                                    int *wdim, float *first_obj_wgt,
                                    int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for border objects; return the next local object 
 *  along the subdomain boundary with a given processor.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  Global ID of the previous object.
 *    LB_ID_PTR local_id        --  Local ID of the previous object.
 *    int nbor_proc             --  processor ID of the neighboring processor.
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *  Output:
 *    LB_ID_PTR next_global_id  --  Global ID of the next object; NULL if no
 *                                  more objects.
 *    LB_ID_PTR next_local_id   --  Local ID of the next object; NULL if no 
 *                                  more objects.
 *    float *next_obj_wgt       --  weight vector for the next object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist (i.e., global_id is
 *                                  the last object).
 */

typedef int LB_NEXT_BORDER_OBJ_FN(void *data, 
                                  int num_gid_entries, int num_lid_entries,
                                  LB_ID_PTR global_id,
                                  LB_ID_PTR local_id, int nbor_proc,
                                  LB_ID_PTR next_global_id,
                                  LB_ID_PTR next_local_id, 
                                  int wdim, float *next_obj_wgt,
                                  int *ierr);

typedef int LB_NEXT_BORDER_OBJ_FORT_FN(void *data, 
                                     int *num_gid_entries, int *num_lid_entries,
                                     LB_ID_PTR global_id,
                                     LB_ID_PTR local_id, int *nbor_proc,
                                     LB_ID_PTR next_global_id,
                                     LB_ID_PTR next_local_id, 
                                     int *wdim, float *next_obj_wgt,
                                     int *ierr);

/*****************************************************************************/
/*
 *  Function to return the size (in bytes) of data to be migrated.
 *  This function is needed only when the application
 *  wants the load-balancer to help migrate the data.  It is used by the
 *  comm.c routines to allocate message buffers.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  the Global ID for the object
 *    LB_ID_PTR local_id        --  the Local ID for the object
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the size of data of the object
 *                                  corresponding to global_id
 */

typedef int LB_OBJ_SIZE_FN(void *data, int num_gid_entries,
            int num_lid_entries, LB_ID_PTR global_id, 
            LB_ID_PTR local_id, int *ierr);

typedef int LB_OBJ_SIZE_FORT_FN(void *data, int *num_gid_entries,
            int *num_lid_entries, LB_ID_PTR global_id, 
            LB_ID_PTR local_id, int *ierr);

/*****************************************************************************/
/*
 *  Function called as a pre-processor to the migration.  This function is 
 *  optional, and is used only when the application wants the load-balancer 
 *  to help migrate the data.  The application can perform any type of 
 *  pre-processing in this function.
 *  Input:  
 *    void *data                  --  pointer to user defined data structure
 *    int num_gid_entries         --  number of array entries of type LB_ID_TYPE
 *                                    in a global ID
 *    int num_lid_entries         --  number of array entries of type LB_ID_TYPE
 *                                    in a local ID
 *    int num_import              --  Number of objects to be imported.
 *    LB_ID_PTR import_global_ids --  Global IDs of objects to be imported.
 *    LB_ID_PTR import_local_ids  --  Local IDs of objects to be imported.
 *    int *import_procs           --  Processor IDs of importing processors.
 *    int num_export              --  Number of objects to be exported.
 *    LB_ID_PTR export_global_ids --  Global IDs of objects to be exported.
 *    LB_ID_PTR export_local_ids  --  Local IDs of objects to be exported.
 *    int *export_procs           --  Processor IDs of processors to receive
 *                                    the objects.
 *  Output:
 *    int *ierr                   --  error code
 */

typedef void LB_PRE_MIGRATE_FN(void *data, 
                               int num_gid_entries, int num_lid_entries,
                               int num_import,
                               LB_ID_PTR import_global_ids,
                               LB_ID_PTR import_local_ids, int *import_procs,
                               int num_export, LB_ID_PTR export_global_ids,
                               LB_ID_PTR export_local_ids, int *export_procs,
                               int *ierr);

typedef void LB_PRE_MIGRATE_FORT_FN(void *data, 
                                    int *num_gid_entries, int *num_lid_entries,
                                    int *num_import,
                                    LB_ID_PTR import_global_ids,
                                    LB_ID_PTR import_local_ids, 
                                    int *import_procs,
                                    int *num_export, 
                                    LB_ID_PTR export_global_ids,
                                    LB_ID_PTR export_local_ids, 
                                    int *export_procs,
                                    int *ierr);

/*****************************************************************************/
/*
 *  Function called between the packing and unpacking phases of data migration.
 *  Within LB_Help_Migrate, the data to be migrated is packed and communicated;
 *  then this function is called (if specified). This function is 
 *  optional, and is used only when the application wants the load-balancer 
 *  to help migrate the data.  The application can perform any type of 
 *  processing in this function.
 *  Input:  
 *    void *data                  --  pointer to user defined data structure
 *    int num_gid_entries         --  number of array entries of type LB_ID_TYPE
 *                                    in a global ID
 *    int num_lid_entries         --  number of array entries of type LB_ID_TYPE
 *                                    in a local ID
 *    int num_import              --  Number of objects to be imported.
 *    LB_ID_PTR import_global_ids --  Global IDs of objects to be imported.
 *    LB_ID_PTR import_local_ids  --  Local IDs of objects to be imported.
 *    int *import_procs           --  Processor IDs of importing processors.
 *    int num_export              --  Number of objects to be exported.
 *    LB_ID_PTR export_global_ids --  Global IDs of objects to be exported.
 *    LB_ID_PTR export_local_ids  --  Local IDs of objects to be exported.
 *    int *export_procs           --  Processor IDs of processors to receive
 *                                    the objects.
 *  Output:
 *    int *ierr                   --  error code
 */

typedef void LB_MID_MIGRATE_FN(void *data, 
                               int num_gid_entries, int num_lid_entries,
                               int num_import,
                               LB_ID_PTR import_global_ids,
                               LB_ID_PTR import_local_ids, int *import_procs,
                               int num_export, LB_ID_PTR export_global_ids,
                               LB_ID_PTR export_local_ids, int *export_procs,
                               int *ierr);

typedef void LB_MID_MIGRATE_FORT_FN(void *data, 
                                    int *num_gid_entries, int *num_lid_entries,
                                    int *num_import,
                                    LB_ID_PTR import_global_ids,
                                    LB_ID_PTR import_local_ids, 
                                    int *import_procs,
                                    int *num_export, 
                                    LB_ID_PTR export_global_ids,
                                    LB_ID_PTR export_local_ids, 
                                    int *export_procs,
                                    int *ierr);

/*****************************************************************************/
/*
 *  Function called as a post-processor to the migration.  This function is 
 *  optional, and is used only when the application wants the load-balancer 
 *  to help migrate the data.  The application can perform any type of 
 *  post-processing in this function.
 *  Input:  
 *    void *data                  --  pointer to user defined data structure
 *    int num_gid_entries         --  number of array entries of type LB_ID_TYPE
 *                                    in a global ID
 *    int num_lid_entries         --  number of array entries of type LB_ID_TYPE
 *                                    in a local ID
 *    int num_import              --  Number of objects to be imported.
 *    LB_ID_PTR import_global_ids --  Global IDs of objects to be imported.
 *    LB_ID_PTR import_local_ids  --  Local IDs of objects to be imported.
 *    int *import_procs           --  Processor IDs of importing processors.
 *    int num_export              --  Number of objects to be exported.
 *    LB_ID_PTR export_global_ids --  Global IDs of objects to be exported.
 *    LB_ID_PTR export_local_ids  --  Local IDs of objects to be exported.
 *    int *export_procs           --  Processor IDs of processors to receive
 *                                    the objects.
 *  Output:
 *    int *ierr                   --  error code
 */

typedef void LB_POST_MIGRATE_FN(void *data, 
                                int num_gid_entries, int num_lid_entries,
                                int num_import,
                                LB_ID_PTR import_global_ids,
                                LB_ID_PTR import_local_ids, int *import_procs,
                                int num_export, LB_ID_PTR export_global_ids,
                                LB_ID_PTR export_local_ids, int *export_procs,
                                int *ierr);

typedef void LB_POST_MIGRATE_FORT_FN(void *data, 
                                     int *num_gid_entries, int *num_lid_entries,
                                     int *num_import,
                                     LB_ID_PTR import_global_ids,
                                     LB_ID_PTR import_local_ids, 
                                     int *import_procs,
                                     int *num_export, 
                                     LB_ID_PTR export_global_ids,
                                     LB_ID_PTR export_local_ids, 
                                     int *export_procs,
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
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  Global ID of the object to be packed.
 *    LB_ID_PTR local_id        --  Local ID of the object to be packed.
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

typedef void LB_PACK_OBJ_FN(void *data, 
                            int num_gid_entries, int num_lid_entries,
                            LB_ID_PTR global_id, LB_ID_PTR local_id,
                            int dest_proc, int size, char *buf, int *ierr);

typedef void LB_PACK_OBJ_FORT_FN(void *data, 
                                 int *num_gid_entries, int *num_lid_entries,
                                 LB_ID_PTR global_id, LB_ID_PTR local_id,
                                 int *dest_proc, int *size,
                                 char *buf, int *ierr);

/*****************************************************************************/
/*
 *  Function to unpack data for an object migrated to a new processor.
 *  This function is needed only when the application wants the load-balancer 
 *  to help migrate the data.  The data is stored in a buffer (char *); the
 *  size of the data for the object is included.
 *  Input:  
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    LB_ID_PTR global_id       --  Global ID of the object to be unpacked.
 *    int size                  --  number of bytes in the buffer for the
 *                                  object.
 *    char *buf                 --  starting address of buffer into which to
 *                                  pack the object.
 *  Output:
 *    int *ierr                 --  error code
 */

typedef void LB_UNPACK_OBJ_FN(void *data, 
                              int num_gid_entries, LB_ID_PTR global_id, 
                              int size, char *buf, int *ierr);

typedef void LB_UNPACK_OBJ_FORT_FN(void *data, int *num_gid_entries, 
                                   LB_ID_PTR global_id, int *size,
                                   char *buf, int *ierr);

/*****************************************************************************/
/*  Function to get the name of the physical processor on which
 *  the current process is running. 
 *
 *  Input:  
 *    void *data                -- pointer to user defined data structure
 *
 *  Output:
 *    char *name                -- name of the processor 
 *    int  *length              -- length of the name
 *    int  *ierr                -- error code
 */

typedef void LB_GET_PROCESSOR_NAME_FN(void *data, char *name, int *length, 
                                      int *ierr);


/*****************************************************************************/
/*
 *  Function to return the number of objects (elements) in the initial coarse
 * grid; used for initialization of the refinement tree.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the number of initial elements.
 */

typedef int LB_NUM_COARSE_OBJ_FN(void *data, int *ierr);

typedef int LB_NUM_COARSE_OBJ_FORT_FN(void *data, int *ierr);

/*****************************************************************************/
/*
 *  Function to return a list of all objects (elements) in the initial coarse
 *  grid.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *  Output:
 *    LB_ID_PTR global_ids      --  array of Global IDs of all coarse objects
 *    LB_ID_PTR local_ids       --  array of Local IDs of all coarse objects
 *    int *assigned             --  array indicating processor assignment.
 *                                  1 if the object is currently
 *                                  assigned to this processor; 0 otherwise.
 *                                  For elements that have been refined, it
 *                                  is ignored.
 *    int *num_vert             --  array containing the number of vertices
 *                                  for each object
 *    int *vertices             --  array containing the vertices for each 
 *                                  object.  If the sum of the number of
 *                                  vertices for objects 0 through i-1 is N,
 *                                  then the vertices for object i are in
 *                                  vertices[N:N+num_vert[i]]
 *    int *in_order             --  1 if the user is providing the objects in
 *                                    the order in which they should be used
 *                                  0 if the order should be determined
 *                                    automatically
 *    int *in_vertex            --  array containing the "in" vertex for each
 *                                  object, if the user provides them.  It is
 *                                  ignored if in_order==0.  For any with the
 *                                  value -1, a vertex will be selected
 *                                  automatically
 *    int *out_vertex           --  array containing the "out" vertex for each
 *                                  object; same provisions as in_vertex
 *    int *ierr                 --  error code
 */

typedef void LB_COARSE_OBJ_LIST_FN(void *data, 
                                   int num_gid_entries, int num_lid_entries,
                                   LB_ID_PTR global_ids, LB_ID_PTR local_ids,
                                   int *assigned,
                                   int *num_vert, int *vertices,
                                   int *in_order, int *in_vertex,
                                   int *out_vertex, int *ierr);

typedef void LB_COARSE_OBJ_LIST_FORT_FN(void *data,
                                    int *num_gid_entries, int *num_lid_entries,
                                    LB_ID_PTR global_ids, 
                                    LB_ID_PTR local_ids,
                                    int *assigned,
                                    int *num_vert, int *vertices,
                                    int *in_order, int *in_vertex,
                                    int *out_vertex, int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for coarse objects; return the first coarse object.
 *  This function should be used with LB_NEXT_COARSE_OBJ_FN.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *  Output:
 *    LB_ID_PTR global_id       --  Global ID of the first coarse object
 *    LB_ID_PTR local_id        --  Local ID of the first coarse object
 *    int *assigned             --  indicates processor assignment.
 *                                  1 if the object is currently
 *                                  assigned to this processor; 0 otherwise.
 *                                  For elements that have been refined, it
 *                                  is ignored.
 *    int *num_vert             --  number of vertices in the first object
 *    int *vertices             --  array containing the vertices of the first
                                    coarse object
 *    int *in_order             --  1 if the user will be providing the elements
 *                                    in the order in which they should be used
 *                                  0 if the order should be determined
 *                                    automatically
 *    int *in_vertex            --  the "in" vertex of the first coarse object.
 *                                  It is ignored if in_order==0.  If the
 *                                  value is -1, a vertex will be selected
 *                                  automatically
 *    int *out_vertex           --  array containing the "out" vertex for the
 *                                  first object; same provisions as in_vertex
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist on the processor.
 */

typedef int LB_FIRST_COARSE_OBJ_FN(void *data, 
                                   int num_gid_entries, int num_lid_entries,
                                   LB_ID_PTR global_id, LB_ID_PTR local_id,
                                   int *assigned,
                                   int *num_vert, int *vertices,
                                   int *in_order, int *in_vertex,
                                   int *out_vertex, int *ierr);

typedef int LB_FIRST_COARSE_OBJ_FORT_FN(void *data,
                                    int *num_gid_entries, int *num_lid_entries,
                                    LB_ID_PTR global_id, LB_ID_PTR local_id,
                                    int *assigned,
                                    int *num_vert, int *vertices,
                                    int *in_order, int *in_vertex,
                                    int *out_vertex, int *ierr);

/*****************************************************************************/
/*
 *  Iterator function for coarse objects; return the next coarse object.
 *  This function should be used with LB_FIRST_COARSE_OBJ_FN.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *  Output:
 *    LB_ID_PTR global_id       --  Global ID of the previous coarse object
 *    LB_ID_PTR local_id        --  Local ID of the previous coarse object
 *    LB_ID_PTR next_global_id  --  Global ID of the next coarse object
 *    LB_ID_PTR next_local_id   --  Local ID of the next coarse object
 *    int *assigned             --  indicates processor assignment.
 *                                  1 if the object is currently
 *                                  assigned to this processor; 0 otherwise.
 *                                  For elements that have been refined, it
 *                                  is ignored.
 *    int *num_vert             --  number of vertices in the next object
 *    int *vertices             --  array containing the vertices of the next
                                    coarse object
 *    int *in_vertex            --  the "in" vertex of the next coarse object.
 *                                  It is ignored if in_order==0 in the call
 *                                  to LB_FIRST_COARSE_OBJ_FN.  If the
 *                                  value is -1, a vertex will be selected
 *                                  automatically
 *    int *out_vertex           --  the "out" vertex for the next object;
 *                                  same provisions as in_vertex
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  1 if a valid object is returned; 0 if
 *                                  no more objects exist on the processor.
 */

typedef int LB_NEXT_COARSE_OBJ_FN(void *data, 
                             int num_gid_entries, int num_lid_entries,
                             LB_ID_PTR global_id, LB_ID_PTR local_id,
                             LB_ID_PTR next_global_id, LB_ID_PTR next_local_id,
                             int *assigned,
                             int *num_vert, int *vertices,
                             int *in_vertex,
                             int *out_vertex, int *ierr);

typedef int LB_NEXT_COARSE_OBJ_FORT_FN(void *data, 
                             int *num_gid_entries, int *num_lid_entries,
                             LB_ID_PTR global_id, LB_ID_PTR local_id,
                             LB_ID_PTR next_global_id, LB_ID_PTR next_local_id,
                             int *assigned,
                             int *num_vert, int *vertices,
                             int *in_vertex,
                             int *out_vertex, int *ierr);

/*****************************************************************************/
/*
 *  Function to return the number of children of an element; used for
 *  building a refinement tree.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  Global ID of the object whose number of
 *                                  children is requested
 *    LB_ID_PTR local_id        --  Local ID of the object whose number of
 *                                  children is requested
 *  Output:
 *    int *ierr                 --  error code
 *  Returned value:
 *    int                       --  the number of children
 */

typedef int LB_NUM_CHILD_FN(void *data, 
                            int num_gid_entries, int num_lid_entries,
                            LB_ID_PTR global_id, LB_ID_PTR local_id,
                            int *ierr);

typedef int LB_NUM_CHILD_FORT_FN(void *data, 
                                 int *num_gid_entries, int *num_lid_entries,
                                 LB_ID_PTR global_id, LB_ID_PTR local_id,
                                 int *ierr);

/*****************************************************************************/
/*
 *  Function to return a list of all children of an object.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR parent_gid      --  Global ID of the object whose children
 *                                  are requested
 *    LB_ID_PTR parent_lid      --  Local ID of the object whose children
 *                                  are requested
 *  Output:
 *    LB_ID_PTR child_gids      --  array of Global IDs of the children
 *    LB_ID_PTR child_lids      --  array of Local IDs of the children
 *    int *assigned             --  array indicating processor assignment.
 *                                  1 if the child object is currently
 *                                  assigned to this processor; 0 otherwise.
 *                                  For elements that have been refined, it
 *                                  is ignored.
 *    int *num_vert             --  array containing the number of vertices
 *                                  for each child
 *    int *vertices             --  array containing the vertices for each 
 *                                  child.  If the sum of the number of
 *                                  vertices for children 0 through i-1 is N,
 *                                  then the vertices for child i are in
 *                                  vertices[N:N+num_vert[i]]
 *    LB_REF_TYPE *ref_type     --  indicates what type of refinement was
 *                                  used to create the children
 *    int *in_vertex            --  array containing the "in" vertex for each
 *                                  child, if the user provides them.  It is
 *                                  ignored if ref_type!=LB_IN_ORDER.  For any
 *                                  with the value -1, a vertex will be selected
 *                                  automatically
 *    int *out_vertex           --  array containing the "out" vertex for each
 *                                  child; same provisions as in_vertex
 *    int *ierr                 --  error code
 */

typedef void LB_CHILD_LIST_FN(void *data, 
                              int num_gid_entries, int num_lid_entries,
                              LB_ID_PTR parent_gid, LB_ID_PTR parent_lid,
                              LB_ID_PTR child_gids, LB_ID_PTR child_lids,
                              int *assigned, int *num_vert, int *vertices,
                              LB_REF_TYPE *ref_type, int *in_vertex,
                              int *out_vertex, int *ierr);

typedef void LB_CHILD_LIST_FORT_FN(void *data, 
                                   int *num_gid_entries, int *num_lid_entries,
                                   LB_ID_PTR parent_gid, LB_ID_PTR parent_lid,
                                   LB_ID_PTR child_gids, LB_ID_PTR child_lids,
                                   int *assigned,
                                   int *num_vert, int *vertices,
                                   LB_REF_TYPE *ref_type, int *in_vertex,
                                   int *out_vertex, int *ierr);

/*****************************************************************************/
/*
 *  Function to return the weight of an object.
 *  Input:
 *    void *data                --  pointer to user defined data structure
 *    int num_gid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a global ID
 *    int num_lid_entries       --  number of array entries of type LB_ID_TYPE
 *                                  in a local ID
 *    LB_ID_PTR global_id       --  Global ID of the object whose weight
 *                                  is requested
 *    LB_ID_PTR local_id        --  Local ID of the object whose weight
 *                                  is requested
 *    int wdim                  --  dimension of object weight, or 0 if
 *                                  the weight is not sought.
 *  Output:
 *    float *obj_wgt            --  weight vector for the object
 *                                  (undefined if wdim=0)
 *    int *ierr                 --  error code
 */

typedef void LB_CHILD_WEIGHT_FN(void *data, 
                                int num_gid_entries, int num_lid_entries,
                                LB_ID_PTR global_id, LB_ID_PTR local_id,
                                int wgt_dim, float *obj_wgt, int *ierr);

typedef void LB_CHILD_WEIGHT_FORT_FN(void *data, 
                                     int *num_gid_entries, int *num_lid_entries,
                                     LB_ID_PTR global_id, LB_ID_PTR local_id,
                                     int *wgt_dim, float *obj_wgt, int *ierr);

/*****************************************************************************/
/*****************************************************************************/
/*******  Functions to set-up Zoltan load-balancing data structure  **********/
/*****************************************************************************/
/*****************************************************************************/

/*
 *  Function to initialize values needed in load balancing tools, and
 *  returns which version of the library this is. If the application
 *  uses MPI, call this function after calling MPI_Init. If the
 *  application does not use MPI, this function calls MPI_Init for
 *  use by Zoltan. This function returns the version of
 *  the Zoltan library.
 *  Input:
 *    int argc                   --  Argument count from main()
 *    char **argv                --  Argument list from main()
 *  Output:
 *    float *ver                 --  Version of Zoltan library
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Initialize(int argc, char **argv, float *ver);

/*****************************************************************************/
/*
 *  Function to create a load balancing structure.  May want more than one
 *  structure if using different decompositions with different techniques.
 *  This function allocates and initializes the structure.
 *  Input:
 *    MPI_Comm communicator      --  MPI Communicator to be used for this
 *                                   load-balancing structure.
 *  Returned value:
 *    struct LB_Struct *         --  Pointer to a Zoltan structure.
 *                                   If there is an error, NULL is returned.
 *                                   Any error in this function should be
 *                                   considered fatal.
 */

extern struct LB_Struct *LB_Create(MPI_Comm communicator);

/*****************************************************************************/
/*
 *  Function to free the space associated with a load balancing structure.
 *  The input pointer is set to NULL when the routine returns.
 *  Input:
 *    struct LB_Struct **         --  Pointer to a Zoltan structure.
 */

extern void LB_Destroy(struct LB_Struct **lb);

/*****************************************************************************/
/*
 *  General function to initialize a given Zoltan callback function.
 *  Input:
 *    struct LB_Struct *lb       --  Pointer to a Zoltan structure.
 *    LB_FN_TYPE fn_type         --  Enum type indicating the function to be
 *                                   set.
 *    void (*fn_ptr)()           --  Pointer to the function to be used in the 
 *                                   assignment.
 *    void *data_ptr             --  Pointer to data that Zoltan will
 *                                   pass as an argument to fn(). May be NULL.
 *  Output:
 *    struct LB_Struct *lb       --  Appropriate field set to value in fn_ptr.
 *  Returned value:
 *    int                        --  Error code
 */
extern int LB_Set_Fn(struct LB_Struct *lb, LB_FN_TYPE fn_type,
                     void (*fn_ptr)(), void *data_ptr);

/*
 *  Functions to initialize specific Zoltan callback functions.  One function
 *  exists for each callback function type, as listed in LB_Fn_Type above.
 *  Use of these specific functions enables stricter type checking of the
 *  callback function types.
 *  Input:
 *    struct LB_Struct *lb       --  Pointer to a Zoltan structure.
 *    FN *fn_ptr                 --  Pointer to the function to be used in the 
 *                                   assignment, where FN is one of the
 *                                   callback function typedef'ed above.
 *    void *data_ptr             --  Pointer to data that Zoltan will
 *                                   pass as an argument to fn(). May be NULL.
 *  Output:
 *    struct LB_Struct *lb       --  Appropriate field set to value in fn_ptr.
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Set_Num_Edges_Fn(struct LB_Struct *lb, 
                               LB_NUM_EDGES_FN *fn_ptr, 
                               void *data_ptr);

extern int LB_Set_Edge_List_Fn(struct LB_Struct *lb, 
                               LB_EDGE_LIST_FN *fn_ptr, 
                               void *data_ptr);

extern int LB_Set_Num_Geom_Fn(struct LB_Struct *lb, 
                              LB_NUM_GEOM_FN *fn_ptr, 
                              void *data_ptr);

extern int LB_Set_Geom_Fn(struct LB_Struct *lb, 
                          LB_GEOM_FN *fn_ptr, 
                          void *data_ptr);

extern int LB_Set_Num_Obj_Fn(struct LB_Struct *lb, 
                             LB_NUM_OBJ_FN *fn_ptr, 
                             void *data_ptr);

extern int LB_Set_Obj_List_Fn(struct LB_Struct *lb, 
                              LB_OBJ_LIST_FN *fn_ptr, 
                              void *data_ptr);

extern int LB_Set_First_Obj_Fn(struct LB_Struct *lb, 
                               LB_FIRST_OBJ_FN *fn_ptr, 
                               void *data_ptr);

extern int LB_Set_Next_Obj_Fn(struct LB_Struct *lb, 
                              LB_NEXT_OBJ_FN *fn_ptr, 
                              void *data_ptr);

extern int LB_Set_Num_Border_Obj_Fn(struct LB_Struct *lb, 
                                    LB_NUM_BORDER_OBJ_FN *fn_ptr,
                                    void *data_ptr);

extern int LB_Set_Border_Obj_List_Fn(struct LB_Struct *lb, 
                                     LB_BORDER_OBJ_LIST_FN *fn_ptr, 
                                     void *data_ptr);

extern int LB_Set_First_Border_Obj_Fn(struct LB_Struct *lb, 
                                      LB_FIRST_BORDER_OBJ_FN *fn_ptr, 
                                      void *data_ptr);

extern int LB_Set_Next_Border_Obj_Fn(struct LB_Struct *lb, 
                                     LB_NEXT_BORDER_OBJ_FN *fn_ptr, 
                                     void *data_ptr);

extern int LB_Set_Pre_Migrate_Fn(struct LB_Struct *lb, 
                                 LB_PRE_MIGRATE_FN *fn_ptr, 
                                 void *data_ptr);

extern int LB_Set_Mid_Migrate_Fn(struct LB_Struct *lb, 
                                 LB_MID_MIGRATE_FN *fn_ptr, 
                                 void *data_ptr);

extern int LB_Set_Post_Migrate_Fn(struct LB_Struct *lb, 
                                  LB_POST_MIGRATE_FN *fn_ptr, 
                                  void *data_ptr);

extern int LB_Set_Obj_Size_Fn(struct LB_Struct *lb, 
                              LB_OBJ_SIZE_FN *fn_ptr, 
                              void *data_ptr);

extern int LB_Set_Pack_Obj_Fn(struct LB_Struct *lb, 
                              LB_PACK_OBJ_FN *fn_ptr, 
                              void *data_ptr);

extern int LB_Set_Unpack_Obj_Fn(struct LB_Struct *lb, 
                                LB_UNPACK_OBJ_FN *fn_ptr, 
                                void *data_ptr);

extern int LB_Set_Num_Coarse_Obj_Fn(struct LB_Struct *lb, 
                                    LB_NUM_COARSE_OBJ_FN *fn_ptr, 
                                    void *data_ptr);

extern int LB_Set_Coarse_Obj_List_Fn(struct LB_Struct *lb, 
                                     LB_COARSE_OBJ_LIST_FN *fn_ptr, 
                                     void *data_ptr);

extern int LB_Set_First_Coarse_Obj_Fn(struct LB_Struct *lb, 
                                      LB_FIRST_COARSE_OBJ_FN *fn_ptr, 
                                      void *data_ptr);

extern int LB_Set_Next_Coarse_Obj_Fn(struct LB_Struct *lb, 
                                     LB_NEXT_COARSE_OBJ_FN *fn_ptr, 
                                     void *data_ptr);

extern int LB_Set_Num_Child_Fn(struct LB_Struct *lb, 
                               LB_NUM_CHILD_FN *fn_ptr, 
                               void *data_ptr);

extern int LB_Set_Child_List_Fn(struct LB_Struct *lb, 
                                LB_CHILD_LIST_FN *fn_ptr, 
                                void *data_ptr);

extern int LB_Set_Child_Weight_Fn(struct LB_Struct *lb, 
                                  LB_CHILD_WEIGHT_FN *fn_ptr, 
                                  void *data_ptr);


/*****************************************************************************/
/*
 *  Function to set the load balancing method to be used.
 *  Input:
 *    struct LB_Struct *lb       --  The load balancing structure to which this
 *                                   method applies.
 *    char *string               --  String specifying the desired method.
 *  Output:
 *    struct LB_Struct *lb       --  Appropriate fields set to designated
 *                                   values.
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Set_Method(struct LB_Struct *lb, char *string);

/*****************************************************************************/
/*
 *  Function to change a parameter value within the bowels of Zoltan.
 *  Default values will be used for all parameters not explicitly altered
 *  by a call to this routine.
 *
 *  Input
 *    struct LB_Struct *lb       --  The load balancing structure to which this
 *                                   parameter alteration applies.
 *    char *name                 --  The name of the parameter to have its
 *                                   value changed.
 *    char *val                  --  The new value of the parameter.
 *
 *  Returned value:
 *    int                        --  Error code
 */

extern int LB_Set_Param(struct LB_Struct *lb, char *name, char *val);

/*****************************************************************************/
/*
 *  Function to invoke the load-balancer.
 *
 *  Input:
 *    struct LB_Struct *lb         --  The load balancing structure containing 
 *                                     info for this load-balancing invocation.
 *  Output:
 *    int *changes                 --  This value tells whether the new 
 *                                     decomposition computed by Zoltan differs 
 *                                     from the one given as input to Zoltan.
 *                                     It can be either a one or a zero:
 *                                     zero - No changes to the decomposition
 *                                            were made by the load-balancing
 *                                            algorithm; migration isn't needed.
 *                                     one  - A new decomposition is suggested
 *                                            by the load-balancer; migration
 *                                            is needed to establish the new
 *                                            decomposition.
 *    int *num_gid_entries         --  number of entries of type LB_ID_TYPE
 *                                     in a global ID
 *    int *num_lid_entries         --  number of entries of type LB_ID_TYPE
 *                                     in a local ID
 *    int *num_import              --  The number of non-local objects in the 
 *                                     processor's new decomposition (i.e.,
 *                                     number of objects to be imported).
 *    LB_ID_PTR *import_global_ids --  Pointer to array of Global IDs for the
 *                                     objects to be imported.
 *    LB_ID_PTR *import_local_ids  --  Pointer to array of Local IDs for the 
 *                                     objects to be imported (local to the
 *                                     exporting processor).
 *    int **import_procs           --  Pointer to array of Processor IDs for the
 *                                     objects to be imported (processor IDs of
 *                                     source processor).
 *    int *num_export              --  The number of local objects that must be
 *                                     exported from the processor to establish
 *                                     the new decomposition.
 *    LB_ID_PTR *export_global_ids --  Pointer to array of Global IDs for the
 *                                     objects to be exported from the current
 *                                     processor.
 *    LB_ID_PTR *export_local_ids  --  Pointer to array of Local IDs for the
 *                                     objects to be exported (local to the
 *                                     current processor).
 *    int **export_procs           --  Pointer to array of Processor IDs for the
 *                                     objects to be exported (processor IDs of
 *                                     destination processors).
 *  Returned value:
 *    int                          --  Error code
 */

extern int LB_Balance(struct LB_Struct *lb, int *changes,
                      int *num_gid_entries, int *num_lid_entries,
                      int *num_import, LB_ID_PTR *import_global_ids,
                      LB_ID_PTR *import_local_ids, int **import_procs,
                      int *num_export, LB_ID_PTR *export_global_ids,
                      LB_ID_PTR *export_local_ids, int **export_procs);

/*****************************************************************************/
/*
 *  Routine to compute the inverse map:  Given, for each processor, a list
 *  of non-local objects assigned to the processor, compute the list of objects
 *  that processor needs to export to other processors to establish the new
 *  decomposition.
 *
 *  Input:
 *    struct LB_Struct *lb         --  Load balancing structure for current 
 *                                     balance.
 *    int num_import               --  Number of non-local objects assigned to 
 *                                     the processor in the new decomposition.
 *    LB_ID_PTR import_global_ids  --  Array of global IDs for non-local objects
 *                                     assigned to this processor in the new
 *                                     decomposition.
 *    LB_ID_PTR import_local_ids   --  Array of local IDs for non-local objects
 *                                     assigned to the processor in the new
 *                                     decomposition.
 *    int *import_procs            --  Array of IDs of processors owning the
 *                                     non-local objects that are assigned to
 *                                     this processor in the new decomposition.
 *  Output:
 *    int *num_export              --  The number of local objects that must be
 *                                     exported from the processor to establish
 *                                     the new decomposition.
 *    LB_ID_PTR *export_global_ids --  Pointer to array of Global IDs for the
 *                                     objects to be exported from the current
 *                                     processor.
 *    LB_ID_PTR *export_local_ids  --  Pointer to array of Local IDs for the
 *                                     objects to be exported (local to the
 *                                     current processor).
 *    int **export_procs           --  Pointer to array of Processor IDs for the
 *                                     objects to be exported (processor IDs of
 *                                     destination processors).
 *  Returned value:
 *    int                          --  Error code
 */


extern int LB_Compute_Destinations(struct LB_Struct *lb,
                                   int num_import, 
                                   LB_ID_PTR import_global_ids,
                                   LB_ID_PTR import_local_ids, 
                                   int *import_procs, 
                                   int *num_export, 
                                   LB_ID_PTR *export_global_ids,
                                   LB_ID_PTR *export_local_ids,
                                   int **export_procs);

/*****************************************************************************/
/*
 *  Routine to help perform migration.  If migration pre-processing routine
 *  (LB_PRE_MIGRATE_FN) is specified, this routine first calls that function.
 *  It then calls a function to obtain the size of the migrating objects
 *  (LB_OBJ_SIZE_FN).  The routine next calls an application-specified
 *  object packing routine (LB_PACK_OBJ_FN) for each object
 *  to be exported.  LB_Help_Migrate then develops the needed communication 
 *  map to move the objects to other processors.  It performs the communication 
 *  according to the map. It then calls a mid-migration processing routine
 *  (LB_MID_MIGRATE_FN) if specified, allowing the application to process 
 *  its own data structures before the imported data is unpacked.
 *  It then calls an application-specified object unpacking 
 *  routine (LB_UNPACK_OBJ_FN) for each object imported.  Finally, a 
 *  post-processing function (LB_POST_MIGRATE_FN) is invoked if specified.
 *
 *  Input:
 *    struct LB_Struct *lb       --  Load balancing structure for current 
 *                                   balance.
 *    int num_import             --  Number of non-local objects assigned to the
 *                                   processor in the new decomposition.
 *    LB_ID_PTR import_global_ids--  Array of global IDs for non-local objects
 *                                   assigned to this processor in the new
 *                                   decomposition.
 *    LB_ID_PTR import_local_ids --  Array of local IDs for non-local objects
 *                                   assigned to the processor in the new
 *                                   decomposition.
 *    int *import_procs          --  Array of processor IDs of processors owning
 *                                   the non-local objects that are assigned to
 *                                   this processor in the new decomposition.
 *    int num_export             --  The number of local objects that need to be
 *                                   exported from the processor to establish
 *                                   the new decomposition.
 *    LB_ID_PTR export_global_ids--  Array of Global IDs for the objects to be
 *                                   exported from the current processor.
 *    LB_ID_PTR export_local_ids --  Array of Local IDs for the objects to be 
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
                           int num_import, LB_ID_PTR import_global_ids,
                           LB_ID_PTR import_local_ids, int *import_procs,
                           int num_export, LB_ID_PTR export_global_ids,
                           LB_ID_PTR export_local_ids, int *export_procs);

/*****************************************************************************/
/*
 *  Routine to free the data arrays returned by LB_Balance.  The arrays
 *  are freed and the pointers are set to NULL.
 *
 *  Input:
 *    LB_ID_PTR *import_global_ids --  Pointer to array of global IDs for 
 *                                     imported objects.
 *    LB_ID_PTR *import_local_ids  --  Pointer to array of local IDs for 
 *                                     imported objects.
 *    int **import_procs           --  Pointer to array of processor IDs of 
 *                                     imported objects.
 *    LB_ID_PTR *export_global_ids --  Pointer to array of global IDs for 
 *                                     exported objects.
 *    LB_ID_PTR *export_local_ids  --  Pointer to array of local IDs for 
 *                                     exported objects.
 *    int **export_procs           --  Pointer to array of destination processor
 *                                     IDs of exported objects.
 *  Returned value:
 *    int                          --  Error code
 */
extern int LB_Free_Data(LB_ID_PTR *import_global_ids, 
                        LB_ID_PTR *import_local_ids,
                        int **import_procs,
                        LB_ID_PTR *export_global_ids, 
                        LB_ID_PTR *export_local_ids,
                        int **export_procs);

/*****************************************************************************/
/* 
 * Routine to determine which processor a new point should be assigned to.
 * Note that this only works of the current partition was produced via a
 * geometric algorithm - currently RCB and RIB.
 * 
 * Input:
 *   lb          -- pointer to lb structure
 *   coords      -- vector of coordinates of new point
 *
 * Output:
 *   proc        -- processor that point should be assigned to
 *
 *  Returned value:
 *    int        --  Error code
 */

extern int LB_Point_Assign(struct LB_Struct *lb, double *coords, int *proc);

/*****************************************************************************/
/* 
 * Routine to determine which processors a bounding box intersects.
 * Note that this only works of the current partition was produced via a
 * geometric algorithm - currently RCB and RIB.
 * 
 * Input:
 *   lb                -- pointer to lb structure
 *   xmin, ymin, zmin  -- lower left corner of bounding box
 *   xmax, ymax, zmax  -- upper right corner of bounding box
 *
 * Output:
 *   procs             -- list of processors that box intersects.  
 *                        Note: application is
 *                            responsible for ensuring sufficient memory.
 *   numprocs          -- number of processors box intersects
 *
 *  Returned value:
 *    int        --  Error code
 */

extern int LB_Box_Assign(struct LB_Struct *lb, double xmin, double ymin,
    double zmin, double xmax, double ymax, double zmax, int *procs,
    int *numprocs);

/*****************************************************************************/
/* 
 * Routine to compute statistics about the current balance/partitioning.
 *
 * Input:
 *   lb          - pointer to lb structure
 *   print_stats - if >0, compute and print max and sum of the metrics
 *
 * Output:
 *   nobj      - number of objects (for each proc)
 *   obj_wgt   - obj_wgt[0:vwgt_dim-1] are the object weights (on each proc)
 *   cut_wgt   - cut size/weight (for each proc)
 *   nboundary - number of boundary objects (for each proc)
 *   nadj      - the number of adjacent procs (for each proc)
 *
 * Returned value:
 *   ierr      - error code
 */

extern int LB_Eval (struct LB_Struct *lb, int print_stats, 
     int *nobj, float *obj_wgt, int *ncuts, float *cut_wgt, 
     int *nboundary, int *nadj);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* !__LBI_CONST_H */
