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


#ifndef __LB_CONST_H
#define __LB_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zoltan.h"

/*
 * Type definitions for functions that depend on 
 * load-balancing method.
 */

struct Zoltan_Struct;

typedef int ZOLTAN_LB_FN(struct Zoltan_Struct *, float *, int *, 
                         ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **,
                         int *, ZOLTAN_ID_PTR *, 
                         ZOLTAN_ID_PTR *, int **, int **);

typedef void ZOLTAN_LB_FREE_DATA_FN(struct Zoltan_Struct *);

typedef int ZOLTAN_LB_POINT_ASSIGN_FN(struct Zoltan_Struct *, double *, int *, 
                                      int *);

typedef int ZOLTAN_LB_BOX_ASSIGN_FN(struct Zoltan_Struct *, 
                                    double, double, double,
                                    double, double, double,
                                    int*, int*, int *, int *);

/*
 *  Define the possible load balancing methods allowed.
 */

typedef enum Zoltan_LB_Method {
  NONE = -1,
  RCB,
  OCTPART,
  PARMETIS,
  JOSTLE,
  REFTREE,
  RIB,
  HSFC,
  HG,
  PHG,
  ZOLTAN_LB_MAX_METHODS          /*  This entry should always be last.      */
} ZOLTAN_LB_METHOD;


/*
 * Values indicating which lists (import, export, both, or none) should
 * be returned by Zoltan_LB_Balance.  ZOLTAN_LB_NO_LISTS must always be zero; 
 * other values should always be greater than zero.
 */
#define ZOLTAN_LB_NO_LISTS 0
#define ZOLTAN_LB_IMPORT_LISTS 1
#define ZOLTAN_LB_EXPORT_LISTS 2
#define ZOLTAN_LB_ALL_LISTS 3

/*
 ******************************************************
 * Define default values for key parameters.
 ******************************************************
 */
#define ZOLTAN_LB_IMBALANCE_TOL_DEF  1.1
#define ZOLTAN_AUTO_MIGRATE_DEF   FALSE
#define ZOLTAN_MIGRATE_ONLY_PROC_CHANGES_DEF  1
#define ZOLTAN_LB_RETURN_LISTS_DEF   ZOLTAN_LB_ALL_LISTS

/* Struct for partition size info. */
struct Zoltan_part_info {
  float Size;          /*  Desired partition size. */
  int Part_id;         /*  Partition id.           */
  int Idx;             /*  Partition weight index. */
  int Global_num;      /*  Global partition numbers? */
};

/* LB_struct. Contains all information about a load balancing "object" */

struct Zoltan_LB_Struct {
  int Part_Info_Len;              /*  Actual length of Part_Info arrays. */
  int Part_Info_Max_Len;          /*  Allocated length of Part_Info arrays. */
  struct Zoltan_part_info *Part_Info; /*  Array of partition size info.  */
  int Num_Global_Parts;           /*  The total number of partitions.
                                      Set in Zoltan_LB_Build_PartDist.       */
  int Num_Global_Parts_Param;     /*  The number of global partitions specified.
                                      If parameter NUM_LOCAL_PARTITIONS or 
                                      NUM_GLOBAL_PARTITIONS is not set,
                                      Num_Global_Parts_Param == Num_Proc.    */
  int Num_Local_Parts_Param;      /*  The number of local partitions specified.
                                      If parameter NUM_LOCAL_PARTITIONS or 
                                      NUM_GLOBAL_PARTITIONS is not set,
                                      Num_Local_Parts_Param == -1.           */
  int Prev_Global_Parts_Param;    /*  The previous values of
                                      Num_Global_Parts_Param.  Stored to 
                                      prevent unnecessary re-creations of 
                                      PartDist. */
  int Prev_Local_Parts_Param;     /*  The previous values of
                                      Num_Local_Parts_Param.  Stored to 
                                      prevent unnecessary re-creations of 
                                      PartDist. */
  int Single_Proc_Per_Part;       /*  Flag indicating whether a partition can
                                      be spread across multiple processors.
                                      Happens only when NUM_GLOBAL_PARTITIONS
                                      is set to be < zz->Num_Proc.           */
  int Remap_Flag;                 /*  Flag indicating whether partitions
                                      should be remapped to reduce data mvmt. */
  int *Remap;                     /*  Remapping array; relabels computed 
                                      partitions to decrease data mvmt. */
  int Return_Lists;               /*  Flag indicating which lists (if any)
                                      should be returned by Zoltan_LB_Balance.*/
  int Uniform_Parts;              /*  Flag indicating whether partitions are
                                      uniformly sized. */
  int *PartDist;                  /*  Array describing distribution of 
                                      partitions to processors.  
                                      If Single_Proc_Per_Part, partition i
                                      is located on processor PartDist[i].
                                      If !Single_Proc_Per_Part, partition i
                                      is located on processors PartDist[i] to
                                      PartDist[i+1]-1. */
  int *ProcDist;                  /*  Array describing distribution of 
                                      processors to partitions.  
                                      If processor i has zero partitions,
                                      ProcDist[i] = -1.  Otherwise,
                                      ProcDist[i] has the lowest partition
                                      number of partitions on processor i.  */
  ZOLTAN_LB_METHOD Method;        /*  Method to be used for load balancing.  */ 
  ZOLTAN_LB_FN *LB_Fn;            /*  Pointer to the function that performs
                                      the load balancing; this ptr is set
                                      based on the method used.              */
  float *Imbalance_Tol;           /*  Tolerance to which to load balance;
                                      Imbalance_Tol = 1.1 implies 10% imbalance
                                      is acceptable, i.e. max/avg = 1.1.     
                                      Imbalance_Tol may be an array of
                                      dimension Obj_Weight_Dim.              */
  int  Imb_Tol_Len;               /*  Length of Imbalance_Tol array.         */
  void *Data_Structure;           /*  Data structure used by the load
                                      balancer; cast by the method routines
                                      to the appropriate data type.          */
  ZOLTAN_LB_FREE_DATA_FN *Free_Structure;
                                  /*  Pointer to function that frees the
                                      Data_Structure memory.                 */
  ZOLTAN_LB_POINT_ASSIGN_FN *Point_Assign;  
                                  /*  Pointer to the function that performs
                                      Point_Assign; this ptr is set based on 
                                      the method used.                       */
  ZOLTAN_LB_BOX_ASSIGN_FN *Box_Assign;      
                                  /*  Pointer to the function that performs
                                      Box_Assign; this ptr is set based on 
                                      the method used.                       */
};

struct Zoltan_Migrate_Struct {
  int Auto_Migrate;                   /*  Flag indicating whether the load
                                          balancer should automatically
                                          help the application
                                          migrate data.  Some applications may
                                          prefer to do it themselves.        */
  int Only_Proc_Changes;              /*  Pack and unpack objects during
                                          migration ONLY if they are assigned
                                          to a new processor.  If partition
                                          number changes but processor does
                                          not, do not pack and unpack.       */
  /*
   *  Pointers to routines that depend on the application.
   */

  ZOLTAN_PRE_MIGRATE_PP_FN *Pre_Migrate_PP;
                                       /* Function that performs application
                                          specific pre-processing (including
                                          partition lists).  Optional
                                          for migration.                */
  ZOLTAN_PRE_MIGRATE_PP_FORT_FN *Pre_Migrate_PP_Fort;
                                       /* Fortran version               */
  void *Pre_Migrate_PP_Data;         /* Ptr to user defined data to be
                                          passed to Pre_Migrate_PP()       */
  ZOLTAN_MID_MIGRATE_PP_FN *Mid_Migrate_PP;
                                       /* Function that performs application
                                          specific processing  (including
                                          partition lists) between packing
                                          and unpacking.  Optional
                                          for migration.                */
  ZOLTAN_MID_MIGRATE_PP_FORT_FN *Mid_Migrate_PP_Fort;
                                       /* Fortran version               */
  void *Mid_Migrate_PP_Data;         /* Ptr to user defined data to be
                                          passed to Mid_Migrate_PP()       */
  ZOLTAN_POST_MIGRATE_PP_FN *Post_Migrate_PP;
                                       /* Function that performs application
                                          specific post-processing (including 
                                          partition lists).  Optional
                                          for migration.                */
  ZOLTAN_POST_MIGRATE_PP_FORT_FN *Post_Migrate_PP_Fort;
                                       /* Fortran version             */
  void *Post_Migrate_PP_Data;        /* Ptr to user defined data to be
                                          passed to Post_Migrate_PP()      */
  ZOLTAN_PRE_MIGRATE_FN *Pre_Migrate;  /* Function that performs application
                                          specific pre-processing.  Optional
                                          for migration.                */
  ZOLTAN_PRE_MIGRATE_FORT_FN *Pre_Migrate_Fort;
                                       /* Fortran version               */
  void *Pre_Migrate_Data;              /* Ptr to user defined data to be
                                          passed to Pre_Migrate()            */
  ZOLTAN_MID_MIGRATE_FN *Mid_Migrate;  /* Function that performs application
                                          specific processing between packing
                                          and unpacking.  Optional
                                          for migration.                */
  ZOLTAN_MID_MIGRATE_FORT_FN *Mid_Migrate_Fort;
                                       /* Fortran version               */
  void *Mid_Migrate_Data;              /* Ptr to user defined data to be
                                          passed to Mid_Migrate()            */
  ZOLTAN_POST_MIGRATE_FN *Post_Migrate;/* Function that performs application
                                          specific post-processing.  Optional
                                          for migration.                */
  ZOLTAN_POST_MIGRATE_FORT_FN *Post_Migrate_Fort;
                                       /* Fortran version             */
  void *Post_Migrate_Data;             /* Ptr to user defined data to be
                                          passed to Post_Migrate()           */
};

/*****************************************************************************/
/* PROTOTYPES */

extern int Zoltan_LB_Set_LB_Method(struct Zoltan_Struct *, char *);
extern void Zoltan_LB_Free_Struct(struct Zoltan_LB_Struct *);
extern int Zoltan_LB_Part_To_Proc(struct Zoltan_Struct *, int, ZOLTAN_ID_PTR);
extern int Zoltan_LB_Proc_To_Part(struct Zoltan_Struct *, int, int *, int *);
extern int Zoltan_LB_Get_Part_Sizes(struct Zoltan_Struct *, int, int, float *);
extern int Zoltan_LB_Build_PartDist(struct Zoltan_Struct *);
extern int Zoltan_LB_Remap(struct Zoltan_Struct *, int *, int, int *, int *,
  int *, int);

/* PARTITIONING FUNCTIONS */
extern ZOLTAN_LB_FN Zoltan_RCB;
extern ZOLTAN_LB_FN Zoltan_Octpart;
extern ZOLTAN_LB_FN Zoltan_ParMetis;
extern ZOLTAN_LB_FN Zoltan_Jostle;
extern ZOLTAN_LB_FN Zoltan_Reftree_Part;
extern ZOLTAN_LB_FN Zoltan_RIB;
extern ZOLTAN_LB_FN Zoltan_HSFC;
extern ZOLTAN_LB_FN Zoltan_HG;
extern ZOLTAN_LB_FN Zoltan_PHG;

/* FREE DATA_STRUCTURE FUNCTIONS */
extern ZOLTAN_LB_FREE_DATA_FN Zoltan_RCB_Free_Structure;
extern ZOLTAN_LB_FREE_DATA_FN Zoltan_RIB_Free_Structure;
extern ZOLTAN_LB_FREE_DATA_FN Zoltan_Oct_Free_Structure;
extern ZOLTAN_LB_FREE_DATA_FN Zoltan_Reftree_Free_Structure;
extern ZOLTAN_LB_FREE_DATA_FN Zoltan_HSFC_Free_Structure;
extern ZOLTAN_LB_FREE_DATA_FN Zoltan_HG_Free_Structure;
extern ZOLTAN_LB_FREE_DATA_FN Zoltan_PHG_Free_Structure;

/* POINT_ASSIGN FUNCTIONS */
extern ZOLTAN_LB_POINT_ASSIGN_FN Zoltan_RB_Point_Assign;
extern ZOLTAN_LB_POINT_ASSIGN_FN Zoltan_HSFC_Point_Assign;

/* BOX_ASSIGN FUNCTIONS */
extern ZOLTAN_LB_BOX_ASSIGN_FN Zoltan_RB_Box_Assign;
extern ZOLTAN_LB_BOX_ASSIGN_FN Zoltan_HSFC_Box_Assign;

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
