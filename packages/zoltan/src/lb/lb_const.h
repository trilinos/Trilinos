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

#include "zoltan.h"

/*
 * Type definitions for load-balancing functions that depend on 
 * load-balancing method.
 */

struct Zoltan_Struct;

typedef int ZOLTAN_LB_FN(struct Zoltan_Struct *, int *, 
                         ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **,
                         int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **);

typedef void ZOLTAN_LB_FREE_DATA_FN(struct Zoltan_Struct *);

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
  BSFC,
  HSFC,
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
#define ZOLTAN_LB_RETURN_LISTS_DEF   ZOLTAN_LB_ALL_LISTS


struct Zoltan_LB_Struct {
  int Return_Lists;            /*  Flag indicating which lists (if any)
                                      should be returned by Zoltan_LB_Balance.*/
  ZOLTAN_LB_METHOD Method;        /*  Method to be used for load balancing.  */
  ZOLTAN_LB_FN *LB_Fn;            /*  Pointer to the function that performs
                                      the load balancing; this ptr is set
                                      based on the method used.              */
  double Imbalance_Tol;           /*  Tolerance to which to load balance;
                                      Imbalance_Tol = 1.1 implies 10% imbalance
                                      is acceptable, i.e. max/avg = 1.1.     */
  void *Data_Structure;           /*  Data structure used by the load
                                      balancer; cast by the method routines
                                      to the appropriate data type.          */
  ZOLTAN_LB_FREE_DATA_FN *Free_Structure;
                                  /*  Pointer to function that frees the
                                      Data_Structure memory.                 */
};

struct Zoltan_Migrate_Struct {
  int Auto_Migrate;                   /*  Flag indicating whether the load
                                          balancer should automatically
                                          help the application
                                          migrate data.  Some applications may
                                          prefer to do it themselves.        */
  /*
   *  Pointers to routines that depend on the application.
   */

  ZOLTAN_PRE_MIGRATE_FN *Pre_Migrate;  /* Function that performs application
                                          specific pre-processing.  Optional
                                          for help-migration.                */
  ZOLTAN_PRE_MIGRATE_FORT_FN *Pre_Migrate_Fort;
                                       /* Fortran version               */
  void *Pre_Migrate_Data;              /* Ptr to user defined data to be
                                          passed to Pre_Migrate()            */
  ZOLTAN_MID_MIGRATE_FN *Mid_Migrate;  /* Function that performs application
                                          specific processing between packing
                                          and unpacking.  Optional
                                          for help-migration.                */
  ZOLTAN_MID_MIGRATE_FORT_FN *Mid_Migrate_Fort;
                                       /* Fortran version               */
  void *Mid_Migrate_Data;              /* Ptr to user defined data to be
                                          passed to Mid_Migrate()            */
  ZOLTAN_POST_MIGRATE_FN *Post_Migrate;/* Function that performs application
                                          specific post-processing.  Optional
                                          for help-migration.                */
  ZOLTAN_POST_MIGRATE_FORT_FN *Post_Migrate_Fort;
                                       /* Fortran version             */
  void *Post_Migrate_Data;             /* Ptr to user defined data to be
                                          passed to Post_Migrate()           */
};

/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */
/* PARTITIONING FUNCTIONS */
extern ZOLTAN_LB_FN Zoltan_RCB;
extern ZOLTAN_LB_FN Zoltan_Octpart;
extern ZOLTAN_LB_FN Zoltan_ParMetis;
extern ZOLTAN_LB_FN Zoltan_Jostle;
extern ZOLTAN_LB_FN Zoltan_Reftree_Part;
extern ZOLTAN_LB_FN Zoltan_RIB;
extern ZOLTAN_LB_FN Zoltan_BSFC;
extern ZOLTAN_LB_FN Zoltan_HSFC;

ZOLTAN_LB_FREE_DATA_FN Zoltan_RCB_Free_Structure;
ZOLTAN_LB_FREE_DATA_FN Zoltan_RIB_Free_Structure;
ZOLTAN_LB_FREE_DATA_FN Zoltan_Oct_Free_Structure;
ZOLTAN_LB_FREE_DATA_FN Zoltan_Reftree_Free_Structure;
ZOLTAN_LB_FREE_DATA_FN Zoltan_HSFC_Free_Structure;

#endif
