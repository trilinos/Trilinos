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

#ifndef __LB_ID_CONST_H
#define __LB_ID_CONST_H

#include "zoltan_id.h"
/*
 *  Macros that access the Zoltan ID manipulation functions.
 *  These macros assume an LB_Struct is available, and access the
 *  appropriate members of LB_Struct in calling the manipulation functions.
 *  Macros are provided for both global IDs (GIDs) and local IDs (LIDs).
 *
 *  The manipulation functions can be called directly; these macros just
 *  make it easier when an LB_Struct is available.
 */

/****************************************************************************/

#define ZOLTAN_LB_INIT_GID(lb,gid) ZOLTAN_INIT_ID((lb)->Num_GID,gid);
#define ZOLTAN_LB_INIT_LID(lb,lid) ZOLTAN_INIT_ID((lb)->Num_LID,lid);

/*
 * Macros for allocating single IDs; IDs are also initialized.
 * If lb->Num_LID is zero, the macro returns NULL.
 */
#define ZOLTAN_LB_MALLOC_GID(lb) ZOLTAN_Malloc_ID((lb)->Num_GID, __FILE__, __LINE__)
#define ZOLTAN_LB_MALLOC_LID(lb) \
    ((lb)->Num_LID \
      ? ZOLTAN_Malloc_ID((lb)->Num_LID, __FILE__, __LINE__) \
      : NULL)

/*
 * Macros for allocating arrays of IDs; arrays are also initialized.
 * If lb->Num_LID is zero, the macro returns NULL.
 */
#define ZOLTAN_LB_MALLOC_GID_ARRAY(lb,num_obj) \
    ZOLTAN_Malloc_ID((num_obj) * (lb)->Num_GID, __FILE__, __LINE__)
#define ZOLTAN_LB_MALLOC_LID_ARRAY(lb,num_obj) \
    ((lb)->Num_LID \
       ? ZOLTAN_Malloc_ID((num_obj) * (lb)->Num_LID, __FILE__, __LINE__) \
       : NULL)

/*
 * Macros for reallocating arrays of IDs.
 */
#define ZOLTAN_LB_REALLOC_GID_ARRAY(lb,ptr,num_obj) \
  (ZOLTAN_ID_PTR) LB_REALLOC(ptr,(num_obj)*(lb)->Num_GID*sizeof(ZOLTAN_ID_TYPE))
#define ZOLTAN_LB_REALLOC_LID_ARRAY(lb,ptr,num_obj) \
  ((lb)->Num_LID \
    ? (ZOLTAN_ID_PTR)LB_REALLOC(ptr, \
                               (num_obj)*(lb)->Num_LID*sizeof(ZOLTAN_ID_TYPE)) \
    : NULL)

/****************************************************************************/
/*
 *  Macros to copy IDs.
 */
#define ZOLTAN_LB_SET_GID(lb,a,b) ZOLTAN_SET_ID((lb)->Num_GID,a,b)
#define ZOLTAN_LB_SET_LID(lb,a,b) ZOLTAN_SET_ID((lb)->Num_LID,a,b)


/****************************************************************************/
/*
 * Macros to print IDs.
 */
#define ZOLTAN_LB_PRINT_GID(lb,a) ZOLTAN_PRINT_ID((lb)->Num_GID,a)
#define ZOLTAN_LB_PRINT_LID(lb,a) ZOLTAN_PRINT_ID((lb)->Num_LID,a)

/****************************************************************************/
/*
 * Macros to compare global IDs. (Comparisons of local IDs are not
 * needed as Zoltan only copies these IDs; it does not use them
 * in its computations.)
 */
#define ZOLTAN_LB_EQ_GID(lb,a,b) ZOLTAN_EQ_ID((lb)->Num_GID,a,b)
#define ZOLTAN_LB_LT_GID(lb,a,b) ZOLTAN_LT_ID((lb)->Num_GID,a,b)
#define ZOLTAN_LB_GT_GID(lb,a,b) ZOLTAN_GT_ID((lb)->Num_GID,a,b)

#endif
