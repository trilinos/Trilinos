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
 *  This file contains the data types and comparison functions
 *  for global and local IDs used by Zoltan.  The basic data type
 *  is ZOLTAN_ID; both global and local IDs use this type.
 *
 *  Global IDs (GIDs) are unique identifiers for objects in the application.  
 *  The Global IDs are used as identifiers within Zoltan as well.
 *
 *  Local IDs (LIDs) are identifiers that are not used by Zoltan.
 *  However, they are stored with objects in Zoltan, and are passed to the 
 *  application query functions.  An application can provide any values it 
 *  wants for local identifiers, and can use them to make access of object 
 *  information in the query routines more efficient.
 */

/****************************************************************************/

#define LB_INIT_GID(lb,gid) ZOLTAN_INIT_ID((lb)->Num_GID,gid);
#define LB_INIT_LID(lb,lid) ZOLTAN_INIT_ID((lb)->Num_LID,lid);

/*
 * Macros for allocating single IDs; IDs are also initialized.
 * If lb->Num_LID is zero, the macro returns NULL.
 */
#define LB_MALLOC_GID(lb) ZOLTAN_Malloc_ID((lb)->Num_GID, __FILE__, __LINE__)
#define LB_MALLOC_LID(lb) \
    ((lb)->Num_LID \
      ? ZOLTAN_Malloc_ID((lb)->Num_LID, __FILE__, __LINE__) \
      : NULL)

/*
 * Macros for allocating arrays of IDs; arrays are also initialized.
 * If lb->Num_LID is zero, the macro returns NULL.
 */
#define LB_MALLOC_GID_ARRAY(lb,num_obj) \
    ZOLTAN_Malloc_ID((num_obj) * (lb)->Num_GID, __FILE__, __LINE__)
#define LB_MALLOC_LID_ARRAY(lb,num_obj) \
    ((lb)->Num_LID \
       ? ZOLTAN_Malloc_ID((num_obj) * (lb)->Num_LID, __FILE__, __LINE__) \
       : NULL)

/*
 * Macros for reallocating arrays of IDs.
 */
#define LB_REALLOC_GID_ARRAY(lb,ptr,num_obj) \
  (ZOLTAN_ID_PTR) LB_REALLOC(ptr,(num_obj)*(lb)->Num_GID*sizeof(ZOLTAN_ID_TYPE))
#define LB_REALLOC_LID_ARRAY(lb,ptr,num_obj) \
  ((lb)->Num_LID \
    ? (ZOLTAN_ID_PTR)LB_REALLOC(ptr, \
                               (num_obj)*(lb)->Num_LID*sizeof(ZOLTAN_ID_TYPE)) \
    : NULL)

/****************************************************************************/
/*
 *  Macros to copy IDs.
 */
#define LB_SET_GID(lb,a,b) ZOLTAN_SET_ID((lb)->Num_GID,a,b)
#define LB_SET_LID(lb,a,b) ZOLTAN_SET_ID((lb)->Num_LID,a,b)


/****************************************************************************/
/*
 * Macros to print IDs.
 */
#define LB_PRINT_GID(lb,a) ZOLTAN_PRINT_ID((lb)->Num_GID,a)
#define LB_PRINT_LID(lb,a) ZOLTAN_PRINT_ID((lb)->Num_LID,a)

/****************************************************************************/
/*
 * Macros to compare global IDs. (Comparisons of local IDs are not
 * needed as Zoltan only copies these IDs; it does not use them
 * in its computations.)
 */
#define LB_EQ_GID(lb,a,b) ZOLTAN_EQ_ID((lb)->Num_GID,a,b)
#define LB_LT_GID(lb,a,b) ZOLTAN_LT_ID((lb)->Num_GID,a,b)
#define LB_GT_GID(lb,a,b) ZOLTAN_GT_ID((lb)->Num_GID,a,b)

#endif
