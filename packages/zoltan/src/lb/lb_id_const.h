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

/*
 *  This file contains the data types and comparison functions
 *  for global and local IDs used by Zoltan.  The basic data type
 *  is LB_ID; both global and local IDs use this type.
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

/*
 * Default value of LB_ID_TYPE; 
 * IDs allocated with LB_MALLOC_GID, LB_MALLOC_LID,
 * LB_MALLOC_GID_ARRAY or LB_MALLOC_LID_ARRAY are initialized with this value.
 */
#define LB_ID_DEFAULT 0

/*
 * Macros for initializing single IDs.
 */

#define LB_INIT_ID(n,id) \
  {int LB_ID_LOOP;       \
  for (LB_ID_LOOP = 0; LB_ID_LOOP < (n); LB_ID_LOOP++)  \
    (id)[LB_ID_LOOP] = LB_ID_DEFAULT;                   \
  }

#define LB_INIT_GID(lb,gid) LB_INIT_ID((lb)->Num_GID,gid);
#define LB_INIT_LID(lb,lid) LB_INIT_ID((lb)->Num_LID,lid);

/*
 * Macros for allocating single IDs; IDs are also initialized.
 * If lb->Num_LID is zero, the macro returns NULL.
 */
#define LB_MALLOC_GID(lb) LB_Malloc_ID((lb)->Num_GID, __FILE__, __LINE__)
#define LB_MALLOC_LID(lb) \
    ((lb)->Num_LID \
      ? LB_Malloc_ID((lb)->Num_LID, __FILE__, __LINE__) \
      : NULL)

/*
 * Macros for allocating arrays of IDs; arrays are also initialized.
 * If lb->Num_LID is zero, the macro returns NULL.
 */
#define LB_MALLOC_GID_ARRAY(lb,num_obj) \
    LB_Malloc_ID((num_obj) * (lb)->Num_GID, __FILE__, __LINE__)
#define LB_MALLOC_LID_ARRAY(lb,num_obj) \
    ((lb)->Num_LID \
       ? LB_Malloc_ID((num_obj) * (lb)->Num_LID, __FILE__, __LINE__) \
       : NULL)

/*
 * Macros for reallocating arrays of IDs.
 */
#define LB_REALLOC_GID_ARRAY(lb,ptr,num_obj) \
  (LB_ID_PTR) LB_REALLOC(ptr, (num_obj) * (lb)->Num_GID * sizeof(LB_ID_TYPE))
#define LB_REALLOC_LID_ARRAY(lb,ptr,num_obj) \
  ((lb)->Num_LID \
     ? (LB_ID_PTR)LB_REALLOC(ptr,(num_obj)*(lb)->Num_LID*sizeof(LB_ID_TYPE)) \
     : NULL)

/****************************************************************************/
/*
 *  Macros to copy IDs.
 */
#define LB_SET_ID(n,a,b)                                            \
   {int LB_ID_LOOP;                                                 \
    for (LB_ID_LOOP = 0; LB_ID_LOOP < (n); LB_ID_LOOP++)            \
      (a)[LB_ID_LOOP] = (b)[LB_ID_LOOP];                            \
   }
#define LB_SET_GID(lb,a,b) LB_SET_ID((lb)->Num_GID,a,b)
#define LB_SET_LID(lb,a,b) LB_SET_ID((lb)->Num_LID,a,b)


/****************************************************************************/
/*
 * Macros to print IDs.
 */
#define LB_PRINT_GID(lb,a) LB_PRINT_ID((lb)->Num_GID,a)
#define LB_PRINT_LID(lb,a) LB_PRINT_ID((lb)->Num_LID,a)

/****************************************************************************/
/*
 * Macros to compare global IDs. (Comparisons of local IDs are not
 * needed as Zoltan only copies these IDs; it does not use them
 * in its computations.)
 */
#define LB_EQ_GID(lb,a,b) LB_EQ_ID((lb)->Num_GID,a,b)
#define LB_LT_GID(lb,a,b) LB_LT_ID((lb)->Num_GID,a,b)
#define LB_GT_GID(lb,a,b) LB_GT_ID((lb)->Num_GID,a,b)


/****************************************************************************/
/*
 *  Prototypes for ID functions in lb_id.c
 */

extern LB_ID_PTR LB_Malloc_ID(int n, char *file, int line);
extern void LB_PRINT_ID(int n, LB_ID_PTR a);
extern int LB_EQ_ID(int n, LB_ID_PTR a, LB_ID_PTR b);
extern int LB_LT_ID(int n, LB_ID_PTR a, LB_ID_PTR b);
extern int LB_GT_ID(int n, LB_ID_PTR a, LB_ID_PTR b);

#endif
