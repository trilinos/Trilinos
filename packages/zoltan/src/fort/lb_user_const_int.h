#ifndef __LB_USER_CONST_H
#define __LB_USER_CONST_H

/*
 *  This file contains the user-defined data types and comparison functions
 *  for global and local IDs used by the application and the
 *  load-balancing library.  Application developers should modify these
 *  data types to match those of identifiers used in their applications.
 *
 *  In this example for Fortran users, LB_GID and LB_LID are both
 *  INTEGER(KIND=LB_INT).
 *
 *  LB_GID are the unique global identifiers for objects in the application.  
 *  The LB_GID are used as identifiers within the load-balancing routines
 *  as well.  Thus, macros defining methods to compare global identifiers 
 *  must be provided.
 *
 *  LB_LID are local identifiers that are not used by the load-balancing
 *  routines.  They are stored with objects in the load-balancing routines,
 *  however, and are passed to the application query functions.  An 
 *  application can provide any values it wants for local identifiers, and
 *  can use them to make access of object information in the query routines
 *  more efficient.
 */

/*
 *  LB_GID and LB_LID data type definitions.
 *  For this example, global IDs (LB_GID) are integers and local IDs (LB_LID)
 *  are also integers.
 */
#define LB_GID int
#define LB_LID int

/*
 * Flags indicating whether or not the IDs are ints.  Set to 1 if the above
 * definition is "int" and 0 otherwise.
 */
#define LB_GID_IS_INT 1
#define LB_LID_IS_INT 1

/*
 *  Macros to copy LB_GIDs and LB_LIDs.
 *  These macros are used by the load-balancing routines to copy LB_GID and
 *  LB_LID values to new storage locations.
 */

#define LB_SET_GID(a,b) (a) = (b)
#define LB_SET_LID(a,b) (a) = (b)

/*
 *  Macros to compare LB_GIDs.
 *  Macros must be provided to test whether two LB_GIDs are equal (EQ),
 *  not equal (NE), less than (LT), less than or equal (LE), 
 *  greater than (GT), and greater than or equal (GE).
 *  Comparison macros are not needed for LB_LIDs as LB_LIDs are not used
 *  within the load-balancing routines.
 */

#define LB_EQ_GID(a,b) ((a) == (b))
#define LB_NE_GID(a,b) ((a) != (b))
#define LB_LT_GID(a,b) ((a) <  (b))
#define LB_LE_GID(a,b) ((a) <= (b))
#define LB_GT_GID(a,b) ((a) >  (b))
#define LB_GE_GID(a,b) ((a) >= (b))

#endif
