#ifndef __LB_USER_CONST_H
#define __LB_USER_CONST_H

/*
 *  This file contains the user-defined data types and comparison functions
 *  for global and local IDs used by the application and the
 *  load-balancing library.  Application developers should modify these
 *  data types to match those of identifiers used in their applications.
 *
 *  In this example for Fortran users, LB_LID is INTEGER(KIND=LB_INT), and
 *  LB_GID is a structure containing a pair of integers, with the second
 *  integer acting as a high order word.  The LIDs are defined in this .h
 *  file, and the GIDs are defined in the corresponding .f90 file.
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
 *  For this example, local IDs (LB_LID) are integers.  global IDs (LB_GID)
 *  is a Fortran 90 structure defined in lb_user_const.f90, and here is a
 *  C type that consumes the same number of bytes.  In general, determine
 *  the number of bytes consumed by the structure and replace "8" with that.
 */
struct LB_fortran_structure_holder {char a[8];};
typedef struct LB_fortran_structure_holder LB_GID;
#define LB_LID int

/*
 * Flags indicating whether or not the IDs are ints.  Set to 1 if the above
 * definition is "int" and 0 otherwise.
 */
#define LB_GID_IS_INT 0
#define LB_LID_IS_INT 1

/*  The following macros for LB_GIDs call the Fortran routines.  The form
 *  of the name mangling will depend on the compiler.
 */

/*
 *  Macros to copy LB_GIDs and LB_LIDs.
 *  These macros are used by the load-balancing routines to copy LB_GID and
 *  LB_LID values to new storage locations.
 */

#define LB_SET_GID(a,b) lb_user_const_MP_lb_set_gid(&(a),&(b))
#define LB_SET_LID(a,b) (a) = (b)

/*
 *  Macros to compare LB_GIDs.
 *  Macros must be provided to test whether two LB_GIDs are equal (EQ),
 *  not equal (NE), less than (LT), less than or equal (LE), 
 *  greater than (GT), and greater than or equal (GE).
 *  Comparison macros are not needed for LB_LIDs as LB_LIDs are not used
 *  within the load-balancing routines.
 */

#define LB_EQ_GID(a,b) lb_user_const_MP_lb_eq_gid(&(a),&(b))
#define LB_NE_GID(a,b) lb_user_const_MP_lb_ne_gid(&(a),&(b))
#define LB_LT_GID(a,b) lb_user_const_MP_lb_lt_gid(&(a),&(b))
#define LB_LE_GID(a,b) lb_user_const_MP_lb_le_gid(&(a),&(b))
#define LB_GT_GID(a,b) lb_user_const_MP_lb_gt_gid(&(a),&(b))
#define LB_GE_GID(a,b) lb_user_const_MP_lb_ge_gid(&(a),&(b))

#endif
