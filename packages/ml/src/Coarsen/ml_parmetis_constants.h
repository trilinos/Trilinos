/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef ML_PARMETIS_CONSTANTS_H
#define ML_PARMETIS_CONSTANTS_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*
   This header defines types used in ML's (Par)METIS interfaces and must be included whenever (Par)METIS data types are used.  This is necessary
   because (Par)METIS does not maintain backwards compatibility between major versions.

   The typedef is based first on ParMETIS, then falls back to METIS.
*/

#if defined(PARMETIS_MAJOR_VERSION)

#  if PARMETIS_MAJOR_VERSION == 3
#     define indextype idxtype
#  elif PARMETIS_MAJOR_VERSION == 4
#     define indextype idx_t
#  endif

#elif defined(METIS_VER_MAJOR)

#  if METIS_VER_MAJOR == 5
#     define indextype idx_t
#  endif

#else

#  define indextype idxtype

#endif

/* METIS 4.x doesn't supply a version macro */
#if !defined(METIS_VER_MAJOR)
# define METIS_VER_MAJOR 4
#endif

#endif /*ifndef ML_PARMETIS_CONSTANTS_H*/
