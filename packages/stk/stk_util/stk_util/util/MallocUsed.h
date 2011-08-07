/*------------------------------------------------------------------------*/
/*                 Copyright 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_UTIL_MALLOCUSED_H
#define STK_UTIL_UTIL_MALLOCUSED_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
size_t malloc_used();
size_t malloc_footprint();
size_t malloc_max_footprint();
#else
inline size_t malloc_used() {
  return 0;
}
inline size_t malloc_footprint() {
  return 0;
}
inline size_t malloc_max_footprint() {
  return 0;
}
#endif

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* STK_UTIL_UTIL_MALLOCUSED_H */
