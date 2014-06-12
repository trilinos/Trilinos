/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/CPUTime.hpp>

#include <sys/resource.h>

namespace stk_classic {

double
cpu_time()
{
#if defined(REDS)
  struct rusage my_rusage;

  ::getrusage(RUSAGE_SELF, &my_rusage);

  double seconds = my_rusage.ru_utime.tv_sec;
  double micro_seconds = my_rusage.ru_utime.tv_usec;
  
  return seconds + micro_seconds*1.0e-6;

#else
  struct rusage my_rusage;

  ::getrusage(RUSAGE_SELF, &my_rusage);

  double seconds = my_rusage.ru_utime.tv_sec + my_rusage.ru_stime.tv_sec;
  double micro_seconds = my_rusage.ru_utime.tv_usec + my_rusage.ru_stime.tv_usec;
  
  return seconds + micro_seconds*1.0e-6;

#endif
}

} // namespace stk_classic
