#include <stk_util/environment/CPUTime.hpp>

#include <sys/resource.h>

namespace stk {

double
cpu_time()
{
#if defined(REDS)
  struct rusage my_rusage;

  ::getrusage(RUSAGE_SELF, &my_rusage);

  double seconds = my_rusage.ru_utime.tv_sec;
  double micro_seconds = my_rusage.ru_utime.tv_usec;
  
  return seconds + micro_seconds*1.0e-6;

#elif ! defined(__PGI)
  struct rusage my_rusage;

  ::getrusage(RUSAGE_SELF, &my_rusage);

  double seconds = my_rusage.ru_utime.tv_sec + my_rusage.ru_stime.tv_sec;
  double micro_seconds = my_rusage.ru_utime.tv_usec + my_rusage.ru_stime.tv_usec;
  
  return seconds + micro_seconds*1.0e-6;
#else
  return 0.0;
#endif
}

} // namespace stk
