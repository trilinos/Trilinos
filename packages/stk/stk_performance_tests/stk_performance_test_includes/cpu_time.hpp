#ifndef samba_performance_cpu_time_hpp
#define samba_performance_cpu_time_hpp

#include <sys/resource.h>

namespace {

inline double cpu_time()
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

} //namespace

#endif //samba_performance_cpu_time_hpp

