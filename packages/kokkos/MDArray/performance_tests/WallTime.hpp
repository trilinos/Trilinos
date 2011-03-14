#ifndef performance_test_WallTime_hpp
#define performance_test_WallTime_hpp

#include <sys/time.h>

namespace timer {

inline double walltime() {
  timeval tp;
  struct timezone tz;
  ::gettimeofday(&tp, &tz);

  return ( tp.tv_sec               //seconds
         + tp.tv_usec * 1.0e-6 );  //milliseconds

}

}

#endif // performance_test_WallTime_hpp
