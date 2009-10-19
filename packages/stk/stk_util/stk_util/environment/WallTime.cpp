#include <stk_util/environment/WallTime.hpp>

#include <sys/time.h>

namespace stk {

double
wall_time()
{
  timeval tp;
  struct timezone tz;
  ::gettimeofday(&tp, &tz);

  double seconds = tp.tv_sec;
  double milliseconds = tp.tv_usec*1.0e-6;
 
  return seconds + milliseconds;
}


double
wall_dtime(double &t)
{
  const double tnew = wall_time();

  const double dt = tnew - t;

  t = tnew ;

  return dt ;
}

} // namespace stk
