/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/WallTime.hpp>

#include <sys/time.h>

namespace stk_classic {

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

} // namespace stk_classic
