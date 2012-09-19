/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/FormatTime.hpp>

#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>

namespace stk {

std::string
formatTime(
  double        time,
  TimeFormat    time_format)
{
  std::stringstream oss;

  if (time < 0.0) {
    time = -time;
    oss << "-";
  }
  
  if ((time_format & TIMEFORMAT_STYLE_MASK) == TIMEFORMAT_SECONDS) {
    if (time_format & TIMEFORMAT_MILLIS)
      oss << std::fixed << std::setprecision(3) << time;
    else
      oss << std::fixed << std::setprecision(0) << time;
  }
  else if ((time_format & TIMEFORMAT_STYLE_MASK) == TIMEFORMAT_HMS) {
    int int_time = int(time);

    if (time >= 3600.0)
      oss << (int_time)/3600 << ':'
             << std::setw(2) << std::setfill('0') << (int_time/60)%60 << ':'
             << std::setw(2) << std::setfill('0') << int_time%60;

    else if (time >= 60.0)
      oss << ((int) (time)/60)%60 << ':'
             << std::setw(2) << std::setfill('0') << int_time%60;


    else
      oss << ((int) time)%60;

    if (time_format & TIMEFORMAT_MILLIS) {
      int milliseconds = int(std::fmod(time, 1.0)*1000.0 + 0.5);

      oss << '.' << std::setw(3) << std::setfill('0') << milliseconds;
    }
  }
  else
    oss << time;

  return oss.str();
}

} // namespace stk
