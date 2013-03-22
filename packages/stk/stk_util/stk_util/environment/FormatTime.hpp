/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_environment_FormatTime_hpp
#define stk_util_environment_FormatTime_hpp

#include <string>

namespace stk {

typedef unsigned long TimeFormat;

enum {
  TIMEFORMAT_NONE       = 0x00,
  TIMEFORMAT_HMS        = 0x01,
  TIMEFORMAT_SECONDS    = 0x02,
  TIMEFORMAT_STYLE_MASK = 0x0F,
  
  TIMEFORMAT_MILLIS     = 0x10
};

std::string formatTime(double time, TimeFormat time_format = TIMEFORMAT_HMS | TIMEFORMAT_MILLIS);

} // namespace stk

#endif // stk_util_environment_FormatTime_hpp
