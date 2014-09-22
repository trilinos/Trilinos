// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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
    int int_time = static_cast<int>(time);

    if (time >= 3600.0)
      oss << (int_time)/3600 << ':'
             << std::setw(2) << std::setfill('0') << (int_time/60)%60 << ':'
             << std::setw(2) << std::setfill('0') << int_time%60;

    else if (time >= 60.0)
      oss << (int_time/60)%60 << ':'
             << std::setw(2) << std::setfill('0') << int_time%60;


    else
      oss << int_time%60;

    if (time_format & TIMEFORMAT_MILLIS) {
      int milliseconds = static_cast<int>(std::fmod(time, 1.0)*1000.0 + 0.5);

      oss << '.' << std::setw(3) << std::setfill('0') << milliseconds;
    }
  }
  else
    oss << time;

  return oss.str();
}

} // namespace stk
