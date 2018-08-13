// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef stk_util_use_cases_UseCaseEnvironement_hpp
#define stk_util_use_cases_UseCaseEnvironement_hpp

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>

#include <stk_util/util/Writer_fwd.hpp>
#include <stk_util/diag/Timer.hpp>

#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <iosfwd>

namespace use_case {

enum LogMask {
  LOG_ALWAYS		= stk::LOG_ALWAYS,
  LOG_TRACE		= stk::LOG_TRACE,
  LOG_TRACE_STATS	= stk::LOG_TRACE_STATS,
  LOG_TRACE_SUB_CALLS	= stk::LOG_TRACE_SUB_CALLS,
  LOG_MEMBERS		= stk::LOG_MEMBERS,

  LOG_SEARCH            = 0x0000010,
  LOG_TRANSFER          = 0x0000020,
  LOG_TIMER             = 0x0000040
};

enum message_type {
  MSG_WARNING = stk::MSG_WARNING,
  MSG_FATAL = stk::MSG_DOOMED,
  MSG_INFORMATION,
  MSG_EXCEPTION
};

enum message_throttle_type {
  MSG_APPLICATION = stk::MSG_APPLICATION,
  MSG_TIME_STEP
};

enum TimerSetMask {
  TIMER_MESH		= 0x00000001,		///< Enable mesh timers
  TIMER_MESH_IO		= 0x00000002,		///< Enable mesh I/O timers
  TIMER_SEARCH		= 0x00000004,		///< Enable search timers
  TIMER_TRANSFER	= 0x00000008,		///< Enable transfer timers
  TIMER_ALL		= 0xFFFFFFFF,		///< Force timer to be active
  TIMER_FORCE		= 0x00000000		///< Force timer to be active
};

stk::diag::TimerSet &timerSet();

stk::diag::Timer &timer();

} // namespace use_case

#endif // stk_util_use_cases_UseCaseEnvironement_hpp
