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

#include <sstream>

#include <stk_util/diag/TimerMetricTraits.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/MallocUsed.h>
#include <stk_util/environment/FormatTime.hpp>
#include <stk_util/environment/FormatMemorySize.hpp>

namespace stk {
namespace diag {

namespace {

int s_timeFormat = TIMEFORMAT_HMS | TIMEFORMAT_MILLIS;

} // namespace <empty>


int
getTimerTimeFormat() 
{
  return s_timeFormat;
}

void
setTimerTimeFormat(
  int           time_format)
{
  s_timeFormat = time_format;
}


MetricTraits<LapCount>::Type
MetricTraits<LapCount>::value_now()
{
  return 1;
}

MetricTraits<CPUTime>::Type
MetricTraits<CPUTime>::value_now()
{
  return stk::cpu_time();
}

MetricTraits<WallTime>::Type
MetricTraits<WallTime>::value_now()
{
  return stk::wall_time();
}

MetricTraits<MPICount>::Type
MetricTraits<MPICount>::value_now()
{
  return 0;
}

MetricTraits<MPIByteCount>::Type
MetricTraits<MPIByteCount>::value_now()
{
  return 0;
}

MetricTraits<HeapAlloc>::Type
MetricTraits<HeapAlloc>::value_now()
{
  return ::malloc_used();
}

std::string
MetricTraits<LapCount>::table_header() {
  return "Count";
}

std::string
MetricTraits<CPUTime>::table_header() {
  return "CPU Time";
}

std::string
MetricTraits<WallTime>::table_header() {
  return "Wall Time";
}

std::string
MetricTraits<MPICount>::table_header() {
  return "MPI Count";
}

std::string
MetricTraits<MPIByteCount>::table_header() {
  return "MPI Byte Count";
}

std::string
MetricTraits<HeapAlloc>::table_header() {
  return "Heap Allocated";
}


std::string
MetricTraits<CPUTime>::format(
  MetricTraits<CPUTime>::Type           time)
{
  return formatTime(time, getTimerTimeFormat());
}


std::string
MetricTraits<WallTime>::format(
  MetricTraits<WallTime>::Type          time)
{
  return formatTime(time, getTimerTimeFormat());
}


std::string
MetricTraits<MPICount>::format(
  MetricTraits<MPICount>::Type          count)
{
  std::stringstream strout;

  strout << count;

  return strout.str();
}


std::string
MetricTraits<MPIByteCount>::format(
  MetricTraits<MPIByteCount>::Type      count)
{
  std::stringstream strout;

  strout << count;

  return strout.str();
}

std::string
MetricTraits<HeapAlloc>::format(
  MetricTraits<HeapAlloc>::Type         count)
{
  return format_memory_size(count);
}

} // namespace diag
} // namespace stk
