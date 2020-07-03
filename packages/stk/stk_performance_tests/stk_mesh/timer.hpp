// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef STK_PERFORMANCE_TIMER_HPP
#define STK_PERFORMANCE_TIMER_HPP

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk {
namespace performance_tests {

class Timer
{
public:
  Timer(MPI_Comm comm)
    : communicator(comm),
      iterationStartTime(0.0),
      cumulativeTime(0.0),
      iterationStartHwm(0),
      iterationStartGpuUsage(0),
      meshOperationHwm(1)
  { }

  void start_timing()
  {
    iterationStartTime = stk::wall_time();
    iterationStartHwm = stk::get_max_hwm_across_procs(communicator);
    iterationStartGpuUsage = stk::get_max_gpu_mem_used_across_procs(communicator);
  }

  void update_timing()
  {
    cumulativeTime += stk::wall_dtime(iterationStartTime);
    size_t currentGpuUsage = stk::get_max_gpu_mem_used_across_procs(communicator) - iterationStartGpuUsage;
    size_t currentHwm = stk::get_max_hwm_across_procs(communicator) - iterationStartHwm;
    meshOperationHwm = std::max(meshOperationHwm, currentHwm + currentGpuUsage);
  }

  void print_timing(unsigned iterationCount)
  {
    double timeAll = stk::get_global_sum(communicator, cumulativeTime);
    stk::print_stats_for_performance_compare(std::cout, timeAll, meshOperationHwm, iterationCount, communicator);
  }

private:
  MPI_Comm communicator;
  double iterationStartTime;
  double cumulativeTime;
  size_t iterationStartHwm;
  size_t iterationStartGpuUsage;
  size_t meshOperationHwm;
};

}}

#endif
