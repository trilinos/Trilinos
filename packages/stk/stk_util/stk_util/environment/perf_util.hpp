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

#ifndef STK_PERF_UTIL_H
#define STK_PERF_UTIL_H

#include <string>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>

namespace stk {

// Functions that take iterators on raw data (names and times) and
// produce results that will get picked-up by the parsing tools that
// produce the performance charts on the STK performance website.

template <typename NameItr, typename TimerItr>
inline
void print_timers(NameItr name_itr, TimerItr timer_itr, int count)
{
  std::cout << std::endl;

  for (int i = 0; i < count; ++i, ++name_itr, ++timer_itr) {
    std::cout << "STKPERF: " << *name_itr << ": " << *timer_itr << std::endl;
  }
}

template <typename NameItr, typename MemoryItr>
inline
void print_memory(NameItr name_itr, MemoryItr memory_itr, int count)
{
  std::cout << std::endl;

  for (int i = 0; i < count; ++i, ++name_itr, ++memory_itr) {
    std::cout << "STKPERF: " << *name_itr << ": " << *memory_itr << " (" << human_bytes(*memory_itr) << ")" << std::endl;
  }
}

// Use this version if providing your own memory numbers
template <typename NameItr, typename TimerItr, typename MemItr, typename MemNameItr>
inline
void print_timers_and_memory(NameItr name_itr, TimerItr timer_itr, int timer_count, MemNameItr mem_name_itr, MemItr mem_itr, int mem_count)
{
  print_timers(name_itr, timer_itr, timer_count);

  print_memory(mem_name_itr, mem_itr, mem_count);
}

// Prefer this version unless you are sure you are not interested in memory.
// This will automatically provide memory numbers for current memory usage
// and high-water mark.
template <typename NameItr, typename TimerItr>
inline
void print_timers_and_memory(NameItr name_itr, TimerItr timer_itr, int count)
{
  static const char* mem_names[2] = {"Current memory", "Memory high water"};
  size_t mem_array[2];
  stk::get_memory_usage(mem_array[0], mem_array[1]);

  print_timers_and_memory(name_itr, timer_itr, count, &mem_names[0], &mem_array[0], 2);
}

// Use this version to produce output compatible with the compare_test_result_timings.py
// It will also do parallel reductions for you.
inline
void parallel_print_time_without_output_and_hwm(MPI_Comm comm, double time_on_this_proc, std::ostream& out = std::cout)
{
  size_t hwm_max = 0, hwm_min = 0, hwm_avg = 0;
  get_memory_high_water_mark_across_processors(comm, hwm_max, hwm_min, hwm_avg);

  size_t avail_max = 0, avail_min = 0, avail_avg = 0;
  get_memory_available_across_processors(comm, avail_max, avail_min, avail_avg);

  double max_time = 0.0, min_time = 0.0, avg_time = 0.0;
  get_max_min_avg(comm, time_on_this_proc, max_time, min_time, avg_time);

  int rank = stk::parallel_machine_rank(comm);

  const double bytes_in_MB = 1024*1024;

  if (rank == 0) {
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Min High-water memory usage " << hwm_min / bytes_in_MB << " MB" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Avg High-water memory usage " << hwm_avg / bytes_in_MB << " MB" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Max High-water memory usage " << hwm_max / bytes_in_MB << " MB\n" << std::endl;

    out << std::setw(6) << std::fixed << std::setprecision(1) << "Min Available memory per processor " << avail_min / bytes_in_MB << " MB" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Avg Available memory per processor " << avail_avg / bytes_in_MB << " MB" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Max Available memory per processor " << avail_max / bytes_in_MB << " MB\n" << std::endl;

    out << std::setw(6) << std::fixed << std::setprecision(4) << "Min No-output time " << min_time << " sec" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(4) << "Avg No-output time " << avg_time << " sec" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(4) << "Max No-output time " << max_time << " sec" << std::endl;
  }
}

inline
size_t get_max_hwm_across_procs(MPI_Comm comm)
{
    size_t now = 0, localHwm = 0, maxHwm = 0;
    stk::get_memory_usage(now, localHwm);
    stk::all_reduce_max(comm, &localHwm, &maxHwm, 1);
    return maxHwm;
}

inline
size_t get_max_gpu_mem_used_across_procs(MPI_Comm comm)
{
  size_t used = 0, free = 0, globalMaxUsed = 0;
  stk::get_gpu_memory_info(used, free);
  stk::all_reduce_max(comm, &used, &globalMaxUsed, 1);
  return globalMaxUsed;
}

inline
double get_max_time_across_procs(double time_on_this_proc, MPI_Comm comm)
{
    double maxTime = 0.0;
    stk::all_reduce_max(comm, &time_on_this_proc, &maxTime, 1);
    return maxTime;
}

template <typename VALUETYPE>
void print_tagged_stat(std::ostream& out, const std::string &tag, VALUETYPE value)
{
    out << std::setw(11) << std::fixed << std::setprecision(6) << tag << value << std::endl;
}

inline
void print_stats(std::ostream& out, double maxTime, size_t maxHwm, unsigned numSteps)
{
    print_tagged_stat(out, "### Total Wall Clock Run Time Used ###: ", maxTime);
    print_tagged_stat(out, "### Total Number of Steps Taken ###: ", numSteps);
    print_tagged_stat(out, "Max High-water memory usage ", maxHwm);
}

inline
void print_stats_for_performance_compare(std::ostream& out, double maxTime, size_t maxHwm, unsigned numSteps, MPI_Comm comm)
{
    if (stk::parallel_machine_rank(comm) == 0)
        print_stats(out, maxTime, maxHwm, numSteps);
}

inline
void parallel_print_time_for_performance_compare(MPI_Comm comm, double time_on_this_proc, std::ostream&out = std::cout)
{
    double maxTime = get_max_time_across_procs(time_on_this_proc, comm);
    size_t maxHwm = get_max_hwm_across_procs(comm);
    print_stats_for_performance_compare(out, maxTime, maxHwm, 1, comm);
}

}

#endif /* STK_PERF_UTIL_H */
