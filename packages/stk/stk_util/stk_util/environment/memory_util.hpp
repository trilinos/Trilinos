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

#ifndef STK_MEMORY_UTIL_H
#define STK_MEMORY_UTIL_H

#include <vector>
#include <string>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk {

  /*
   *  return current memory usage ('now', resident set size) and high-water-mark ('hwm') in bytes.
   *
   *  Note: Only systems with PROCFS (i.e., Linux) report the high-water-mark of the
   *  resident set size for a process. For those systems, *hwm may be larger than any
   *  observed value of *now. On all other systems, *hwm is the maximum value of *now
   *  that was observed through invoking this function so the true high-water-mark may
   *  be larger than that reported here, depending on your timing.
   */
  size_t get_memory_usage_now();
  void get_gpu_memory_info(size_t& used, size_t& free);
  void get_memory_usage(size_t & now, size_t & hwm);
  void get_memory_available(size_t & avail);
  void get_processor_count(std::vector<int> &procinfo);
  void get_memory_high_water_mark_across_processors(MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg);
  void get_memory_available_across_processors(MPI_Comm comm, size_t& avail_max, size_t& avail_min, size_t& avail_avg);
  void get_current_memory_usage_across_processors(MPI_Comm comm, size_t& curr_max, size_t& curr_min, size_t& curr_avg);

  template <typename T>
  inline
  void get_max_min_avg(MPI_Comm comm, T this_proc, T& max, T& min, T& avg)
  {
    int num_procs = stk::parallel_machine_size(comm);

    stk::all_reduce_max(comm, &this_proc, &max, 1);
    stk::all_reduce_min(comm, &this_proc, &min, 1);
    stk::all_reduce_sum(comm, &this_proc, &avg, 1);
    avg /= static_cast<T>(num_procs);
  }
}
#endif /* STK_MEMORY_UTIL_H */
