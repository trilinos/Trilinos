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
  void get_memory_usage(size_t & now, size_t & hwm);
  void get_memory_available(size_t & avail);
  void get_processor_count(std::vector<int> &procinfo);
  void get_memory_high_water_mark_across_processors(MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg);
  void get_memory_available_across_processors(MPI_Comm comm, size_t& avail_max, size_t& avail_min, size_t& avail_avg);

  template <typename T>
  inline
  void get_max_min_avg(MPI_Comm comm, T this_proc, T& max, T& min, T& avg)
  {
    int num_procs = stk::parallel_machine_size(comm);
    T this_proc_average = this_proc / num_procs;

    stk::all_reduce_max(comm, &this_proc, &max, 1);
    stk::all_reduce_min(comm, &this_proc, &min, 1);
    stk::all_reduce_sum(comm, &this_proc_average, &avg, 1);
  }
}
#endif /* STK_MEMORY_UTIL_H */
