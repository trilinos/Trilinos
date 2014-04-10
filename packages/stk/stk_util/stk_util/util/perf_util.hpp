#ifndef STK_PERF_UTIL_H
#define STK_PERF_UTIL_H

#include <string>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <mpi.h>
#include <iomanip>

#include <stk_util/util/memory_util.hpp>

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

  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  const double bytes_in_MB = 1024*1024;

  if (rank == 0) {
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Min High-water memory usage " << hwm_min / bytes_in_MB << " MB" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Avg High-water memory usage " << hwm_avg / bytes_in_MB << " MB" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Max High-water memory usage " << hwm_max / bytes_in_MB << " MB\n" << std::endl;

    out << std::setw(6) << std::fixed << std::setprecision(1) << "Min Available memory per processor " << avail_min / bytes_in_MB << " MB" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Avg Available memory per processor " << avail_avg / bytes_in_MB << " MB" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Max Available memory per processor " << avail_max / bytes_in_MB << " MB\n" << std::endl;

    out << std::setw(6) << std::fixed << std::setprecision(1) << "Min No-output time " << min_time << " sec" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Avg No-output time " << avg_time << " sec" << std::endl;
    out << std::setw(6) << std::fixed << std::setprecision(1) << "Max No-output time " << max_time << " sec" << std::endl;
  }
}

}

#endif /* STK_PERF_UTIL_H */
