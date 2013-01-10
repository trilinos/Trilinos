#ifndef STK_PERF_UTIL_H
#define STK_PERF_UTIL_H

#include <string>
#include <algorithm>
#include <cassert>
#include <iostream>

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

// Use this version if providing your own memory numbers
template <typename NameItr, typename TimerItr, typename MemItr, typename MemNameItr>
inline
void print_timers_and_memory(NameItr name_itr, TimerItr timer_itr, int timer_count, MemNameItr mem_name_itr, MemItr mem_itr, int mem_count)
{
  print_timers(name_itr, timer_itr, timer_count);

  print_timers(mem_name_itr, mem_itr, mem_count);
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

}

#endif /* STK_PERF_UTIL_H */
