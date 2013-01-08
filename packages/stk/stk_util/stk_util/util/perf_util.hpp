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

// Prefer this version unless you are sure you are not interested in memory.
template <typename NameItr, typename TimerItr>
inline
void print_timers_and_memory(NameItr name_itr, TimerItr timer_itr, int count)
{
  print_timers(name_itr, timer_itr, count);

  size_t mem_now = 0, mem_hwm = 0;
  stk::get_memory_usage(mem_now, mem_hwm);

  std::cout << "STKPERF: Current memory: " << mem_now << std::endl;
  std::cout << "STKPERF: Memory high water: " << mem_hwm << std::endl;
}

}

#endif /* STK_PERF_UTIL_H */
