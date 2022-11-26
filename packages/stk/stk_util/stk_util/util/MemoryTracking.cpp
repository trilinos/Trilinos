
#include "stk_util/util/MemoryTracking.hpp"
#include <stddef.h>
#include <cstddef>
#include <cstdlib>
#include <utility>
#include <cstdio>
#include <algorithm>
#include <iostream>

#ifdef STK_MEMORY_TRACKING

static constexpr size_t allocOverhead = 16;
static size_t total_bytes_allocated = 0;
static size_t high_water_mark_bytes_baseline = 0;
static size_t high_water_mark_bytes = 0;
static size_t high_water_mark_in_ptrs_baseline = 0;
static size_t high_water_mark_in_ptrs = 0;

static const size_t max_num_ptrs = 20000000;
static size_t num_ptrs = 0;
static std::pair<void*,size_t> ptr_sizes[max_num_ptrs];

void add_bytes(size_t nbytes)
{
  total_bytes_allocated += nbytes;
  size_t current_water_mark_bytes = total_bytes_allocated > high_water_mark_bytes_baseline ? total_bytes_allocated - high_water_mark_bytes_baseline : 0;
  high_water_mark_bytes = std::max(high_water_mark_bytes, current_water_mark_bytes);
}

void add_ptr(void* ptr, size_t sz)
{
  if (num_ptrs >= max_num_ptrs) {
      std::printf("too many ptrs allocated!");
      std::abort();
  }
  ptr_sizes[num_ptrs++] = std::make_pair(ptr,sz);
  size_t current_water_mark_ptrs = num_ptrs > high_water_mark_in_ptrs_baseline ? num_ptrs - high_water_mark_in_ptrs_baseline : 0;
  high_water_mark_in_ptrs = std::max(high_water_mark_in_ptrs, current_water_mark_ptrs);
}

size_t remove_ptr_and_return_size(void*ptr) {
  int num = num_ptrs;
  for(int i=num-1; i>=0; --i) {
    if (ptr_sizes[i].first == ptr) {
      size_t size = ptr_sizes[i].second;
      ptr_sizes[i] = ptr_sizes[num-1];
      --num_ptrs;
      return size;
    }
  }
  return 0;
}

void* operator new(std::size_t sz) {
  void* ptr = std::malloc(sz);
  add_bytes(sz+allocOverhead);
  add_ptr(ptr,sz);
  return ptr;
}

void* operator new  ( std::size_t sz, std::align_val_t /* al */)
{
  void* ptr = std::malloc(sz);
  add_bytes(sz+allocOverhead);
  add_ptr(ptr,sz);
  return ptr;
}

void operator delete(void* ptr) throw()
{
  total_bytes_allocated -= (remove_ptr_and_return_size(ptr)+allocOverhead);
  std::free(ptr);
}

void operator delete(void* ptr, size_t /*sz*/) throw()
{
  total_bytes_allocated -= (remove_ptr_and_return_size(ptr)+allocOverhead);
  std::free(ptr);
}

void operator delete  ( void* ptr, std::align_val_t /*al*/ ) noexcept
{
  total_bytes_allocated -= (remove_ptr_and_return_size(ptr)+allocOverhead);
  std::free(ptr);
}

namespace stk {

size_t get_total_bytes_currently_allocated()
{
    return total_bytes_allocated;
}

size_t get_current_num_ptrs()
{
  return num_ptrs;
}

size_t get_high_water_mark_in_bytes()
{
    return high_water_mark_bytes;
}

size_t get_high_water_mark_in_ptrs()
{
    return high_water_mark_in_ptrs;
}

void reset_high_water_mark_in_bytes()
{
  high_water_mark_bytes_baseline = high_water_mark_bytes;
  high_water_mark_bytes = 0;
}

void reset_high_water_mark_in_ptrs()
{
  high_water_mark_in_ptrs_baseline = high_water_mark_in_ptrs;
  high_water_mark_in_ptrs = 0;
}

}//namespace stk

#endif

