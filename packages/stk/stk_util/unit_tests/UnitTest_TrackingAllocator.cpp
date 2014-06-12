
#include <stddef.h>                     // for size_t
#include <gtest/gtest.h>
#include <stk_util/util/TrackingAllocator.hpp>  // for tracking_allocator
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ
#include "stk_util/util/AllocatorMemoryUsage.hpp"

TEST( tracking_allocator, vector )
{
#ifndef __IBMCPP__
  typedef stk::tracking_allocator<int> allocator;
  typedef stk::allocator_memory_usage<void> usage;

  usage::reset();

  EXPECT_EQ( usage::peak_memory(), 0u);
  EXPECT_EQ( usage::current_memory(), 0u);
  EXPECT_EQ( usage::num_allocations(), 0u);

  {
    std::vector<int,allocator> vec;

    EXPECT_EQ( usage::peak_memory(), 0u);
    EXPECT_EQ( usage::current_memory(), 0u);
    EXPECT_EQ( usage::num_allocations(), 0u);
  }


  {
    std::vector<int,allocator> vec(1000);

    EXPECT_EQ( usage::peak_memory(), sizeof(int)*1000u);
    EXPECT_EQ( usage::current_memory(), sizeof(int)*1000u);
    EXPECT_EQ( usage::num_allocations(), 1u);
  }

  {
    std::vector<int,allocator> vec(100);

    EXPECT_EQ( usage::peak_memory(), sizeof(int)*1000u);
    EXPECT_EQ( usage::current_memory(), sizeof(int)*100u);
    EXPECT_EQ( usage::num_allocations(), 2u);
  }

  {
    std::vector<int,allocator> vec(1000);

    EXPECT_EQ( usage::peak_memory(), sizeof(int)*1000u);
    EXPECT_EQ( usage::current_memory(), sizeof(int)*1000u);
    EXPECT_EQ( usage::num_allocations(), 3u);

    vec.resize(2000);

    size_t capacity = vec.capacity();
    EXPECT_EQ( usage::peak_memory(), sizeof(int)*(1000u+capacity));
    EXPECT_EQ( usage::current_memory(), sizeof(int)*capacity);
    EXPECT_EQ( usage::num_allocations(), 4u);

    {
      std::vector<int,allocator> tmp;
      vec.swap(tmp);
    }

    EXPECT_EQ( usage::peak_memory(), sizeof(int)*(1000u+capacity));
    EXPECT_EQ( usage::current_memory(), 0u);
    EXPECT_EQ( usage::num_allocations(), 4u);
  }
#endif
}
