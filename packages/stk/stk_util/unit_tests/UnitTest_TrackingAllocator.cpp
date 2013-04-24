#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/TrackingAllocator.hpp>

#include <iostream>
#include <vector>

STKUNIT_UNIT_TEST( tracking_allocator, vector )
{
  typedef stk::tracking_allocator<int> allocator;

  EXPECT_EQ( allocator::peak_memory(), 0u);
  EXPECT_EQ( allocator::current_memory(), 0u);
  EXPECT_EQ( allocator::num_allocations(), 0u);

  {
    std::vector<int,allocator> vec;

    EXPECT_EQ( allocator::peak_memory(), 0u);
    EXPECT_EQ( allocator::current_memory(), 0u);
    EXPECT_EQ( allocator::num_allocations(), 0u);
  }


  {
    std::vector<int,allocator> vec(1000);

    EXPECT_EQ( allocator::peak_memory(), sizeof(int)*1000u);
    EXPECT_EQ( allocator::current_memory(), sizeof(int)*1000u);
    EXPECT_EQ( allocator::num_allocations(), 1u);
  }

  {
    std::vector<int,allocator> vec(100);

    EXPECT_EQ( allocator::peak_memory(), sizeof(int)*1000u);
    EXPECT_EQ( allocator::current_memory(), sizeof(int)*100u);
    EXPECT_EQ( allocator::num_allocations(), 2u);
  }

  {
    std::vector<int,allocator> vec(1000);

    EXPECT_EQ( allocator::peak_memory(), sizeof(int)*1000u);
    EXPECT_EQ( allocator::current_memory(), sizeof(int)*1000u);
    EXPECT_EQ( allocator::num_allocations(), 3u);

    vec.resize(2000);

    size_t capacity = vec.capacity();
    EXPECT_EQ( allocator::peak_memory(), sizeof(int)*(1000u+capacity));
    EXPECT_EQ( allocator::current_memory(), sizeof(int)*capacity);
    EXPECT_EQ( allocator::num_allocations(), 4u);

    {
      std::vector<int,allocator> tmp;
      vec.swap(tmp);
    }

    EXPECT_EQ( allocator::peak_memory(), sizeof(int)*(1000u+capacity));
    EXPECT_EQ( allocator::current_memory(), 0u);
    EXPECT_EQ( allocator::num_allocations(), 4u);
  }
}
