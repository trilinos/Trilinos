// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
