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

#include "gtest/gtest.h"
#include "stk_util/environment/memory_util.hpp"  // for get_gpu_memory_info
#include <Kokkos_Core.hpp>
#include <cstddef>                               // for size_t
#include "stk_util/ngp/NgpSpaces.hpp"

#if defined (KOKKOS_ENABLE_CUDA) || defined (KOKKOS_ENABLE_HIP)
TEST(GPUMemoryInfo, singleAllocationOnGPU)
{
  size_t initialUsed, initialFree, initialTotal;
  size_t finalUsed, finalFree, finalTotal;

  stk::get_gpu_memory_info(initialUsed, initialFree);
  initialTotal = initialUsed + initialFree;

  size_t allocationSize = 8000;
  void* allocatedMemory = nullptr;

  do {
    Kokkos::kokkos_free<stk::ngp::MemSpace>(allocatedMemory);
    allocationSize *= 10;
    allocatedMemory = Kokkos::kokkos_malloc<stk::ngp::MemSpace>(allocationSize);

    stk::get_gpu_memory_info(finalUsed, finalFree);
    finalTotal = finalUsed + finalFree;
  } while (finalUsed - initialUsed == 0u);

  EXPECT_GT(initialUsed, 0u);
  EXPECT_GT(initialFree, 0u);
  EXPECT_GT(finalUsed, 0u);
  EXPECT_GT(finalFree, 0u);

  EXPECT_EQ(initialTotal, finalTotal);
  EXPECT_EQ(finalUsed - initialUsed, initialFree - finalFree);

  Kokkos::kokkos_free<stk::ngp::MemSpace>(allocatedMemory);
}
#else
TEST(GPUMemoryInfo, alwaysZeroOnCPU)
{
  size_t used;
  size_t free;

  stk::get_gpu_memory_info(used, free);
  EXPECT_EQ(used, 0u);
  ASSERT_EQ(free, 0u);
}
#endif
