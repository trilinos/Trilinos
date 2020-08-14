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

#include <stddef.h>                     // for size_t, ptrdiff_t
#include <stk_util/environment/memory_util.hpp>  // for cpu_time
#include <gtest/gtest.h>
#include <iostream>

//inline T* kokkos_malloc_on_device(const std::string& debuggingName) 
//{
//  return static_cast<T*>(Kokkos::kokkos_malloc<Kokkos::CudaSpace>(debuggingName, sizeof(T)));
//}
//
//inline void kokkos_free_on_device(void* ptr) 
//{ 
//  Kokkos::kokkos_free<Kokkos::CudaSpace>(ptr); 
//}


TEST(GPUMemoryInfo, alwaysZeroOnCPU)
{
#ifndef __CUDACC__
  size_t used;
  size_t free;

  stk::get_gpu_memory_info(used, free);
  EXPECT_EQ(used, 0u);
  ASSERT_EQ(free, 0u);
#endif
}

TEST(GPUMemoryInfo, singleAllocationOnGPU)
{
#ifdef __CUDACC__
  size_t initialUsed, initialFree, initialTotal;
  size_t finalUsed, finalFree, finalTotal;

  constexpr size_t allocationSize = sizeof(double) * 100000000;

  cudaSetDevice(0);

  stk::get_gpu_memory_info(initialUsed, initialFree);
  initialTotal = initialUsed + initialFree;
  std::cout << "InitialUsed: " << initialUsed << "\tInitialFree: " << initialFree << "\tInitialTotal: " << initialTotal << std::endl;

  double* allocatedMemory = static_cast<double*>(Kokkos::kokkos_malloc<Kokkos::CudaSpace>(allocationSize));

  stk::get_gpu_memory_info(finalUsed, finalFree);
  finalTotal = finalUsed + finalFree;
  std::cout << "FinalUsed: " << finalUsed << "\tFinalFree: " << finalFree << "\tFinalTotal: " << finalTotal << std::endl;

  EXPECT_EQ(initialTotal, finalTotal);

  EXPECT_GT(finalUsed - initialUsed, allocationSize);
  EXPECT_GT(initialFree - finalFree, allocationSize);
  EXPECT_EQ(finalUsed - initialUsed, initialFree - finalFree);

  Kokkos::kokkos_free<Kokkos::CudaSpace>(allocatedMemory);
#endif
}
