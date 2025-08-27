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

#include <cstddef>

#include "gtest/gtest.h"
#include "stk_util/util/AlignedAllocator.hpp"

TEST(HostAlignedAllocator, address_divisible_by_alignment)
{
  constexpr size_t alignment = 32;
  non_std::HostAlignedAllocator<std::byte, alignment> alloc;
  auto* ptr = alloc.allocate(alignment);
  auto addr = reinterpret_cast<uintptr_t>(ptr);
  EXPECT_EQ(0U, addr % alignment);
  alloc.deallocate(ptr, 0);
}

TEST(HostAlignedAllocator, allocate_pads_requested_size_to_multiple_of_alignment)
{
  // according to [cppreference](https://en.cppreference.com/w/cpp/memory/c/aligned_alloc),
  // it is up to the implementation on whether the requested allocation must be a strict
  // multiple of the alignment.
  //
  // At least with address sanitizer included in llvm-14.0.6, an error is emitted if the
  // requested allocation is not a multiple of the alignment. Note that this is an ASAN error
  // and **not** an exception. The padding is broken if an ASAN build fails this test.
  constexpr size_t alignment = 32;
  non_std::HostAlignedAllocator<std::byte, alignment> alloc;
  std::byte* ptr = nullptr;
  EXPECT_NO_THROW(ptr = alloc.allocate(7));
  alloc.deallocate(ptr, 0);
}

// TEST(HostAlignedAllocator, alignment_of_non_pow2_is_compile_error)
// {
//   non_std::HostAlignedAllocator<std::byte, 7> alloc;
// }
