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
#include "stk_util/util/StridedArray.hpp"
#include <vector>

TEST( StridedArray, ptr_and_size)
{
  std::vector<int> vec = {1, 3, 5, 2};
  stk::util::StridedArray<int> stridedArray(vec.data(), vec.size());

  EXPECT_EQ(vec.size(), stridedArray.size());

  for(unsigned i=0; i<stridedArray.size(); ++i) {
    EXPECT_EQ(vec[i],  stridedArray[i]);
  }
}

TEST( StridedArray, pair_iter)
{
  std::vector<int> vec = {1, 3, 5, 2};
  stk::util::StridedArray<int> stridedArray(stk::PairIter<int*>(vec.data(), vec.data()+vec.size()));

  EXPECT_EQ(vec.size(), stridedArray.size());

  for(unsigned i=0; i<stridedArray.size(); ++i) {
    EXPECT_EQ(vec[i],  stridedArray[i]);
  }
}

TEST( StridedArray, comparison)
{
  std::vector<int> vec = {1, 2, 3, 4};
  stk::util::StridedArray<int> stridedArray1(vec.data(), vec.size());
  stk::util::StridedArray<int> stridedArray2(vec.data(), vec.size());

  EXPECT_EQ(stridedArray1, stridedArray2);

  std::vector<int> otherVec = {1, 4, 3, 2};
  stk::util::StridedArray<int> otherStridedArray(otherVec.data(), otherVec.size());

  EXPECT_NE(stridedArray1, otherStridedArray);

  std::vector<int> smallerVec = {1, 2};
  stk::util::StridedArray<int> smallerStridedArray(smallerVec.data(), smallerVec.size());

  EXPECT_NE(stridedArray1, smallerStridedArray);
}

#ifdef STK_ENABLE_GPU
void run_comparison_test_on_device()
{
  constexpr int N1 = 3;
  constexpr int N2 = 5;
  Kokkos::View<int**> v1("v1", N1, N2);

  Kokkos::parallel_for(N1, KOKKOS_LAMBDA(const int& i) {
    for(int j=0; j<N2; ++j) {
      v1(i,j) = 1 + i + j;
    }
  });

  int result = 0;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const int& i, int& localResult) {
    int stride = N1;
    int wrongStride = 2;
    stk::util::StridedArray<int> sa(v1.data(), N2, stride);
    localResult = static_cast<int>(sa.size())==N2 ? 0 : 1;
    localResult += sa[1]==v1(0,1) ? 0 : 3;
    stk::util::StridedArray<int> sa2(v1.data(), N2, stride);
    localResult += (sa == sa2) ? 0 : 7;

    stk::util::StridedArray<int> saWrong(v1.data(), N2, wrongStride);
    localResult += (sa != saWrong) ? 0 : 15;
  }, result);

  EXPECT_EQ(0, result);
}

TEST(StridedArray, comparison_on_device)
{
  run_comparison_test_on_device();
}

#endif

