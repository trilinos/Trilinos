/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include "KokkosKernels_TestUtils.hpp"

namespace Test {
template <class ScalarType, class LayoutType, class ExeSpaceType>
void doCscMat(size_t m, size_t n, ScalarType min_val, ScalarType max_val) {
  auto expected_min    = ScalarType(1.0);
  int64_t expected_nnz = 0;
  RandCscMat<ScalarType, LayoutType, ExeSpaceType> cm(m, n, min_val, max_val);

  for (int64_t i = 0; i < cm.get_nnz(); ++i)
    ASSERT_GE(cm(i), expected_min) << cm.info;

  for (int64_t j = 0; j < cm.get_n(); ++j) {
    for (int64_t i = 0; i < cm.get_col_len(j); ++i)
      ASSERT_FLOAT_EQ(cm(cm.get_col_start(j) + i), cm(expected_nnz + i))
          << cm.info;
    expected_nnz += cm.get_col_len(j);
  }
  ASSERT_EQ(cm.get_nnz(), expected_nnz) << cm.info;

  // No need to check data here. Kokkos unit-tests deep_copy.
  auto vals = cm.get_vals();
  ASSERT_EQ(vals.extent(0), cm.get_nnz() + 1) << cm.info;

  auto row_ids = cm.get_row_ids();
  ASSERT_EQ(row_ids.extent(0), cm.get_n() * cm.get_m() + 1) << cm.info;

  auto col_map = cm.get_col_map();
  ASSERT_EQ(col_map.extent(0), cm.get_n() + 1);
}

template <class ExeSpaceType>
void doAllCscMat(size_t m, size_t n) {
  int min = 1, max = 10;

  // Verify that CscMax is constructed properly.
  doCscMat<float, Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doCscMat<float, Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);

  doCscMat<double, Kokkos::LayoutLeft, ExeSpaceType>(m, n, min, max);
  doCscMat<double, Kokkos::LayoutRight, ExeSpaceType>(m, n, min, max);

  // Verify that CscMax can be instantiated with complex types.
  RandCscMat<Kokkos::complex<float>, Kokkos::LayoutLeft, ExeSpaceType> cmcf(
      m, n, min, max);
  RandCscMat<Kokkos::complex<double>, Kokkos::LayoutRight, ExeSpaceType> cmcd(
      m, n, min, max);
}

// Test randomly generated csc matrices
TEST_F(TestCategory, sparse_randcscmat) {
  // Square cases
  for (int dim = 1; dim < 1024; dim *= 4) doAllCscMat<TestExecSpace>(dim, dim);

  // Non-square cases
  for (int dim = 1; dim < 1024; dim *= 4) {
    doAllCscMat<TestExecSpace>(dim * 3, dim);
    doAllCscMat<TestExecSpace>(dim, dim * 3);
  }
}
}  // namespace Test