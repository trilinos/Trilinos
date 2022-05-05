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

/// \file Test_Common_Transpose.hpp

#ifndef KOKKOSKERNELS_TRANSPOSE_HPP
#define KOKKOSKERNELS_TRANSPOSE_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <KokkosKernels_SparseUtils.hpp>
#include <KokkosKernels_Sorting.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

template <typename size_type, typename V>
struct ExactCompare {
  ExactCompare(const V& v1_, const V& v2_) : v1(v1_), v2(v2_) {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type i, size_type& ldiffs) const {
    if (v1(i) != v2(i)) ldiffs++;
  }

  V v1;
  V v2;
};

template <typename exec_space>
void testTranspose(int numRows, int numCols, bool doValues) {
  using range_pol  = Kokkos::RangePolicy<exec_space>;
  using scalar_t   = default_scalar;
  using lno_t      = default_lno_t;
  using size_type  = default_size_type;
  using mem_space  = typename exec_space::memory_space;
  using device_t   = Kokkos::Device<exec_space, mem_space>;
  using crsMat_t   = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t,
                                                    void, size_type>;
  using c_rowmap_t = typename crsMat_t::row_map_type;
  using c_entries_t = typename crsMat_t::index_type;
  using c_values_t  = typename crsMat_t::values_type;
  using rowmap_t    = typename crsMat_t::row_map_type::non_const_type;
  using entries_t   = typename crsMat_t::index_type::non_const_type;
  using values_t    = typename crsMat_t::values_type::non_const_type;
  size_type nnz     = 10 * numRows;
  // Generate a matrix that has 0 entries in some rows
  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, 3 * 10, numRows / 2);
  // compute the transpose while unsorted, then transpose again
  rowmap_t t_rowmap("Rowmap^T", numCols + 1);  // this view is initialized to 0
  entries_t t_entries(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries^T"),
      input_mat.graph.entries.extent(0));
  values_t t_values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values^T"),
                    input_mat.values.extent(0));
  rowmap_t tt_rowmap("Rowmap^T^T",
                     numRows + 1);  // this view is initialized to 0
  entries_t tt_entries(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries^T^T"),
      input_mat.graph.entries.extent(0));
  values_t tt_values(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values^T"),
      input_mat.values.extent(0));
  if (doValues) {
    KokkosKernels::Impl::transpose_matrix<c_rowmap_t, c_entries_t, c_values_t,
                                          rowmap_t, entries_t, values_t,
                                          rowmap_t, exec_space>(
        numRows, numCols, input_mat.graph.row_map, input_mat.graph.entries,
        input_mat.values, t_rowmap, t_entries, t_values);
    KokkosKernels::Impl::transpose_matrix<rowmap_t, entries_t, values_t,
                                          rowmap_t, entries_t, values_t,
                                          rowmap_t, exec_space>(
        numCols, numRows, t_rowmap, t_entries, t_values, tt_rowmap, tt_entries,
        tt_values);
  } else {
    KokkosKernels::Impl::transpose_graph<c_rowmap_t, c_entries_t, rowmap_t,
                                         entries_t, rowmap_t, exec_space>(
        numRows, numCols, input_mat.graph.row_map, input_mat.graph.entries,
        t_rowmap, t_entries);
    KokkosKernels::Impl::transpose_graph<rowmap_t, entries_t, rowmap_t,
                                         entries_t, rowmap_t, exec_space>(
        numCols, numRows, t_rowmap, t_entries, tt_rowmap, tt_entries);
  }
  // Sort both the transpose-transpose, and the original matrix (to compare
  // directly)
  KokkosKernels::sort_crs_matrix(input_mat);
  KokkosKernels::sort_crs_matrix<exec_space, c_rowmap_t, entries_t, values_t>(
      tt_rowmap, tt_entries, tt_values);
  // The views should now be exactly identical, since they represent the same
  // matrix and are sorted
  size_type rowmapDiffs;
  Kokkos::parallel_reduce(
      range_pol(0, numRows + 1),
      ExactCompare<size_type, c_rowmap_t>(input_mat.graph.row_map, tt_rowmap),
      rowmapDiffs);
  size_type entriesDiffs;
  Kokkos::parallel_reduce(
      range_pol(0, input_mat.nnz()),
      ExactCompare<size_type, c_entries_t>(input_mat.graph.entries, tt_entries),
      entriesDiffs);
  EXPECT_EQ(size_type(0), rowmapDiffs);
  EXPECT_EQ(size_type(0), entriesDiffs);
  if (doValues) {
    size_type valuesDiffs;
    Kokkos::parallel_reduce(
        range_pol(0, input_mat.nnz()),
        ExactCompare<size_type, values_t>(input_mat.values, tt_values),
        valuesDiffs);
    EXPECT_EQ(size_type(0), valuesDiffs);
  }
}

TEST_F(TestCategory, common_transpose_matrix) {
  // Test both matrix and graph transpose with various sizes
  testTranspose<TestExecSpace>(100, 100, true);
  testTranspose<TestExecSpace>(500, 50, true);
  testTranspose<TestExecSpace>(50, 500, true);
  testTranspose<TestExecSpace>(4000, 2000, true);
  testTranspose<TestExecSpace>(2000, 4000, true);
  testTranspose<TestExecSpace>(2000, 2000, true);
}

TEST_F(TestCategory, common_transpose_graph) {
  testTranspose<TestExecSpace>(100, 100, false);
  testTranspose<TestExecSpace>(500, 50, false);
  testTranspose<TestExecSpace>(50, 500, false);
  testTranspose<TestExecSpace>(4000, 2000, false);
  testTranspose<TestExecSpace>(2000, 4000, false);
  testTranspose<TestExecSpace>(2000, 2000, false);
}

#endif
