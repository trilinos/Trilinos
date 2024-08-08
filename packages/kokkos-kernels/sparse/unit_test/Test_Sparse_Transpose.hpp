//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

/// \file Test_Common_Transpose.hpp

#ifndef KOKKOSKERNELS_TRANSPOSE_HPP
#define KOKKOSKERNELS_TRANSPOSE_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <KokkosSparse_Utils.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>
#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_SortCrs.hpp>

template <typename size_type, typename V>
struct ExactCompare {
  ExactCompare(const V& v1_, const V& v2_) : v1(v1_), v2(v2_) {}

  KOKKOS_INLINE_FUNCTION void operator()(size_type i, size_type& ldiffs) const {
    if (v1(i) != v2(i)) ldiffs++;
  }

  V v1;
  V v2;
};

template <typename device_t>
void testTranspose(int numRows, int numCols, bool doValues) {
  using exec_space  = typename device_t::execution_space;
  using range_pol   = Kokkos::RangePolicy<exec_space>;
  using scalar_t    = default_scalar;
  using lno_t       = default_lno_t;
  using size_type   = default_size_type;
  using crsMat_t    = typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type>;
  using c_rowmap_t  = typename crsMat_t::row_map_type;
  using c_entries_t = typename crsMat_t::index_type;
  using c_values_t  = typename crsMat_t::values_type;
  using rowmap_t    = typename crsMat_t::row_map_type::non_const_type;
  using entries_t   = typename crsMat_t::index_type::non_const_type;
  using values_t    = typename crsMat_t::values_type::non_const_type;
  size_type nnz     = (numRows * numCols > 0) ? 10 * numRows : 0;
  // Generate a matrix that has 0 entries in some rows
  crsMat_t input_mat =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows, numCols, nnz, 3 * 10, numRows / 2);
  // compute the transpose while unsorted, then transpose again
  rowmap_t t_rowmap("Rowmap^T", numCols + 1);  // this view is initialized to 0
  entries_t t_entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries^T"), input_mat.graph.entries.extent(0));
  values_t t_values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values^T"), input_mat.values.extent(0));
  rowmap_t tt_rowmap("Rowmap^T^T",
                     numRows + 1);  // this view is initialized to 0
  entries_t tt_entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries^T^T"),
                       input_mat.graph.entries.extent(0));
  values_t tt_values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values^T"), input_mat.values.extent(0));
  if (doValues) {
    KokkosSparse::Impl::transpose_matrix<c_rowmap_t, c_entries_t, c_values_t, rowmap_t, entries_t, values_t, rowmap_t,
                                         exec_space>(numRows, numCols, input_mat.graph.row_map, input_mat.graph.entries,
                                                     input_mat.values, t_rowmap, t_entries, t_values);
    KokkosSparse::Impl::transpose_matrix<rowmap_t, entries_t, values_t, rowmap_t, entries_t, values_t, rowmap_t,
                                         exec_space>(numCols, numRows, t_rowmap, t_entries, t_values, tt_rowmap,
                                                     tt_entries, tt_values);
  } else {
    KokkosSparse::Impl::transpose_graph<c_rowmap_t, c_entries_t, rowmap_t, entries_t, rowmap_t, exec_space>(
        numRows, numCols, input_mat.graph.row_map, input_mat.graph.entries, t_rowmap, t_entries);
    KokkosSparse::Impl::transpose_graph<rowmap_t, entries_t, rowmap_t, entries_t, rowmap_t, exec_space>(
        numCols, numRows, t_rowmap, t_entries, tt_rowmap, tt_entries);
  }
  // Sort both the transpose-transpose, and the original matrix (to compare
  // directly)
  KokkosSparse::sort_crs_matrix(input_mat);
  KokkosSparse::sort_crs_matrix<exec_space, c_rowmap_t, entries_t, values_t>(tt_rowmap, tt_entries, tt_values);
  // The views should now be exactly identical, since they represent the same
  // matrix and are sorted
  size_type rowmapDiffs;
  Kokkos::parallel_reduce(range_pol(0, numRows + 1),
                          ExactCompare<size_type, c_rowmap_t>(input_mat.graph.row_map, tt_rowmap), rowmapDiffs);
  size_type entriesDiffs;
  Kokkos::parallel_reduce(range_pol(0, input_mat.nnz()),
                          ExactCompare<size_type, c_entries_t>(input_mat.graph.entries, tt_entries), entriesDiffs);
  EXPECT_EQ(size_type(0), rowmapDiffs);
  EXPECT_EQ(size_type(0), entriesDiffs);
  if (doValues) {
    size_type valuesDiffs;
    Kokkos::parallel_reduce(range_pol(0, input_mat.nnz()),
                            ExactCompare<size_type, values_t>(input_mat.values, tt_values), valuesDiffs);
    EXPECT_EQ(size_type(0), valuesDiffs);
  }
}

template <class bsrMat_t>
void CompareBsrMatrices(bsrMat_t& A, bsrMat_t& B) {
  using exec_space  = typename bsrMat_t::execution_space;
  using range_pol   = Kokkos::RangePolicy<exec_space>;
  using size_type   = default_size_type;
  using c_rowmap_t  = typename bsrMat_t::row_map_type;
  using c_entries_t = typename bsrMat_t::index_type;
  using values_t    = typename bsrMat_t::values_type::non_const_type;

  // The views should now be exactly identical, since they represent the same
  // matrix and are sorted

  size_type rowmapDiffs;
  Kokkos::parallel_reduce(range_pol(0, A.numRows() + 1),
                          ExactCompare<size_type, c_rowmap_t>(A.graph.row_map, B.graph.row_map), rowmapDiffs);

  size_type entriesDiffs;
  Kokkos::parallel_reduce(range_pol(0, A.nnz()), ExactCompare<size_type, c_entries_t>(A.graph.entries, B.graph.entries),
                          entriesDiffs);

  EXPECT_EQ(size_type(0), rowmapDiffs);
  EXPECT_EQ(size_type(0), entriesDiffs);

  size_type valuesDiffs;
  Kokkos::parallel_reduce(range_pol(0, A.nnz() * A.blockDim() * A.blockDim()),
                          ExactCompare<size_type, values_t>(A.values, B.values), valuesDiffs);
  EXPECT_EQ(size_type(0), valuesDiffs);
}

template <typename device_t>
void testTransposeBsrRef() {
  using scalar_t  = default_scalar;
  using lno_t     = default_lno_t;
  using size_type = default_size_type;
  using bsrMat_t  = typename KokkosSparse::Experimental::BsrMatrix<scalar_t, lno_t, device_t, void, size_type>;
  using rowmap_t  = typename bsrMat_t::row_map_type::non_const_type;
  using entries_t = typename bsrMat_t::index_type::non_const_type;
  using values_t  = typename bsrMat_t::values_type::non_const_type;

  const int numRows    = 4;
  const int nnz        = 7;
  const int block_size = 2;

  // Coming up with a BsrMatrix
  bsrMat_t A;
  {
    rowmap_t row_map("row map", numRows + 1);
    entries_t entries("entries", nnz);
    values_t values("values", nnz * block_size * block_size);

    const size_type row_mapPtr[] = {0, 2, 3, 5, 7};
    const lno_t entriesPtr[]     = {2, 3, 1, 0, 1, 1, 3};
    const scalar_t valuesPtr[]   = {0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3, 3.0, 3.1,
                                    3.2, 3.3, 4.0, 4.1, 4.2, 4.3, 5.0, 5.1, 5.2, 5.3, 6.0, 6.1, 6.2, 6.3};

    typename rowmap_t::HostMirror::const_type row_map_h(row_mapPtr, numRows + 1);
    typename entries_t::HostMirror::const_type entries_h(entriesPtr, nnz);
    typename values_t::HostMirror::const_type values_h(valuesPtr, nnz * block_size * block_size);

    Kokkos::deep_copy(row_map, row_map_h);
    Kokkos::deep_copy(entries, entries_h);
    Kokkos::deep_copy(values, values_h);

    A = bsrMat_t("A", numRows, numRows, nnz, values, row_map, entries, block_size);
  }

  // Constructing the transpose of A manually
  bsrMat_t At_ref;
  {
    rowmap_t row_map("row map", numRows + 1);
    entries_t entries("entries", nnz);
    values_t values("values", nnz * block_size * block_size);

    const size_type row_mapPtr[] = {0, 1, 4, 5, 7};
    const lno_t entriesPtr[]     = {2, 1, 2, 3, 0, 0, 3};
    const scalar_t valuesPtr[]   = {3.0, 3.2, 3.1, 3.3, 2.0, 2.2, 2.1, 2.3, 4.0, 4.2, 4.1, 4.3, 5.0, 5.2,
                                    5.1, 5.3, 0.0, 0.2, 0.1, 0.3, 1.0, 1.2, 1.1, 1.3, 6.0, 6.2, 6.1, 6.3};

    typename rowmap_t::HostMirror::const_type row_map_h(row_mapPtr, numRows + 1);
    typename entries_t::HostMirror::const_type entries_h(entriesPtr, nnz);
    typename values_t::HostMirror::const_type values_h(valuesPtr, nnz * block_size * block_size);

    Kokkos::deep_copy(row_map, row_map_h);
    Kokkos::deep_copy(entries, entries_h);
    Kokkos::deep_copy(values, values_h);

    At_ref = bsrMat_t("A", numRows, numRows, nnz, values, row_map, entries, block_size);
  }

  bsrMat_t At = KokkosSparse::Impl::transpose_bsr_matrix(A);
  KokkosSparse::sort_bsr_matrix(At);

  CompareBsrMatrices(At, At_ref);
}

template <typename device_t>
void testTransposeBsr(int numRows, int numCols, int blockSize) {
  using scalar_t    = default_scalar;
  using lno_t       = default_lno_t;
  using size_type   = default_size_type;
  using exec_space  = typename device_t::execution_space;
  using bsrMat_t    = typename KokkosSparse::Experimental::BsrMatrix<scalar_t, lno_t, device_t, void, size_type>;
  using c_rowmap_t  = typename bsrMat_t::row_map_type;
  using c_entries_t = typename bsrMat_t::index_type;
  using c_values_t  = typename bsrMat_t::values_type;
  using rowmap_t    = typename bsrMat_t::row_map_type::non_const_type;
  using entries_t   = typename bsrMat_t::index_type::non_const_type;
  using values_t    = typename bsrMat_t::values_type::non_const_type;

  // Generate a matrix that has 0 entries in some rows
  size_type nnz = (numRows * numCols > 0) ? 10 * numRows : 0;
  bsrMat_t A =
      KokkosSparse::Impl::kk_generate_sparse_matrix<bsrMat_t>(blockSize, numRows, numCols, nnz, 3, numRows / 4);

  // compute the transpose while unsorted, then transpose again
  rowmap_t t_rowmap("Rowmap^T", numCols + 1);  // this view is initialized to 0
  entries_t t_entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries^T"), A.graph.entries.extent(0));
  values_t t_values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values^T"), A.values.extent(0));
  rowmap_t tt_rowmap("Rowmap^T^T",
                     numRows + 1);  // this view is initialized to 0
  entries_t tt_entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Entries^T^T"), A.graph.entries.extent(0));
  values_t tt_values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Values^T"), A.values.extent(0));

  KokkosSparse::Impl::transpose_bsr_matrix<c_rowmap_t, c_entries_t, c_values_t, rowmap_t, entries_t, values_t,
                                           exec_space>(numRows, numCols, blockSize, A.graph.row_map, A.graph.entries,
                                                       A.values, t_rowmap, t_entries, t_values);

  KokkosSparse::Impl::transpose_bsr_matrix<rowmap_t, entries_t, values_t, rowmap_t, entries_t, values_t, exec_space>(
      numCols, numRows, blockSize, t_rowmap, t_entries, t_values, tt_rowmap, tt_entries, tt_values);
  bsrMat_t Att("Att", numRows, numCols, nnz, tt_values, tt_rowmap, tt_entries, blockSize);

  // Sort both the transpose-transpose, and the original matrix (to compare
  // directly)
  KokkosSparse::sort_bsr_matrix(A);

  KokkosSparse::sort_bsr_matrix(Att);

  CompareBsrMatrices(A, Att);
}

TEST_F(TestCategory, sparse_transpose_matrix) {
  // Test both matrix and graph transpose with various sizes
  testTranspose<TestDevice>(0, 0, true);
  testTranspose<TestDevice>(100, 0, true);
  testTranspose<TestDevice>(0, 100, true);
  testTranspose<TestDevice>(100, 100, true);
  testTranspose<TestDevice>(500, 50, true);
  testTranspose<TestDevice>(50, 500, true);
  testTranspose<TestDevice>(4000, 2000, true);
  testTranspose<TestDevice>(2000, 4000, true);
  testTranspose<TestDevice>(2000, 2000, true);
}

TEST_F(TestCategory, sparse_transpose_graph) {
  testTranspose<TestDevice>(0, 0, false);
  testTranspose<TestDevice>(100, 0, false);
  testTranspose<TestDevice>(0, 100, false);
  testTranspose<TestDevice>(100, 100, false);
  testTranspose<TestDevice>(500, 50, false);
  testTranspose<TestDevice>(50, 500, false);
  testTranspose<TestDevice>(4000, 2000, false);
  testTranspose<TestDevice>(2000, 4000, false);
  testTranspose<TestDevice>(2000, 2000, false);
}

TEST_F(TestCategory, sparse_transpose_bsr_matrix) {
  testTransposeBsrRef<TestDevice>();
  // Test bsrMatrix transpose with various sizes
  testTransposeBsr<TestDevice>(0, 0, 5);
  testTransposeBsr<TestDevice>(100, 0, 5);
  testTransposeBsr<TestDevice>(0, 100, 5);
  testTransposeBsr<TestDevice>(100, 100, 3);
  testTransposeBsr<TestDevice>(500, 50, 5);
  testTransposeBsr<TestDevice>(50, 500, 16);
  testTransposeBsr<TestDevice>(4000, 2000, 3);
  testTransposeBsr<TestDevice>(2000, 4000, 3);
  testTransposeBsr<TestDevice>(2000, 2000, 5);
}

#endif
