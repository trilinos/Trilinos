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

/// \file Test_Sparse_SortCrs.hpp
/// \brief Tests for sort_crs_matrix and sort_crs_graph in
/// KokkosSparse_SortCrs.hpp

#ifndef KOKKOSSPARSE_SORTCRSTEST_HPP
#define KOKKOSSPARSE_SORTCRSTEST_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <KokkosKernels_Utils.hpp>
#include "KokkosSparse_IOUtils.hpp"
#include <KokkosSparse_SortCrs.hpp>
#include <KokkosKernels_default_types.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_Complex.hpp>
#include <cstdlib>

template <typename exec_space>
void testSortCRS(default_lno_t numRows, default_lno_t numCols,
                 default_size_type nnz, bool doValues, bool doStructInterface) {
  using scalar_t  = default_scalar;
  using lno_t     = default_lno_t;
  using size_type = default_size_type;
  using mem_space = typename exec_space::memory_space;
  using device_t  = Kokkos::Device<exec_space, mem_space>;
  using crsMat_t =
      KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type>;
  using rowmap_t  = typename crsMat_t::row_map_type;
  using entries_t = typename crsMat_t::index_type;
  using values_t  = typename crsMat_t::values_type;
  // Create a random matrix on device
  // IMPORTANT: kk_generate_sparse_matrix does not sort the rows, if it did this
  // wouldn't test anything
  crsMat_t A = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(
      numRows, numCols, nnz, 2, numCols / 2);
  auto rowmap  = A.graph.row_map;
  auto entries = A.graph.entries;
  auto values  = A.values;
  Kokkos::View<size_type*, Kokkos::HostSpace> rowmapHost("rowmap host",
                                                         numRows + 1);
  Kokkos::View<lno_t*, Kokkos::HostSpace> entriesHost("sorted entries host",
                                                      nnz);
  Kokkos::View<scalar_t*, Kokkos::HostSpace> valuesHost("sorted values host",
                                                        nnz);
  Kokkos::deep_copy(rowmapHost, rowmap);
  Kokkos::deep_copy(entriesHost, entries);
  Kokkos::deep_copy(valuesHost, values);
  struct ColValue {
    ColValue() {}
    ColValue(lno_t c, scalar_t v) : col(c), val(v) {}
    bool operator<(const ColValue& rhs) const { return col < rhs.col; }
    bool operator==(const ColValue& rhs) const {
      return col == rhs.col && val == rhs.val;
    }
    lno_t col;
    scalar_t val;
  };
  // sort one row at a time on host using STL.
  {
    for (lno_t i = 0; i < numRows; i++) {
      std::vector<ColValue> rowCopy;
      for (size_type j = rowmapHost(i); j < rowmapHost(i + 1); j++)
        rowCopy.emplace_back(entriesHost(j), valuesHost(j));
      std::sort(rowCopy.begin(), rowCopy.end());
      // write sorted row back
      for (size_t j = 0; j < rowCopy.size(); j++) {
        entriesHost(rowmapHost(i) + j) = rowCopy[j].col;
        valuesHost(rowmapHost(i) + j)  = rowCopy[j].val;
      }
    }
  }
  // call the actual sort routine being tested
  if (doValues) {
    if (doStructInterface) {
      KokkosSparse::sort_crs_matrix(A);
    } else {
      KokkosSparse::sort_crs_matrix<exec_space, rowmap_t, entries_t, values_t>(
          A.graph.row_map, A.graph.entries, A.values);
    }
  } else {
    if (doStructInterface) {
      KokkosSparse::sort_crs_graph(A.graph);
    } else {
      KokkosSparse::sort_crs_graph<exec_space, rowmap_t, entries_t>(
          A.graph.row_map, A.graph.entries);
    }
  }
  // Copy to host and compare
  Kokkos::View<lno_t*, Kokkos::HostSpace> entriesOut("sorted entries host",
                                                     nnz);
  Kokkos::View<scalar_t*, Kokkos::HostSpace> valuesOut("sorted values host",
                                                       nnz);
  Kokkos::deep_copy(entriesOut, entries);
  Kokkos::deep_copy(valuesOut, values);
  for (size_type i = 0; i < nnz; i++) {
    EXPECT_EQ(entriesHost(i), entriesOut(i))
        << "Sorted column indices are wrong!";
    if (doValues) {
      EXPECT_EQ(valuesHost(i), valuesOut(i)) << "Sorted values are wrong!";
    }
  }
}

template <typename exec_space>
void testSortCRSUnmanaged(bool doValues, bool doStructInterface) {
  // This test is about bug #960.
  using scalar_t  = default_scalar;
  using lno_t     = default_lno_t;
  using size_type = default_size_type;
  using mem_space = typename exec_space::memory_space;
  using device_t  = Kokkos::Device<exec_space, mem_space>;
  using crsMat_t =
      KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                              size_type>;
  using crsMat_Managed_t =
      KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type>;
  using rowmap_t      = typename crsMat_t::row_map_type;
  using entries_t     = typename crsMat_t::index_type;
  using values_t      = typename crsMat_t::values_type;
  const lno_t numRows = 50;
  const lno_t numCols = numRows;
  size_type nnz       = numRows * 5;
  // Create a random matrix on device
  // IMPORTANT: kk_generate_sparse_matrix does not sort the rows, if it did this
  // wouldn't test anything
  crsMat_Managed_t A_managed =
      KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_Managed_t>(
          numRows, numCols, nnz, 2, numCols / 2);
  crsMat_t A(A_managed);
  auto rowmap  = A.graph.row_map;
  auto entries = A.graph.entries;
  auto values  = A.values;
  if (doValues) {
    if (doStructInterface) {
      KokkosSparse::sort_crs_matrix(A);
    } else {
      KokkosSparse::sort_crs_matrix<exec_space, rowmap_t, entries_t, values_t>(
          A.graph.row_map, A.graph.entries, A.values);
    }
  } else {
    if (doStructInterface) {
      KokkosSparse::sort_crs_graph(A.graph);
    } else {
      KokkosSparse::sort_crs_graph<exec_space, rowmap_t, entries_t>(
          A.graph.row_map, A.graph.entries);
    }
  }
}

template <typename exec_space>
void testSortAndMerge() {
  using size_type = default_size_type;
  using lno_t     = default_lno_t;
  using scalar_t  = default_scalar;
  using mem_space = typename exec_space::memory_space;
  using device_t  = Kokkos::Device<exec_space, mem_space>;
  using crsMat_t =
      KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type>;
  using rowmap_t  = typename crsMat_t::row_map_type::non_const_type;
  using entries_t = typename crsMat_t::index_type;
  using values_t  = typename crsMat_t::values_type;
  using Kokkos::HostSpace;
  using Kokkos::MemoryTraits;
  using Kokkos::Unmanaged;
  // Create a small CRS matrix on host
  std::vector<size_type> inRowmap = {0, 4, 4, 5, 7, 10};
  std::vector<lno_t> inEntries    = {
      4, 3, 5, 3,  // row 0
                   // row 1 has no entries
      6,           // row 2
      2, 2,        // row 3
      0, 1, 2      // row 4
  };
  // note: choosing values that can be represented exactly by float
  std::vector<scalar_t> inValues = {
      1.5, 4, 1, -3,  // row 0
                      // row 1
      2,              // row 2
      -1, -2,         // row 3
      0, 3.5, -2.25   // row 4
  };
  lno_t nrows   = 5;
  lno_t ncols   = 7;
  size_type nnz = inEntries.size();
  Kokkos::View<size_type*, HostSpace, MemoryTraits<Unmanaged>> hostInRowmap(
      inRowmap.data(), nrows + 1);
  Kokkos::View<lno_t*, HostSpace, MemoryTraits<Unmanaged>> hostInEntries(
      inEntries.data(), nnz);
  Kokkos::View<scalar_t*, HostSpace, MemoryTraits<Unmanaged>> hostInValues(
      inValues.data(), nnz);
  rowmap_t devInRowmap("", nrows + 1);
  entries_t devInEntries("", nnz);
  values_t devInValues("", nnz);
  Kokkos::deep_copy(devInRowmap, hostInRowmap);
  Kokkos::deep_copy(devInEntries, hostInEntries);
  Kokkos::deep_copy(devInValues, hostInValues);
  crsMat_t input("Input", nrows, ncols, nnz, devInValues, devInRowmap,
                 devInEntries);
  crsMat_t output = KokkosSparse::sort_and_merge_matrix(input);
  exec_space().fence();
  EXPECT_EQ(output.numRows(), nrows);
  EXPECT_EQ(output.numCols(), ncols);
  auto outRowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                       output.graph.row_map);
  auto outEntries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                        output.graph.entries);
  auto outValues =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output.values);
  // Expect 2 merges to have taken place
  std::vector<size_type> goldRowmap = {0, 3, 3, 4, 5, 8};
  std::vector<lno_t> goldEntries    = {
      3, 4, 5,  // row 0
                // row 1 has no entries
      6,        // row 2
      2,        // row 3
      0, 1, 2   // row 4
  };
  // note: choosing values that can be represented exactly by float
  std::vector<scalar_t> goldValues = {
      1, 1.5, 1,     // row 0
                     // row 1
      2,             // row 2
      -3,            // row 3
      0, 3.5, -2.25  // row 4
  };
  EXPECT_EQ(goldRowmap.size(), outRowmap.extent(0));
  EXPECT_EQ(goldEntries.size(), outEntries.extent(0));
  EXPECT_EQ(goldValues.size(), outValues.extent(0));
  EXPECT_EQ(goldValues.size(), output.nnz());
  for (lno_t i = 0; i < nrows + 1; i++) EXPECT_EQ(goldRowmap[i], outRowmap(i));
  for (size_type i = 0; i < output.nnz(); i++) {
    EXPECT_EQ(goldEntries[i], outEntries(i));
    EXPECT_EQ(goldValues[i], outValues(i));
  }
}

TEST_F(TestCategory, common_sort_crsgraph) {
  for (int doStructInterface = 0; doStructInterface < 2; doStructInterface++) {
    testSortCRS<TestExecSpace>(10, 10, 20, false, doStructInterface);
    testSortCRS<TestExecSpace>(100, 100, 2000, false, doStructInterface);
    testSortCRS<TestExecSpace>(1000, 1000, 30000, false, doStructInterface);
    testSortCRSUnmanaged<TestExecSpace>(false, doStructInterface);
  }
}

TEST_F(TestCategory, common_sort_crsmatrix) {
  for (int doStructInterface = 0; doStructInterface < 2; doStructInterface++) {
    testSortCRS<TestExecSpace>(10, 10, 20, true, doStructInterface);
    testSortCRS<TestExecSpace>(100, 100, 2000, true, doStructInterface);
    testSortCRS<TestExecSpace>(1000, 1000, 30000, true, doStructInterface);
    testSortCRSUnmanaged<TestExecSpace>(true, doStructInterface);
  }
}

TEST_F(TestCategory, common_sort_crs_longrows) {
  testSortCRS<TestExecSpace>(1, 50000, 10000, false, false);
  testSortCRS<TestExecSpace>(1, 50000, 10000, true, false);
}

TEST_F(TestCategory, common_sort_merge_crsmatrix) {
  testSortAndMerge<TestExecSpace>();
}

#endif  // KOKKOSSPARSE_SORTCRSTEST_HPP
