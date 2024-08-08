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

#ifndef KOKKOSSPARSE_REMOVECRSZEROS_HPP
#define KOKKOSSPARSE_REMOVECRSZEROS_HPP

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_Utils.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace TestRemoveCrsMatrixZeros {

// Simple, sequential implementation of zero-removal to compare against
template <typename Matrix>
Matrix removeMatrixZerosReference(const Matrix& A) {
  using Offset     = typename Matrix::non_const_size_type;
  using Ordinal    = typename Matrix::ordinal_type;
  using Scalar     = typename Matrix::value_type;
  using KAT        = Kokkos::ArithTraits<Scalar>;
  auto rowmapHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto entriesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  auto valuesHost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
  // First, create the filtered rowmap (the CrsMatrix constructor taking host
  // pointers does expect rowmap to be in Ordinal)
  Ordinal filteredNNZ                 = 0;
  std::vector<Ordinal> filteredRowmap = {0};  // first row begins at 0
  for (Ordinal i = 0; i < A.numRows(); i++) {
    for (Offset j = rowmapHost(i); j < rowmapHost(i + 1); j++) {
      if (valuesHost(j) != KAT::zero()) {
        filteredNNZ++;
      }
    }
    filteredRowmap.push_back(filteredNNZ);
  }
  // Then allocate and fill in the filtered entries and values
  std::vector<Ordinal> filteredEntries;
  std::vector<Scalar> filteredValues;
  for (Offset i = 0; i < A.nnz(); i++) {
    if (valuesHost(i) != KAT::zero()) {
      filteredEntries.push_back(entriesHost(i));
      filteredValues.push_back(valuesHost(i));
    }
  }
  // Copy all the views back to device and construct matrix
  return Matrix("A filtered", A.numRows(), A.numCols(), filteredNNZ, filteredValues.data(), filteredRowmap.data(),
                filteredEntries.data());
}

template <typename Matrix>
Matrix loadMatrixFromVectors(int numRows, int numCols, const std::vector<int>& rowmapRawInt,
                             const std::vector<int>& entriesRawInt, const std::vector<double>& valuesRawDouble) {
  using Offset  = typename Matrix::non_const_size_type;
  using Ordinal = typename Matrix::ordinal_type;
  using Scalar  = typename Matrix::value_type;
  // The CrsMatrix constructor taking host pointers expects rowmap to be in
  // Ordinal
  std::vector<Ordinal> rowmapRaw;
  std::vector<Ordinal> entriesRaw;
  std::vector<Scalar> valuesRaw;
  for (auto val : rowmapRawInt) rowmapRaw.push_back(val);
  for (auto val : entriesRawInt) entriesRaw.push_back(val);
  for (auto val : valuesRawDouble) valuesRaw.push_back(Scalar(val));
  Offset nnz = rowmapRaw.size() ? rowmapRaw[numRows] : 0;
  return Matrix("A", numRows, numCols, nnz, valuesRaw.data(), rowmapRaw.data(), entriesRaw.data());
}

template <typename Matrix>
void getTestInput(int test, Matrix& A, Matrix& Afiltered_ref) {
  using Offset                = typename Matrix::size_type;
  bool haveHardcodedReference = true;
  switch (test) {
    case 0: {
      // No entries, but nonzero dimensions.
      std::vector<int> rowmap = {0, 0, 0, 0, 0};
      std::vector<int> entries;
      std::vector<double> values;
      A             = loadMatrixFromVectors<Matrix>(4, 4, rowmap, entries, values);
      Afiltered_ref = loadMatrixFromVectors<Matrix>(4, 4, rowmap, entries, values);
      break;
    }
    case 1: {
      // Some empty rows, and some zero values
      std::vector<int> rowmap        = {0, 0, 3, 3, 5};
      std::vector<int> entries       = {0, 1, 3, 1, 2};
      std::vector<double> values     = {1, 3, 0, 0, 2};
      A                              = loadMatrixFromVectors<Matrix>(4, 4, rowmap, entries, values);
      std::vector<int> rowmapFilt    = {0, 0, 2, 2, 3};
      std::vector<int> entriesFilt   = {0, 1, 2};
      std::vector<double> valuesFilt = {1, 3, 2};
      Afiltered_ref                  = loadMatrixFromVectors<Matrix>(4, 4, rowmapFilt, entriesFilt, valuesFilt);
      break;
    }
    case 2: {
      // Zero-row matrix, length-0 rowmap
      typename Matrix::row_map_type rowmap;
      typename Matrix::index_type entries;
      typename Matrix::values_type values;
      A             = Matrix("A empty", 0, 0, 0, values, rowmap, entries);
      Afiltered_ref = A;
      break;
    }
    case 3: {
      // Zero-row matrix, length-1 rowmap
      std::vector<int> rowmap = {0};
      std::vector<int> entries;
      std::vector<double> values;
      A             = loadMatrixFromVectors<Matrix>(0, 0, rowmap, entries, values);
      Afiltered_ref = A;
      break;
    }
    case 4: {
      // A row of all zeros that will be filtered
      std::vector<int> rowmap        = {0, 3, 6};
      std::vector<int> entries       = {0, 1, 2, 3, 4, 5};
      std::vector<double> values     = {0, 0, 0, 1, 1, 1};
      A                              = loadMatrixFromVectors<Matrix>(2, 6, rowmap, entries, values);
      std::vector<int> rowmapFilt    = {0, 0, 3};
      std::vector<int> entriesFilt   = {3, 4, 5};
      std::vector<double> valuesFilt = {1, 1, 1};
      Afiltered_ref                  = loadMatrixFromVectors<Matrix>(2, 6, rowmapFilt, entriesFilt, valuesFilt);
      break;
    }
    case 5: {
      // One zero in each row that will be filtered
      std::vector<int> rowmap        = {0, 2, 4, 7};
      std::vector<int> entries       = {0, 1, 1, 2, 0, 1, 2};
      std::vector<double> values     = {0, 1, 1, 0, 0, 3, -3};
      A                              = loadMatrixFromVectors<Matrix>(3, 3, rowmap, entries, values);
      std::vector<int> rowmapFilt    = {0, 1, 2, 4};
      std::vector<int> entriesFilt   = {1, 1, 1, 2};
      std::vector<double> valuesFilt = {1, 1, 3, -3};
      Afiltered_ref                  = loadMatrixFromVectors<Matrix>(3, 3, rowmapFilt, entriesFilt, valuesFilt);
      break;
    }
    case 6: {
      // First and last rows empty
      std::vector<int> rowmap        = {0, 0, 2, 2};
      std::vector<int> entries       = {0, 1};
      std::vector<double> values     = {0, 3.14};
      A                              = loadMatrixFromVectors<Matrix>(3, 2, rowmap, entries, values);
      std::vector<int> rowmapFilt    = {0, 0, 1, 1};
      std::vector<int> entriesFilt   = {1};
      std::vector<double> valuesFilt = {3.14};
      Afiltered_ref                  = loadMatrixFromVectors<Matrix>(3, 2, rowmapFilt, entriesFilt, valuesFilt);
      break;
    }
    case 7: {
      // First and last rows nonempty, but will be empty after filtering
      std::vector<int> rowmap        = {0, 2, 4, 6};
      std::vector<int> entries       = {0, 1, 1, 2, 0, 3};
      std::vector<double> values     = {0, 0, 1, -1, 0, 0};
      A                              = loadMatrixFromVectors<Matrix>(3, 4, rowmap, entries, values);
      std::vector<int> rowmapFilt    = {0, 0, 2, 2};
      std::vector<int> entriesFilt   = {1, 2};
      std::vector<double> valuesFilt = {1, -1};
      Afiltered_ref                  = loadMatrixFromVectors<Matrix>(3, 4, rowmapFilt, entriesFilt, valuesFilt);
      break;
    }
    case 8: {
      // Large, random matrix with 30% of values converted to zero
      Offset nnz      = 40 * 10000;
      A               = KokkosSparse::Impl::kk_generate_sparse_matrix<Matrix>(10000, 10000, nnz, 10, 5000);
      auto valuesHost = Kokkos::create_mirror_view(A.values);
      Kokkos::deep_copy(valuesHost, A.values);
      for (Offset i = 0; i < A.nnz(); i++) {
        if (rand() % 10 < 3) valuesHost(i) = 0.0;
      }
      Kokkos::deep_copy(A.values, valuesHost);
      Afiltered_ref          = removeMatrixZerosReference(A);
      haveHardcodedReference = false;
      break;
    }
    case 9: {
      // Large, sparser random matrix with 99% of values converted to zero
      Offset nnz      = 10 * 40000;
      A               = KokkosSparse::Impl::kk_generate_sparse_matrix<Matrix>(40000, 40000, nnz, 10, 10000);
      auto valuesHost = Kokkos::create_mirror_view(A.values);
      Kokkos::deep_copy(valuesHost, A.values);
      for (Offset i = 0; i < A.nnz(); i++) {
        if (rand() % 100 != 99) valuesHost(i) = 0.0;
      }
      Kokkos::deep_copy(A.values, valuesHost);
      Afiltered_ref          = removeMatrixZerosReference(A);
      haveHardcodedReference = false;
      break;
    }
    default: throw std::invalid_argument("Test case number of out bounds");
  }
  // If we have a hardcoded reference, check that the reference impl is correct
  // on this case
  if (haveHardcodedReference) {
    Matrix Afiltered_refimpl           = removeMatrixZerosReference(A);
    bool referenceImplMatchesHardcoded = Test::is_same_matrix<Matrix, TestDevice>(Afiltered_ref, Afiltered_refimpl);
    ASSERT_TRUE(referenceImplMatchesHardcoded) << "Test case " << test << ": reference impl gave wrong answer!";
  }
}

}  // namespace TestRemoveCrsMatrixZeros

void testRemoveCrsMatrixZeros(int testCase) {
  using namespace TestRemoveCrsMatrixZeros;
  using Matrix = KokkosSparse::CrsMatrix<default_scalar, default_lno_t, TestDevice, void, default_size_type>;
  Matrix A, Afiltered_ref;
  getTestInput<Matrix>(testCase, A, Afiltered_ref);
  Matrix Afiltered_actual = KokkosSparse::removeCrsMatrixZeros(A);
  bool matches            = Test::is_same_matrix<Matrix, TestDevice>(Afiltered_actual, Afiltered_ref);
  EXPECT_TRUE(matches) << "Test case " << testCase << ": matrix with zeros filtered out does not match reference.";
}

TEST_F(TestCategory, sparse_remove_crs_zeros) {
  for (int testCase = 0; testCase < 10; testCase++) testRemoveCrsMatrixZeros(testCase);
}

#endif  // KOKKOSSPARSE_REMOVECRSZEROS_HPP
