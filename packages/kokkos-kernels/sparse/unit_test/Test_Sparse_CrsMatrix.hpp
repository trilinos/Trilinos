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

// #include "KokkosKernels_ETIHelperMacros.h"
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <stdexcept>
#include "KokkosSparse_CrsMatrix.hpp"
#include "Kokkos_ArithTraits.hpp"

// #ifndef kokkos_complex_double
// #define kokkos_complex_double Kokkos::complex<double>
// #define kokkos_complex_float Kokkos::complex<float>
// #endif

typedef Kokkos::complex<double> kokkos_complex_double;
typedef Kokkos::complex<float> kokkos_complex_float;

namespace Test {  // anonymous

using std::cerr;
using std::endl;

// Create a test sparse matrix A.
//
// Identify the matrix to create by number (whichMatrix).  The
// following lists the valid options for whichMatrix:
//
// 0: A square 10 x 10 nonsymmetric diagonally dominant sparse matrix.
//
// \param ptr [out] Array of row offsets, of length numRows+1.
// \param ind [out] Array of column indices, of length nnz.
// \param val [out] Array of entries (values), of length nnz.
// \param numRows [out] The number of rows in the matrix.
// \param numCols [out] The number of columns in the matrix.
// \param nnz [out] The number of stored entries in the matrix.
// \param whichMatrix [in] The index of the matrix to create.
template <typename crsMat_t>
void makeSparseMatrix(typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type &ptr,
                      typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type &ind,
                      typename crsMat_t::values_type::non_const_type &val, typename crsMat_t::ordinal_type &numRows,
                      typename crsMat_t::ordinal_type &numCols, typename crsMat_t::size_type &nnz,
                      const int whichMatrix) {
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type ptr_type;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type ind_type;
  typedef typename crsMat_t::values_type::non_const_type val_type;
  typedef typename crsMat_t::ordinal_type lno_t;
  typedef typename crsMat_t::size_type size_type;
  typedef typename crsMat_t::value_type scalar_t;

  using Kokkos::HostSpace;
  using Kokkos::MemoryUnmanaged;
  using Kokkos::View;

  if (whichMatrix == 0) {
    numRows                  = 10;
    numCols                  = 10;
    nnz                      = 21;
    const size_type ptrRaw[] = {0, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21};
    const lno_t indRaw[]     = {0, 1, 9, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 1, 9};

    const scalar_t valRaw[] = {1.0, 4.0, 0.5, 0.5,  5.0, 1.0,  6.0, 1.5,  7.0, 2.0, 8.0,
                               2.5, 9.0, 3.0, 10.0, 3.5, 11.0, 4.0, 12.0, 4.5, 13.0};

    // Create the output Views.
    ptr = ptr_type("ptr", numRows + 1);
    ind = ind_type("ind", nnz);
    val = val_type("val", nnz);

    // Wrap the above three arrays in unmanaged Views, so we can use deep_copy.
    typename ptr_type::HostMirror::const_type ptrIn(ptrRaw, numRows + 1);
    typename ind_type::HostMirror::const_type indIn(indRaw, nnz);
    typename val_type::HostMirror::const_type valIn(valRaw, nnz);

    Kokkos::deep_copy(ptr, ptrIn);
    Kokkos::deep_copy(ind, indIn);
    Kokkos::deep_copy(val, valIn);
  } else {  // whichMatrix != 0
    std::ostringstream os;
    os << "Invalid whichMatrix value " << whichMatrix << ".  Valid value(s) include " << 0 << ".";
    throw std::invalid_argument(os.str());
  }
}

// Return the Kokkos::CrsMatrix corresponding to makeSparseMatrix().
template <typename crsMat_t>
crsMat_t makeCrsMatrix() {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::ordinal_type lno_t;
  typedef typename crsMat_t::size_type size_type;

  lno_view_t ptr;
  lno_nnz_view_t ind;
  scalar_view_t val;
  lno_t numRows;
  lno_t numCols;
  size_type nnz;

  const int whichMatrix = 0;
  makeSparseMatrix<crsMat_t>(ptr, ind, val, numRows, numCols, nnz, whichMatrix);
  return crsMat_t("A", numRows, numCols, nnz, val, ptr, ind);
}

}  // namespace Test
// Create a Kokkos::CrsMatrix.  This mainly tests that the class
// compiles.  However, it does need to initialize the MemorySpace's
// default execution space, because it allocates Views and calls
// deep_copy a few times.
template <typename scalar_t, typename lno_t, typename size_type, typename device>
void testCrsMatrix() {
  using namespace Test;

  typedef KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crs_matrix_type;
  crs_matrix_type A = makeCrsMatrix<crs_matrix_type>();
  // mfh 28 Sep 2013: Use A in some way, so the compiler can't
  // optimize it away completely.  This forces the compiler to
  // compile CrsMatrix, which is the whole point of this test.
  // printf ("A is %d by %d\n", A.numRows (), A.numCols ());
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void testCrsMatrixRawConstructor() {
  int nrows = 5;
  // note: last 2 columns will be empty.
  // This makes sure the ncols provided to constructor is preserved.
  int ncols = 7;
  int nnz   = 9;
  // NOTE: this is not a mistake, the raw ptr constructor takes rowmap as
  // ordinal.
  std::vector<lno_t> rowmap  = {0, 0, 2, 5, 6, 9};
  std::vector<lno_t> entries = {3, 4, 0, 1, 2, 2, 0, 3, 4};
  std::vector<scalar_t> values;
  for (int i = 0; i < nnz; i++) values.push_back(Kokkos::ArithTraits<scalar_t>::one() * (1.0 * rand() / RAND_MAX));
  KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> A("A", nrows, ncols, nnz, values.data(),
                                                                      rowmap.data(), entries.data());
  EXPECT_EQ(A.numRows(), nrows);
  EXPECT_EQ(A.numCols(), ncols);
  EXPECT_EQ(A.nnz(), nnz);
  // verify rowmap, entries, values: should all be identical to original raw
  // arrays (except the rowmap elements are now size_type)
  auto checkRowmap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto checkEntries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  auto checkValues  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
  for (int i = 0; i < nrows + 1; i++) EXPECT_EQ(checkRowmap(i), (size_type)rowmap[i]);
  for (int i = 0; i < nnz; i++) {
    EXPECT_EQ(checkEntries(i), entries[i]);
    EXPECT_EQ(checkValues(i), values[i]);
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void testCrsMatrixHostMirror() {
  using namespace Test;
  using crs_matrix      = KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using crs_matrix_host = typename crs_matrix::HostMirror;
  using crs_graph       = typename crs_matrix::StaticCrsGraphType;
  using crs_graph_host  = typename crs_graph::HostMirror;
  crs_matrix A          = makeCrsMatrix<crs_matrix>();
  typename crs_matrix::values_type::HostMirror valuesHost("values host", A.nnz());
  typename crs_matrix::row_map_type::HostMirror rowmapHost("rowmap host", A.numRows() + 1);
  typename crs_matrix::index_type::HostMirror entriesHost("entries host", A.nnz());
  crs_graph_host graphHost(entriesHost, rowmapHost);
  // Test the two CrsMatrix constructors that take the StaticCrsGraph
  crs_matrix_host Ahost1("Ahost1", graphHost, A.numCols());
  crs_matrix_host Ahost2("Ahost2", A.numCols(), valuesHost, graphHost);
  // Test deep copy constructor (can copy between any two spaces)
  {
    crs_matrix Bdev("B device", Ahost1);
    crs_matrix_host Bhost("B host", A);
  }
  // Test the empty (0x0, 0 entries) case - zero-length rowmap.
  typename crs_graph::row_map_type::non_const_type zeroRowmap;
  typename crs_graph::entries_type zeroEntries;
  typename crs_matrix::values_type zeroValues;
  crs_matrix zero("ZeroRow", 0, 0, 0, zeroValues, zeroRowmap, zeroEntries);
  crs_matrix_host zeroHost("zero1Host", zero);
  EXPECT_EQ(zeroHost.numRows(), 0);
  EXPECT_EQ(zeroHost.numCols(), 0);
  EXPECT_EQ(zeroHost.nnz(), 0);
  EXPECT_EQ(zeroHost.graph.row_map.extent(0), 0);
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                                     \
  TEST_F(TestCategory, sparse##_##crsmatrix##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) {             \
    testCrsMatrix<SCALAR, ORDINAL, OFFSET, DEVICE>();                                                   \
    testCrsMatrixRawConstructor<SCALAR, ORDINAL, OFFSET, DEVICE>();                                     \
  }                                                                                                     \
  TEST_F(TestCategory, sparse##_##crsmatrix_host_mirror##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    testCrsMatrixHostMirror<SCALAR, ORDINAL, OFFSET, DEVICE>();                                         \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
