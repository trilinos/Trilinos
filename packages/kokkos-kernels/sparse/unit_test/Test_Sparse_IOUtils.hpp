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
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_Utils.hpp"
#include "Test_vector_fixtures.hpp"

#include <fstream>

namespace Test {

struct TestIOUtils {
  using size_type = size_t;
  using lno_t     = int;
  using scalar_t  = double;

  using exe_space   = Kokkos::DefaultHostExecutionSpace;
  using mem_space   = typename exe_space::memory_space;
  using host_device = Kokkos::Device<exe_space, mem_space>;

  using RowMapType  = Kokkos::View<size_type*, host_device>;
  using EntriesType = Kokkos::View<lno_t*, host_device>;
  using ValuesType  = Kokkos::View<scalar_t*, host_device>;

  using sp_matrix_type = KokkosSparse::CrsMatrix<scalar_t, lno_t, host_device, void, size_type>;

  static std::vector<std::vector<scalar_t>> get_sym_fixture() {
    std::vector<std::vector<scalar_t>> A = {
        {11.00, 12.00, 13.00, 14.00, 15.00, 16.00}, {12.00, 2.00, 0.00, 0.00, 0.00, 0.00},
        {13.00, 0.00, 0.00, 0.00, 0.00, 0.00},      {14.00, 0.00, 0.00, 4.00, 0.00, 0.00},
        {15.00, 0.00, 0.00, 0.00, 5.00, 0.00},      {16.00, 0.00, 0.00, 0.00, 0.00, 6.00}};
    return A;
  }

  static std::vector<std::vector<scalar_t>> get_asym_fixture() {
    std::vector<std::vector<scalar_t>> A = {{1.00, 0.00, 0.00, 9.00, 0.00, 0.00}, {0.00, 2.00, 0.00, 0.00, 0.00, 0.00},
                                            {0.00, 0.00, 0.00, 0.00, 0.00, 8.00}, {0.00, 0.00, 0.00, 4.00, 0.00, 0.00},
                                            {0.00, 7.00, 0.00, 0.00, 5.00, 0.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 6.00}};
    return A;
  }

  static void compare_matrices(const sp_matrix_type& A1, const sp_matrix_type& A2) {
    // Compare matrices
    auto row_map1 = A1.graph.row_map;
    auto entries1 = A1.graph.entries;
    auto values1  = A1.values;
    auto row_map2 = A2.graph.row_map;
    auto entries2 = A2.graph.entries;
    auto values2  = A2.values;
    ASSERT_EQ(row_map1.size(), row_map2.size());
    ASSERT_EQ(entries1.size(), entries2.size());
    ASSERT_EQ(values1.size(), values2.size());
    ASSERT_EQ(values1.size(), entries1.size());
    for (size_type i = 0; i < row_map1.size(); ++i) {
      EXPECT_EQ(row_map1(i), row_map2(i));
    }
    for (size_type i = 0; i < entries1.size(); ++i) {
      EXPECT_EQ(entries1(i), entries2(i));
      EXPECT_EQ(values1(i), values2(i));
    }
  }

  template <typename RowMapView, typename EntriesView, typename ValuesView>
  static void write_as_hb(const RowMapView& row_map, const EntriesView& entries, const ValuesView& values,
                          const std::string& filename, const char mtx_type) {
    std::ofstream out(filename);
    size_type nrows = row_map.size() - 1;
    size_type nnz   = entries.size();

    out << "1SYMMETRIC MATRIX, FE APPROXIMATION TO BIHARMONIC OPERATOR ON "
           "BEAM.     NOS1    \n";  // Title is inaccurate, but doesn't matter
    out << "             3             1             1             1           "
           "  0          \n";
    out << "R" << mtx_type << "A                        " << nrows << "             " << nrows << "             " << nnz
        << "             0          \n";
    out << "(16I5)          (16I5)          (5E16.8)                           "
           "             \n";
    for (size_type row_idx = 0; row_idx < nrows + 1; ++row_idx) {
      out << row_map(row_idx) + 1 << " ";
    }
    out << "\n";
    for (size_type n = 0; n < nnz; ++n) {
      out << entries[n] + 1 << " ";
    }
    out << "\n";
    for (size_type n = 0; n < nnz; ++n) {
      out << values[n] << " ";
    }
    out << "\n";

    out.close();
  }

  template <typename RowMapView, typename EntriesView, typename ValuesView>
  static void write_as_mtx(const RowMapView& row_map, const EntriesView& entries, const ValuesView& values,
                           const std::string& filename, const char mtx_type) {
    std::ofstream out(filename);
    size_type nrows = row_map.size() - 1;

    std::map<char, std::string> type_name_map = {
        {'U', "general"}, {'S', "symmetric"}, {'H', "hermitian"}, {'Z', "skew-symmetric"}};
    std::string type_name = type_name_map[mtx_type];

    out << "%%MatrixMarket matrix coordinate real " << type_name << "\n";
    out << nrows << " " << nrows << " " << entries.size() << "\n";
    for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
      const size_type row_nnz_begin = row_map(row_idx);
      const size_type row_nnz_end   = row_map(row_idx + 1);
      for (size_type row_nnz = row_nnz_begin; row_nnz < row_nnz_end; ++row_nnz) {
        const auto col_idx   = entries(row_nnz);
        const scalar_t value = values(row_nnz);
        out << row_idx + 1 << " " << col_idx + 1 << " " << value << "\n";
      }
    }

    out.close();
  }

  static void full_test(const std::vector<std::vector<scalar_t>>& fixture, const std::string& filename_root,
                        const char mtx_type) {
    RowMapType row_map;
    EntriesType entries;
    ValuesType values;
    compress_matrix(row_map, entries, values, fixture);
    sp_matrix_type A("A", row_map.size() - 1, row_map.size() - 1, values.extent(0), values, row_map, entries);
    const bool is_symmetric = mtx_type != 'U';
    std::string hb_file     = filename_root + ".hb";
    std::string mtx_file    = filename_root + ".mtx";

    if (is_symmetric) {
      sp_matrix_type L = KokkosSparse::Impl::kk_get_lower_triangle(A, NULL, false, 4, true, true);
      auto lrow_map    = L.graph.row_map;
      auto lentries    = L.graph.entries;
      auto lvalues     = L.values;

      write_as_hb(lrow_map, lentries, lvalues, hb_file, mtx_type);
      write_as_mtx(lrow_map, lentries, lvalues, mtx_file, mtx_type);
    } else {
      write_as_hb(row_map, entries, values, hb_file, mtx_type);
      write_as_mtx(row_map, entries, values, mtx_file, mtx_type);
    }

    auto Ahb  = KokkosSparse::Impl::read_kokkos_crst_matrix<sp_matrix_type>(hb_file.c_str());
    auto Amtx = KokkosSparse::Impl::read_kokkos_crst_matrix<sp_matrix_type>(mtx_file.c_str());
    if (mtx_type == 'Z') {
      compare_matrices(Ahb, Amtx);
    } else {
      compare_matrices(Ahb, A);
      compare_matrices(Amtx, A);
    }
  }

  static void test() {
    const std::string filename_root = "test_sparse_ioutils";
    auto sym_fix                    = get_sym_fixture();
    auto asym_fix                   = get_asym_fixture();
    full_test(asym_fix, filename_root + "_asym", 'U');
    full_test(sym_fix, filename_root + "_sym", 'S');
    full_test(sym_fix, filename_root + "_herm", 'H');
    full_test(sym_fix, filename_root + "_skew", 'Z');
  }
};

// Test randomly generated Cs matrices
TEST_F(TestCategory, sparse_ioutils) { TestIOUtils::test(); }

}  // namespace Test
