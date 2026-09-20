// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef TEST_SPARSE_UTILS_HPP
#define TEST_SPARSE_UTILS_HPP

#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_SortCrs.hpp"

namespace TestUtils {

template <typename crsMat_t, typename vector_t>
vector_t create_random_y_vector(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Y VECTOR"), crsMat.numRows());
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}

template <typename crsMat_t, typename vector_t>
vector_t create_random_y_vector_mv(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Y VECTOR"), crsMat.numRows(), x_vector.extent(1));
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}

// Test whether two CRS matrices are identical (shape, #nonzeros, entries, and values to within tolerance),
// including the order of entries within each row.
template <typename crsMat_t, typename device>
bool is_same_matrix(crsMat_t output_mat_actual, crsMat_t output_mat_reference) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  size_t nrows_actual    = output_mat_actual.numRows();
  size_t ncols_actual    = output_mat_actual.numCols();
  size_t nentries_actual = output_mat_actual.graph.entries.extent(0);
  size_t nvals_actual    = output_mat_actual.values.extent(0);

  size_t nrows_reference    = output_mat_reference.numRows();
  size_t ncols_reference    = output_mat_reference.numCols();
  size_t nentries_reference = output_mat_reference.graph.entries.extent(0);
  size_t nvals_reference    = output_mat_reference.values.extent(0);

  if (nrows_actual != nrows_reference || ncols_actual != ncols_reference) {
    std::cout << "dimensions (actual):" << nrows_actual << 'x' << ncols_actual
              << ", dimensions (reference): " << nrows_reference << 'x' << ncols_reference << '\n';
    return false;
  }
  if (nentries_actual != nentries_reference) {
    std::cout << "nentries_actual:" << nentries_actual << " nentries_reference:" << nentries_reference << std::endl;
    return false;
  }
  if (nvals_actual != nvals_reference) {
    std::cout << "nvals_actual:" << nvals_actual << " nvals_reference:" << nvals_reference << std::endl;
    return false;
  }

  bool is_identical = true;
  // Special case: a matrix with 0 rows can have a rowmap of length 0 or 1.
  // Treat these as equivalent.
  bool zero_row_equivalent = false;
  if (nrows_reference == 0) {
    auto rm1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output_mat_actual.graph.row_map);
    auto rm2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output_mat_reference.graph.row_map);
    if (rm1.extent_int(0) == 0 && rm2.extent_int(0) == 1) {
      // Make sure the one element of rm2 is 0
      zero_row_equivalent = !rm2(0);
    } else if (rm1.extent_int(0) == 1 && rm2.extent_int(0) == 0) {
      // Make sure the one element of rm1 is 0
      zero_row_equivalent = !rm1(0);
    }
  }
  if (!zero_row_equivalent) {
    is_identical =
        KokkosKernels::Impl::kk_is_identical_view<typename graph_t::row_map_type, typename graph_t::row_map_type,
                                                  typename lno_view_t::value_type, typename device::execution_space>(
            output_mat_actual.graph.row_map, output_mat_reference.graph.row_map, 0);
  }

  if (!is_identical) {
    std::cout << "rowmaps are different." << std::endl;
    std::cout << "Actual rowmap:\n";
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.graph.row_map, true);
    std::cout << "Correct rowmap:\n";
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.graph.row_map, true);
    return false;
  }

  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
                                                typename device::execution_space>(
          output_mat_actual.graph.entries, output_mat_reference.graph.entries, 0);

  if (!is_identical) {
    std::cout << "entries are different." << std::endl;
    auto rowmap_actual = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output_mat_actual.graph.row_map);
    auto rowmap_reference =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output_mat_reference.graph.row_map);
    auto entries_actual = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output_mat_actual.graph.entries);
    auto entries_reference =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output_mat_reference.graph.entries);

    size_t first_mismatch = nentries_actual;
    for (size_t i = 0; i < nentries_actual; ++i) {
      if (entries_actual(i) != entries_reference(i)) {
        first_mismatch = i;
        break;
      }
    }

    if (first_mismatch < nentries_actual) {
      size_t mismatch_row = 0;
      for (size_t r = 0; r + 1 < rowmap_actual.extent(0); ++r) {
        const size_t row_begin = static_cast<size_t>(rowmap_actual(r));
        const size_t row_end   = static_cast<size_t>(rowmap_actual(r + 1));
        if (row_begin <= first_mismatch && first_mismatch < row_end) {
          mismatch_row = r;
          break;
        }
      }

      size_t row_begin_actual = rowmap_actual(mismatch_row);
      size_t row_end_actual   = rowmap_actual(mismatch_row + 1);
      size_t row_begin_ref    = rowmap_reference(mismatch_row);
      size_t row_end_ref      = rowmap_reference(mismatch_row + 1);

      std::cout << "  first mismatch index: " << first_mismatch << '\n';
      std::cout << "  mismatch row: " << mismatch_row << '\n';
      std::cout << "  actual row range: [" << row_begin_actual << ", " << row_end_actual << ")\n";
      std::cout << "  reference row range: [" << row_begin_ref << ", " << row_end_ref << ")\n";

      std::cout << "  actual row entries: ";
      for (size_t i = row_begin_actual; i < row_end_actual; ++i) {
        std::cout << entries_actual(i) << (i + 1 < row_end_actual ? ' ' : '\n');
      }

      std::cout << "  reference row entries: ";
      for (size_t i = row_begin_ref; i < row_end_ref; ++i) {
        std::cout << entries_reference(i) << (i + 1 < row_end_ref ? ' ' : '\n');
      }
    }
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.graph.entries);
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.graph.entries);

    // Check whether the two matrices are at least functionally equivalent
    // under SpMV; this helps distinguish ordering-only mismatches from
    // correctness bugs.
    using scalar_t = typename scalar_view_t::value_type;
    using x_view_t = Kokkos::View<scalar_t *, typename crsMat_t::device_type>;
    x_view_t x("debug_x", output_mat_actual.numCols());
    scalar_t randStart, randEnd;
    KokkosKernels::Impl::getRandomBounds(10.0, randStart, randEnd);
    Kokkos::Random_XorShift64_Pool<typename crsMat_t::execution_space> pool(424242);
    Kokkos::fill_random(x, pool, randStart, randEnd);
    auto y_actual  = create_random_y_vector(output_mat_actual, x);
    auto y_ref     = create_random_y_vector(output_mat_reference, x);
    using mag_t    = typename KokkosKernels::ArithTraits<scalar_t>::mag_type;
    mag_t spmv_eps = std::is_same<mag_t, float>::value ? 3.7e-3 : 1e-7;
    bool spmv_same = KokkosKernels::Impl::kk_is_relatively_identical_view<decltype(y_actual), decltype(y_ref), mag_t,
                                                                          typename device::execution_space>(
        y_actual, y_ref, spmv_eps);
    std::cout << "  SpMV functional equivalence: " << (spmv_same ? "PASS" : "FAIL") << std::endl;
    return false;
  }

  typedef typename KokkosKernels::ArithTraits<typename scalar_view_t::non_const_value_type>::mag_type eps_type;
  eps_type eps = std::is_same<eps_type, float>::value ? 3.7e-3 : 1e-7;

  is_identical = KokkosKernels::Impl::kk_is_relatively_identical_view<scalar_view_t, scalar_view_t, eps_type,
                                                                      typename device::execution_space>(
      output_mat_actual.values, output_mat_reference.values, eps);

  if (!is_identical) {
    std::cout << "values are different." << std::endl;
    KokkosKernels::Impl::kk_print_1Dview(output_mat_actual.values);
    KokkosKernels::Impl::kk_print_1Dview(output_mat_reference.values);

    return false;
  }
  return true;
}

}  // namespace TestUtils

#endif  // TEST_SPARSE_UTILS_HPP
