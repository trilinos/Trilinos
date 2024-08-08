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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <string>
#include <stdexcept>

#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include <KokkosKernels_IOUtils.hpp>
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_spiluk.hpp"
#include "KokkosSparse_crs_to_bsr_impl.hpp"
#include "KokkosSparse_bsr_to_crs_impl.hpp"
#include "KokkosSparse_LUPrec.hpp"
#include "KokkosSparse_gmres.hpp"

#include "Test_vector_fixtures.hpp"

#include <tuple>
#include <random>

using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

using kokkos_complex_double = Kokkos::complex<double>;
using kokkos_complex_float  = Kokkos::complex<float>;

// Comment this out to do focussed debugging
#define TEST_SPILUK_FULL_CHECKS

// Test verbosity level. 0 = none, 1 = print residuals, 2 = print L,U
#define TEST_SPILUK_VERBOSE_LEVEL 0

// #define TEST_SPILUK_TINY_TEST

namespace Test {

#ifdef TEST_SPILUK_TINY_TEST
template <typename scalar_t>
std::vector<std::vector<scalar_t>> get_fixture() {
  std::vector<std::vector<scalar_t>> A = {
      {10.00, 1.00, 0.00, 0.00}, {0.00, 11.00, 0.00, 0.00}, {0.00, 2.00, 12.00, 0.00}, {5.00, 0.00, 3.00, 13.00}};
  return A;
}
#else
template <typename scalar_t>
std::vector<std::vector<scalar_t>> get_fixture() {
  std::vector<std::vector<scalar_t>> A = {
      {10.00, 0.00, 0.30, 0.00, 0.00, 0.60, 0.00, 0.00, 0.00}, {0.00, 11.00, 0.00, 0.00, 0.00, 0.00, 0.70, 0.00, 0.00},
      {0.00, 0.00, 12.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00}, {5.00, 0.00, 0.00, 13.00, 1.00, 0.00, 0.00, 0.00, 0.00},
      {4.00, 0.00, 0.00, 0.00, 14.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 3.00, 0.00, 0.00, 0.00, 15.00, 0.00, 0.00, 0.00},
      {0.00, 0.00, 7.00, 0.00, 0.00, 0.00, 16.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 6.00, 5.00, 0.00, 0.00, 17.00, 0.00},
      {0.00, 0.00, 0.00, 2.00, 2.50, 0.00, 0.00, 0.00, 18.00}};
  return A;
}
#endif

template <typename MatrixType, typename CRS, typename std::enable_if<is_crs_matrix<MatrixType>::value>::type* = nullptr>
MatrixType get_A(CRS A_unblocked, const size_t) {
  return A_unblocked;
}

template <typename MatrixType, typename CRS, typename std::enable_if<is_bsr_matrix<MatrixType>::value>::type* = nullptr>
MatrixType get_A(CRS A_unblocked, const size_t block_size) {
  // Convert to BSR
  MatrixType A(A_unblocked, block_size);

  return A;
}

template <typename MatrixType, typename RowMapType, typename EntriesType, typename ValuesType,
          typename std::enable_if<is_crs_matrix<MatrixType>::value>::type* = nullptr>
MatrixType make_matrix(const char* name, const RowMapType& row_map, const EntriesType& entries,
                       const ValuesType& values, const size_t) {
  const auto nrows = row_map.extent(0) - 1;
  return MatrixType(name, nrows, nrows, values.extent(0), values, row_map, entries);
}

template <typename MatrixType, typename RowMapType, typename EntriesType, typename ValuesType,
          typename std::enable_if<is_bsr_matrix<MatrixType>::value>::type* = nullptr>
MatrixType make_matrix(const char* name, const RowMapType& row_map, const EntriesType& entries,
                       const ValuesType& values, const size_t block_size) {
  const auto nrows = row_map.extent(0) - 1;
  return MatrixType(name, nrows, nrows, values.extent(0), values, row_map, entries, block_size);
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
struct SpilukTest {
  using RowMapType  = Kokkos::View<size_type*, device>;
  using EntriesType = Kokkos::View<lno_t*, device>;
  using ValuesType  = Kokkos::View<scalar_t*, device>;
  using AT          = Kokkos::ArithTraits<scalar_t>;
  using mag_t       = typename Kokkos::ArithTraits<scalar_t>::mag_type;

  using RowMapType_hostmirror  = typename RowMapType::HostMirror;
  using EntriesType_hostmirror = typename EntriesType::HostMirror;
  using ValuesType_hostmirror  = typename ValuesType::HostMirror;
  using execution_space        = typename device::execution_space;
  using memory_space           = typename device::memory_space;
  using range_policy           = Kokkos::RangePolicy<execution_space>;

  static constexpr double EPS = std::is_same<mag_t, double>::value ? 1e-7 : 1e-4;

  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, execution_space,
                                                                        memory_space, memory_space>;

  using Crs = CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using Bsr = BsrMatrix<scalar_t, lno_t, device, void, size_type>;

  template <typename AType, typename LType, typename UType>
  static typename AT::mag_type check_result_impl(const AType& A, const LType& L, const UType& U, const size_type nrows,
                                                 const size_type block_size = 1) {
    const scalar_t ZERO = scalar_t(0);
    const scalar_t ONE  = scalar_t(1);
    const scalar_t MONE = scalar_t(-1);

    // Create a reference view e set to all 1's
    ValuesType e_one("e_one", nrows * block_size);
    Kokkos::deep_copy(e_one, ONE);

    // Create two views for spmv results
    ValuesType bb("bb", nrows * block_size);
    ValuesType bb_tmp("bb_tmp", nrows * block_size);

    // Compute norm2(L*U*e_one - A*e_one)/norm2(A*e_one)
    KokkosSparse::spmv("N", ONE, A, e_one, ZERO, bb);

    typename AT::mag_type bb_nrm = KokkosBlas::nrm2(bb);

    KokkosSparse::spmv("N", ONE, U, e_one, ZERO, bb_tmp);
    KokkosSparse::spmv("N", ONE, L, bb_tmp, MONE, bb);

    typename AT::mag_type diff_nrm = KokkosBlas::nrm2(bb);

    return diff_nrm / bb_nrm;
  }

  static bool is_triangular(const RowMapType& drow_map, const EntriesType& dentries, bool check_lower) {
    const size_type nrows = drow_map.extent(0) - 1;

    auto row_map = Kokkos::create_mirror_view(drow_map);
    auto entries = Kokkos::create_mirror_view(dentries);
    Kokkos::deep_copy(row_map, drow_map);
    Kokkos::deep_copy(entries, dentries);

    for (size_type row = 0; row < nrows; ++row) {
      const size_type row_nnz_begin = row_map(row);
      const size_type row_nnz_end   = row_map(row + 1);
      for (size_type nnz = row_nnz_begin; nnz < row_nnz_end; ++nnz) {
        const size_type col = entries(nnz);
        if (col > row && check_lower) {
          return false;
        } else if (col < row && !check_lower) {
          return false;
        }
      }
    }
    return true;
  }

  template <bool UseBlocks>
  static void check_result(const RowMapType& row_map, const EntriesType& entries, const ValuesType& values,
                           const RowMapType& L_row_map, const EntriesType& L_entries, const ValuesType& L_values,
                           const RowMapType& U_row_map, const EntriesType& U_entries, const ValuesType& U_values,
                           const lno_t fill_lev, const size_type block_size = 1) {
    using sp_matrix_type = std::conditional_t<UseBlocks, Bsr, Crs>;

    KK_REQUIRE(UseBlocks || (block_size == 1));

    // Checking
    const auto nrows = row_map.extent(0) - 1;
    auto A           = make_matrix<sp_matrix_type>("A_Mtx", row_map, entries, values, block_size);
    auto L           = make_matrix<sp_matrix_type>("L_Mtx", L_row_map, L_entries, L_values, block_size);
    auto U           = make_matrix<sp_matrix_type>("U_Mtx", U_row_map, U_entries, U_values, block_size);

    EXPECT_TRUE(is_triangular(L_row_map, L_entries, true));
    EXPECT_TRUE(is_triangular(U_row_map, U_entries, false));

    const auto result = check_result_impl(A, L, U, nrows, block_size);
    if (TEST_SPILUK_VERBOSE_LEVEL > 0) {
      std::cout << "For nrows=" << nrows << ", fill_level=" << fill_lev;
      if (UseBlocks) {
        std::cout << ", block_size=" << block_size;
      } else {
        std::cout << ", unblocked";
      }
      std::cout << " had residual: " << result << std::endl;
    }
    if (TEST_SPILUK_VERBOSE_LEVEL > 1) {
      std::cout << "L result" << std::endl;
      print_matrix(decompress_matrix(L_row_map, L_entries, L_values, block_size));
      std::cout << "U result" << std::endl;
      print_matrix(decompress_matrix(U_row_map, U_entries, U_values, block_size));
    }

    if (fill_lev > 1) {
      EXPECT_LT(result, 1e-4);
    }
  }

  template <bool UseBlocks, int TeamSize = -1>
  static std::tuple<RowMapType, EntriesType, ValuesType, RowMapType, EntriesType, ValuesType> run_and_check_spiluk(
      KernelHandle& kh, const RowMapType& row_map, const EntriesType& entries, const ValuesType& values,
      SPILUKAlgorithm alg, const lno_t fill_lev, const size_type block_size = 1) {
    KK_REQUIRE(UseBlocks || (block_size == 1));

    const size_type block_items = block_size * block_size;
    const size_type nrows       = row_map.extent(0) - 1;
    kh.create_spiluk_handle(alg, nrows, 40 * nrows, 40 * nrows, !UseBlocks ? 0 : block_size);

    auto spiluk_handle = kh.get_spiluk_handle();
    if (TeamSize != -1) {
      spiluk_handle->set_team_size(TeamSize);
    }

    // Allocate L and U as outputs
    RowMapType L_row_map("L_row_map", nrows + 1);
    EntriesType L_entries("L_entries", spiluk_handle->get_nnzL());
    RowMapType U_row_map("U_row_map", nrows + 1);
    EntriesType U_entries("U_entries", spiluk_handle->get_nnzU());

    spiluk_symbolic(&kh, fill_lev, row_map, entries, L_row_map, L_entries, U_row_map, U_entries);

    Kokkos::fence();

    Kokkos::resize(L_entries, spiluk_handle->get_nnzL());
    Kokkos::resize(U_entries, spiluk_handle->get_nnzU());
    ValuesType L_values("L_values", spiluk_handle->get_nnzL() * block_items);
    ValuesType U_values("U_values", spiluk_handle->get_nnzU() * block_items);

    spiluk_numeric(&kh, fill_lev, row_map, entries, values, L_row_map, L_entries, L_values, U_row_map, U_entries,
                   U_values);

    Kokkos::fence();

    check_result<UseBlocks>(row_map, entries, values, L_row_map, L_entries, L_values, U_row_map, U_entries, U_values,
                            fill_lev, block_size);

    kh.destroy_spiluk_handle();

#ifdef TEST_SPILUK_FULL_CHECKS
    // If block_size is 1, results should exactly match unblocked results
    if (block_size == 1 && UseBlocks) {
      const auto [L_row_map_u, L_entries_u, L_values_u, U_row_map_u, U_entries_u, U_values_u] =
          run_and_check_spiluk<false, TeamSize>(kh, row_map, entries, values, alg, fill_lev);

      EXPECT_NEAR_KK_1DVIEW(L_row_map, L_row_map_u, EPS);
      EXPECT_NEAR_KK_1DVIEW(L_entries, L_entries_u, EPS);
      EXPECT_NEAR_KK_1DVIEW(L_values, L_values_u, EPS);
      EXPECT_NEAR_KK_1DVIEW(U_row_map, U_row_map_u, EPS);
      EXPECT_NEAR_KK_1DVIEW(U_entries, U_entries_u, EPS);
      EXPECT_NEAR_KK_1DVIEW(U_values, U_values_u, EPS);
    }

    // Check that team size = 1 produces same result
    if (TeamSize != 1) {
      const auto [L_row_map_ts1, L_entries_ts1, L_values_ts1, U_row_map_ts1, U_entries_ts1, U_values_ts1] =
          run_and_check_spiluk<UseBlocks, 1>(kh, row_map, entries, values, alg, fill_lev, block_size);

      EXPECT_NEAR_KK_1DVIEW(L_row_map, L_row_map_ts1, EPS);
      EXPECT_NEAR_KK_1DVIEW(L_entries, L_entries_ts1, EPS);
      EXPECT_NEAR_KK_1DVIEW(L_values, L_values_ts1, EPS);
      EXPECT_NEAR_KK_1DVIEW(U_row_map, U_row_map_ts1, EPS);
      EXPECT_NEAR_KK_1DVIEW(U_entries, U_entries_ts1, EPS);
      EXPECT_NEAR_KK_1DVIEW(U_values, U_values_ts1, EPS);
    }
#endif

    return std::make_tuple(L_row_map, L_entries, L_values, U_row_map, U_entries, U_values);
  }

  static void run_test_spiluk() {
    std::vector<std::vector<scalar_t>> A = get_fixture<scalar_t>();

    if (TEST_SPILUK_VERBOSE_LEVEL > 1) {
      std::cout << "A input" << std::endl;
      print_matrix(A);
    }

    RowMapType row_map;
    EntriesType entries;
    ValuesType values;

    compress_matrix(row_map, entries, values, A);

    const lno_t fill_lev = 2;

    KernelHandle kh;

    run_and_check_spiluk<false>(kh, row_map, entries, values, SPILUKAlgorithm::SEQLVLSCHD_TP1, fill_lev);
  }

  static void run_test_spiluk_blocks() {
    std::vector<std::vector<scalar_t>> A = get_fixture<scalar_t>();

    if (TEST_SPILUK_VERBOSE_LEVEL > 1) {
      std::cout << "A input" << std::endl;
      print_matrix(A);
    }

    RowMapType row_map, brow_map;
    EntriesType entries, bentries;
    ValuesType values, bvalues;

    compress_matrix(row_map, entries, values, A);

    const size_type nrows      = A.size();
    const size_type nnz        = values.extent(0);
    const lno_t fill_lev       = 2;
    const size_type block_size = nrows % 2 == 0 ? 2 : 3;
    ASSERT_EQ(nrows % block_size, 0);

    KernelHandle kh;

    Crs crs("crs for block spiluk test", nrows, nrows, nnz, values, row_map, entries);

    std::vector<size_type> block_sizes = {1, block_size};

    for (auto block_size_itr : block_sizes) {
      Bsr bsr(crs, block_size_itr);

      // Pull out views from BSR
      Kokkos::resize(brow_map, bsr.graph.row_map.extent(0));
      Kokkos::resize(bentries, bsr.graph.entries.extent(0));
      Kokkos::resize(bvalues, bsr.values.extent(0));
      Kokkos::deep_copy(brow_map, bsr.graph.row_map);
      Kokkos::deep_copy(bentries, bsr.graph.entries);
      Kokkos::deep_copy(bvalues, bsr.values);

      run_and_check_spiluk<true>(kh, brow_map, bentries, bvalues, SPILUKAlgorithm::SEQLVLSCHD_TP1, fill_lev,
                                 block_size_itr);
    }
  }

  static void run_test_spiluk_scale() {
    // Create a diagonally dominant sparse matrix to test:
    constexpr auto nrows         = 5000;
    constexpr auto diagDominance = 2;

    size_type nnz = 10 * nrows;
    auto A        = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<Crs>(nrows, nrows, nnz, 0,
                                                                                    lno_t(0.01 * nrows), diagDominance);

    KokkosSparse::sort_crs_matrix(A);

    // Pull out views from CRS
    RowMapType row_map("row_map", A.graph.row_map.extent(0));
    EntriesType entries("entries", A.graph.entries.extent(0));
    ValuesType values("values", A.values.extent(0));
    Kokkos::deep_copy(row_map, A.graph.row_map);
    Kokkos::deep_copy(entries, A.graph.entries);
    Kokkos::deep_copy(values, A.values);

    for (lno_t fill_lev = 0; fill_lev < 4; ++fill_lev) {
      KernelHandle kh;

      run_and_check_spiluk<false>(kh, row_map, entries, values, SPILUKAlgorithm::SEQLVLSCHD_TP1, fill_lev);
    }
  }

  static void run_test_spiluk_scale_blocks() {
    // Create a diagonally dominant sparse matrix to test:
    constexpr auto nrows         = 5000;
    constexpr auto diagDominance = 2;

    RowMapType brow_map;
    EntriesType bentries;
    ValuesType bvalues;

    // const size_type block_size = 10;

    size_type nnz = 10 * nrows;
    auto A        = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<Crs>(nrows, nrows, nnz, 0,
                                                                                    lno_t(0.01 * nrows), diagDominance);

    KokkosSparse::sort_crs_matrix(A);

    std::vector<size_type> block_sizes = {1, 2, 4, 10};

    for (auto block_size : block_sizes) {
      // Convert to BSR
      Bsr bsr(A, block_size);

      // Pull out views from BSR
      Kokkos::resize(brow_map, bsr.graph.row_map.extent(0));
      Kokkos::resize(bentries, bsr.graph.entries.extent(0));
      Kokkos::resize(bvalues, bsr.values.extent(0));
      Kokkos::deep_copy(brow_map, bsr.graph.row_map);
      Kokkos::deep_copy(bentries, bsr.graph.entries);
      Kokkos::deep_copy(bvalues, bsr.values);

      for (lno_t fill_lev = 0; fill_lev < 4; ++fill_lev) {
        KernelHandle kh;

        run_and_check_spiluk<true>(kh, brow_map, bentries, bvalues, SPILUKAlgorithm::SEQLVLSCHD_TP1, fill_lev,
                                   block_size);
      }
    }
  }

  static void run_test_spiluk_streams(SPILUKAlgorithm test_algo, int nstreams) {
    // Workaround for OpenMP: skip tests if concurrency < nstreams because of
    // not enough resource to partition
    bool run_streams_test = true;
#ifdef KOKKOS_ENABLE_OPENMP
    if (std::is_same<typename device::execution_space, Kokkos::OpenMP>::value) {
      int exec_concurrency = execution_space().concurrency();
      if (exec_concurrency < nstreams) {
        run_streams_test = false;
        std::cout << "  Skip stream test: concurrency = " << exec_concurrency << std::endl;
      }
    }
#endif
    if (!run_streams_test) return;

    std::vector<int> weights(nstreams, 1);
    std::vector<execution_space> instances = Kokkos::Experimental::partition_space(execution_space(), weights);

    std::vector<KernelHandle> kh_v(nstreams);
    std::vector<KernelHandle*> kh_ptr_v(nstreams);
    std::vector<RowMapType> A_row_map_v(nstreams);
    std::vector<EntriesType> A_entries_v(nstreams);
    std::vector<ValuesType> A_values_v(nstreams);
    std::vector<RowMapType> L_row_map_v(nstreams);
    std::vector<EntriesType> L_entries_v(nstreams);
    std::vector<ValuesType> L_values_v(nstreams);
    std::vector<RowMapType> U_row_map_v(nstreams);
    std::vector<EntriesType> U_entries_v(nstreams);
    std::vector<ValuesType> U_values_v(nstreams);

    std::vector<std::vector<scalar_t>> Afix = get_fixture<scalar_t>();

    RowMapType row_map;
    EntriesType entries;
    ValuesType values;

    compress_matrix(row_map, entries, values, Afix);

    const size_type nrows = Afix.size();
    const size_type nnz   = values.extent(0);

    RowMapType_hostmirror hrow_map("hrow_map", nrows + 1);
    EntriesType_hostmirror hentries("hentries", nnz);
    ValuesType_hostmirror hvalues("hvalues", nnz);

    Kokkos::deep_copy(hrow_map, row_map);
    Kokkos::deep_copy(hentries, entries);
    Kokkos::deep_copy(hvalues, values);

    typename KernelHandle::const_nnz_lno_t fill_lev = 2;

    for (int i = 0; i < nstreams; i++) {
      // Allocate A as input
      A_row_map_v[i] = RowMapType("A_row_map", nrows + 1);
      A_entries_v[i] = EntriesType("A_entries", nnz);
      A_values_v[i]  = ValuesType("A_values", nnz);

      // Copy from host to device
      Kokkos::deep_copy(A_row_map_v[i], hrow_map);
      Kokkos::deep_copy(A_entries_v[i], hentries);
      Kokkos::deep_copy(A_values_v[i], hvalues);

      // Create handle
      kh_v[i] = KernelHandle();
      kh_v[i].create_spiluk_handle(test_algo, nrows, 4 * nrows, 4 * nrows);
      kh_ptr_v[i] = &kh_v[i];

      auto spiluk_handle = kh_v[i].get_spiluk_handle();

      // Allocate L and U as outputs
      L_row_map_v[i] = RowMapType("L_row_map", nrows + 1);
      L_entries_v[i] = EntriesType("L_entries", spiluk_handle->get_nnzL());
      U_row_map_v[i] = RowMapType("U_row_map", nrows + 1);
      U_entries_v[i] = EntriesType("U_entries", spiluk_handle->get_nnzU());

      // Symbolic phase
      spiluk_symbolic(kh_ptr_v[i], fill_lev, A_row_map_v[i], A_entries_v[i], L_row_map_v[i], L_entries_v[i],
                      U_row_map_v[i], U_entries_v[i], nstreams);

      Kokkos::fence();

      Kokkos::resize(L_entries_v[i], spiluk_handle->get_nnzL());
      Kokkos::resize(U_entries_v[i], spiluk_handle->get_nnzU());
      L_values_v[i] = ValuesType("L_values", spiluk_handle->get_nnzL());
      U_values_v[i] = ValuesType("U_values", spiluk_handle->get_nnzU());
    }  // Done handle creation and spiluk_symbolic on all streams

    // Numeric phase
    spiluk_numeric_streams(instances, kh_ptr_v, fill_lev, A_row_map_v, A_entries_v, A_values_v, L_row_map_v,
                           L_entries_v, L_values_v, U_row_map_v, U_entries_v, U_values_v);

    for (int i = 0; i < nstreams; i++) instances[i].fence();

    // Checking
    for (int i = 0; i < nstreams; i++) {
      check_result<false>(A_row_map_v[i], A_entries_v[i], A_values_v[i], L_row_map_v[i], L_entries_v[i], L_values_v[i],
                          U_row_map_v[i], U_entries_v[i], U_values_v[i], fill_lev);

      kh_v[i].destroy_spiluk_handle();
    }
  }

  static void run_test_spiluk_streams_blocks(SPILUKAlgorithm test_algo, int nstreams) {
    // Workaround for OpenMP: skip tests if concurrency < nstreams because of
    // not enough resource to partition
    bool run_streams_test = true;
#ifdef KOKKOS_ENABLE_OPENMP
    if (std::is_same<typename device::execution_space, Kokkos::OpenMP>::value) {
      int exec_concurrency = execution_space().concurrency();
      if (exec_concurrency < nstreams) {
        run_streams_test = false;
        std::cout << "  Skip stream test: concurrency = " << exec_concurrency << std::endl;
      }
    }
#endif
    if (!run_streams_test) return;

    std::vector<int> weights(nstreams, 1);
    std::vector<execution_space> instances = Kokkos::Experimental::partition_space(execution_space(), weights);

    std::vector<KernelHandle> kh_v(nstreams);
    std::vector<KernelHandle*> kh_ptr_v(nstreams);
    std::vector<RowMapType> A_row_map_v(nstreams);
    std::vector<EntriesType> A_entries_v(nstreams);
    std::vector<ValuesType> A_values_v(nstreams);
    std::vector<RowMapType> L_row_map_v(nstreams);
    std::vector<EntriesType> L_entries_v(nstreams);
    std::vector<ValuesType> L_values_v(nstreams);
    std::vector<RowMapType> U_row_map_v(nstreams);
    std::vector<EntriesType> U_entries_v(nstreams);
    std::vector<ValuesType> U_values_v(nstreams);

    std::vector<std::vector<scalar_t>> Afix = get_fixture<scalar_t>();

    RowMapType row_map, brow_map;
    EntriesType entries, bentries;
    ValuesType values, bvalues;

    compress_matrix(row_map, entries, values, Afix);

    const size_type nrows       = Afix.size();
    const size_type block_size  = nrows % 2 == 0 ? 2 : 3;
    const size_type block_items = block_size * block_size;
    ASSERT_EQ(nrows % block_size, 0);

    // Convert to BSR
    Crs crs("crs for block spiluk test", nrows, nrows, values.extent(0), values, row_map, entries);
    Bsr bsr(crs, block_size);

    // Pull out views from BSR
    Kokkos::resize(brow_map, bsr.graph.row_map.extent(0));
    Kokkos::resize(bentries, bsr.graph.entries.extent(0));
    Kokkos::resize(bvalues, bsr.values.extent(0));
    Kokkos::deep_copy(brow_map, bsr.graph.row_map);
    Kokkos::deep_copy(bentries, bsr.graph.entries);
    Kokkos::deep_copy(bvalues, bsr.values);

    const size_type bnrows = brow_map.extent(0) - 1;
    const size_type bnnz   = bentries.extent(0);

    RowMapType_hostmirror hrow_map("hrow_map", bnrows + 1);
    EntriesType_hostmirror hentries("hentries", bnnz);
    ValuesType_hostmirror hvalues("hvalues", bnnz * block_items);

    Kokkos::deep_copy(hrow_map, brow_map);
    Kokkos::deep_copy(hentries, bentries);
    Kokkos::deep_copy(hvalues, bvalues);

    typename KernelHandle::const_nnz_lno_t fill_lev = 2;

    for (int i = 0; i < nstreams; i++) {
      // Allocate A as input
      A_row_map_v[i] = RowMapType("A_row_map", bnrows + 1);
      A_entries_v[i] = EntriesType("A_entries", bnnz);
      A_values_v[i]  = ValuesType("A_values", bnnz * block_items);

      // Copy from host to device
      Kokkos::deep_copy(A_row_map_v[i], hrow_map);
      Kokkos::deep_copy(A_entries_v[i], hentries);
      Kokkos::deep_copy(A_values_v[i], hvalues);

      // Create handle
      kh_v[i] = KernelHandle();
      kh_v[i].create_spiluk_handle(test_algo, bnrows, 4 * bnrows, 4 * bnrows, block_size);
      kh_ptr_v[i] = &kh_v[i];

      auto spiluk_handle = kh_v[i].get_spiluk_handle();

      // Allocate L and U as outputs
      L_row_map_v[i] = RowMapType("L_row_map", bnrows + 1);
      L_entries_v[i] = EntriesType("L_entries", spiluk_handle->get_nnzL());
      U_row_map_v[i] = RowMapType("U_row_map", bnrows + 1);
      U_entries_v[i] = EntriesType("U_entries", spiluk_handle->get_nnzU());

      // Symbolic phase
      spiluk_symbolic(kh_ptr_v[i], fill_lev, A_row_map_v[i], A_entries_v[i], L_row_map_v[i], L_entries_v[i],
                      U_row_map_v[i], U_entries_v[i], nstreams);

      Kokkos::fence();

      Kokkos::resize(L_entries_v[i], spiluk_handle->get_nnzL());
      Kokkos::resize(U_entries_v[i], spiluk_handle->get_nnzU());
      L_values_v[i] = ValuesType("L_values", spiluk_handle->get_nnzL() * block_items);
      U_values_v[i] = ValuesType("U_values", spiluk_handle->get_nnzU() * block_items);
    }  // Done handle creation and spiluk_symbolic on all streams

    // Numeric phase
    spiluk_numeric_streams(instances, kh_ptr_v, fill_lev, A_row_map_v, A_entries_v, A_values_v, L_row_map_v,
                           L_entries_v, L_values_v, U_row_map_v, U_entries_v, U_values_v);

    for (int i = 0; i < nstreams; i++) instances[i].fence();

    // Checking
    for (int i = 0; i < nstreams; i++) {
      check_result<true>(A_row_map_v[i], A_entries_v[i], A_values_v[i], L_row_map_v[i], L_entries_v[i], L_values_v[i],
                         U_row_map_v[i], U_entries_v[i], U_values_v[i], fill_lev, block_size);

      kh_v[i].destroy_spiluk_handle();
    }
  }

  template <bool UseBlocks>
  static void run_test_spiluk_precond() {
    // Test using spiluk as a preconditioner
    // Does (LU)^inv Ax = (LU)^inv b converge faster than solving Ax=b?

    // Create a diagonally dominant sparse matrix to test:
    using sp_matrix_type = std::conditional_t<UseBlocks, Bsr, Crs>;

    constexpr auto nrows         = 5000;
    constexpr auto m             = 15;
    constexpr auto diagDominance = 2;
    constexpr auto tol           = 1e-5;
    constexpr bool verbose       = false;

    if (UseBlocks) {
      // Skip test if not on host. block trsv only works on host
      static constexpr bool is_host = std::is_same<execution_space, typename Kokkos::DefaultHostExecutionSpace>::value;
      if (!is_host) {
        return;
      }
    }

    RowMapType brow_map;
    EntriesType bentries;
    ValuesType bvalues;

    size_type nnz    = 10 * nrows;
    auto A_unblocked = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<Crs>(
        nrows, nrows, nnz, 0, lno_t(0.01 * nrows), diagDominance);

    KokkosSparse::sort_crs_matrix(A_unblocked);

    std::vector<size_type> block_sizes_blocked   = {1, 2, 4, 10};
    std::vector<size_type> block_sizes_unblocked = {1};
    std::vector<size_type> block_sizes           = UseBlocks ? block_sizes_blocked : block_sizes_unblocked;

    for (auto block_size : block_sizes) {
      // Convert to BSR if block enabled
      auto A = get_A<sp_matrix_type>(A_unblocked, block_size);

      // Pull out views from BSR
      Kokkos::resize(brow_map, A.graph.row_map.extent(0));
      Kokkos::resize(bentries, A.graph.entries.extent(0));
      Kokkos::resize(bvalues, A.values.extent(0));
      Kokkos::deep_copy(brow_map, A.graph.row_map);
      Kokkos::deep_copy(bentries, A.graph.entries);
      Kokkos::deep_copy(bvalues, A.values);

      // Make kernel handles
      KernelHandle kh;
      kh.create_gmres_handle(m, tol);
      auto gmres_handle = kh.get_gmres_handle();
      gmres_handle->set_verbose(verbose);
      using GMRESHandle = typename std::remove_reference<decltype(*gmres_handle)>::type;

      for (lno_t fill_lev = 0; fill_lev < 4; ++fill_lev) {
        const auto [L_row_map, L_entries, L_values, U_row_map, U_entries, U_values] = run_and_check_spiluk<UseBlocks>(
            kh, brow_map, bentries, bvalues, SPILUKAlgorithm::SEQLVLSCHD_TP1, fill_lev, block_size);

        // Create L, U
        auto L = make_matrix<sp_matrix_type>("L_Mtx", L_row_map, L_entries, L_values, block_size);
        auto U = make_matrix<sp_matrix_type>("U_Mtx", U_row_map, U_entries, U_values, block_size);

        // Set initial vectors:
        ValuesType X("X", nrows);    // Solution and initial guess
        ValuesType Wj("Wj", nrows);  // For checking residuals at end.
        ValuesType B(Kokkos::view_alloc(Kokkos::WithoutInitializing, "B"),
                     nrows);  // right-hand side vec
        // Make rhs ones so that results are repeatable:
        Kokkos::deep_copy(B, 1.0);

        int num_iters_plain(0), num_iters_precond(0);

        // Solve Ax = b
        {
          gmres(&kh, A, B, X);

          // Double check residuals at end of solve:
          float_t nrmB = KokkosBlas::nrm2(B);
          KokkosSparse::spmv("N", 1.0, A, X, 0.0, Wj);  // wj = Ax
          KokkosBlas::axpy(-1.0, Wj, B);                // b = b-Ax.
          float_t endRes = KokkosBlas::nrm2(B) / nrmB;

          const auto conv_flag = gmres_handle->get_conv_flag_val();
          num_iters_plain      = gmres_handle->get_num_iters();

          EXPECT_GT(num_iters_plain, 0);
          EXPECT_LT(endRes, gmres_handle->get_tol());
          EXPECT_EQ(conv_flag, GMRESHandle::Flag::Conv);

          if (TEST_SPILUK_VERBOSE_LEVEL > 0) {
            std::cout << "Without LUPrec, with block_size=" << block_size << ", converged in " << num_iters_plain
                      << " steps with endres=" << endRes << std::endl;
          }
        }

        // Solve Ax = b with LU preconditioner.
        {
          gmres_handle->reset_handle(m, tol);
          gmres_handle->set_verbose(verbose);

          // Make precond.
          KokkosSparse::Experimental::LUPrec<sp_matrix_type, KernelHandle> myPrec(L, U, UseBlocks ? block_size : 0);

          // reset X for next gmres call
          Kokkos::deep_copy(X, 0.0);

          gmres(&kh, A, B, X, &myPrec);

          // Double check residuals at end of solve:
          float_t nrmB = KokkosBlas::nrm2(B);
          KokkosSparse::spmv("N", 1.0, A, X, 0.0, Wj);  // wj = Ax
          KokkosBlas::axpy(-1.0, Wj, B);                // b = b-Ax.
          float_t endRes = KokkosBlas::nrm2(B) / nrmB;

          const auto conv_flag = gmres_handle->get_conv_flag_val();
          num_iters_precond    = gmres_handle->get_num_iters();

          EXPECT_LT(endRes, gmres_handle->get_tol());
          EXPECT_EQ(conv_flag, GMRESHandle::Flag::Conv);
          EXPECT_LT(num_iters_precond, num_iters_plain);

          if (TEST_SPILUK_VERBOSE_LEVEL > 0) {
            std::cout << "With LUPrec, with block_size=" << block_size << ", and fill_level=" << fill_lev
                      << ", converged in " << num_iters_precond << " steps with endres=" << endRes << std::endl;
          }
        }
      }
    }
  }
};

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_spiluk() {
  using TestStruct = Test::SpilukTest<scalar_t, lno_t, size_type, device>;
  TestStruct::run_test_spiluk();
  TestStruct::run_test_spiluk_blocks();
  TestStruct::run_test_spiluk_scale();
  TestStruct::run_test_spiluk_scale_blocks();
  TestStruct::template run_test_spiluk_precond<false>();
  TestStruct::template run_test_spiluk_precond<true>();
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_spiluk_streams() {
  using TestStruct = Test::SpilukTest<scalar_t, lno_t, size_type, device>;

  TestStruct::run_test_spiluk_streams(SPILUKAlgorithm::SEQLVLSCHD_TP1, 1);
  TestStruct::run_test_spiluk_streams(SPILUKAlgorithm::SEQLVLSCHD_TP1, 2);
  TestStruct::run_test_spiluk_streams(SPILUKAlgorithm::SEQLVLSCHD_TP1, 3);
  TestStruct::run_test_spiluk_streams(SPILUKAlgorithm::SEQLVLSCHD_TP1, 4);

  TestStruct::run_test_spiluk_streams_blocks(SPILUKAlgorithm::SEQLVLSCHD_TP1, 1);
  TestStruct::run_test_spiluk_streams_blocks(SPILUKAlgorithm::SEQLVLSCHD_TP1, 2);
  TestStruct::run_test_spiluk_streams_blocks(SPILUKAlgorithm::SEQLVLSCHD_TP1, 3);
  TestStruct::run_test_spiluk_streams_blocks(SPILUKAlgorithm::SEQLVLSCHD_TP1, 4);
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                      \
  TEST_F(TestCategory, sparse##_##spiluk##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_spiluk<SCALAR, ORDINAL, OFFSET, DEVICE>();                                      \
    test_spiluk_streams<SCALAR, ORDINAL, OFFSET, DEVICE>();                              \
  }

#define NO_TEST_COMPLEX

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
#undef NO_TEST_COMPLEX
