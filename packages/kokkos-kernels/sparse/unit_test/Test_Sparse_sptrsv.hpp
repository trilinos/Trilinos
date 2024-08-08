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

#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#include "KokkosSparse_sptrsv.hpp"
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
#include "KokkosSparse_sptrsv_supernode.hpp"
#endif

#include <gtest/gtest.h>

using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Impl;
using namespace KokkosKernels::Experimental;

using kokkos_complex_double = Kokkos::complex<double>;
using kokkos_complex_float  = Kokkos::complex<float>;

namespace Test {

template <typename scalar_t, typename lno_t, typename size_type, typename device>
struct SptrsvTest {
  // Define useful types
  using RowMapType      = Kokkos::View<size_type *, device>;
  using EntriesType     = Kokkos::View<lno_t *, device>;
  using ValuesType      = Kokkos::View<scalar_t *, device>;
  using execution_space = typename device::execution_space;
  using memory_space    = typename device::memory_space;
  using KernelHandle    = KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, execution_space,
                                                                        memory_space, memory_space>;

  using Crs = CrsMatrix<scalar_t, lno_t, device, void, size_type>;
  using Bsr = BsrMatrix<scalar_t, lno_t, device, void, size_type>;

  using crs_graph_t = typename Crs::StaticCrsGraphType;

  using range_policy_t = Kokkos::RangePolicy<execution_space>;

  static std::vector<std::vector<scalar_t>> get_5x5_ut_ones_fixture() {
    std::vector<std::vector<scalar_t>> A = {{1.00, 0.00, 1.00, 0.00, 0.00},
                                            {0.00, 1.00, 0.00, 0.00, 1.00},
                                            {0.00, 0.00, 1.00, 1.00, 1.00},
                                            {0.00, 0.00, 0.00, 1.00, 1.00},
                                            {0.00, 0.00, 0.00, 0.00, 1.00}};
    return A;
  }

  static std::vector<std::vector<scalar_t>> get_6x6_ut_ones_fixture() {
    std::vector<std::vector<scalar_t>> A = {{1.00, 1.00, 0.00, 0.00, 0.00, 0.00}, {0.00, 1.00, 0.00, 0.00, 0.00, 1.00},
                                            {0.00, 0.00, 1.00, 1.00, 0.00, 1.00}, {0.00, 0.00, 0.00, 1.00, 0.00, 1.00},
                                            {0.00, 0.00, 0.00, 0.00, 1.00, 1.00}, {0.00, 0.00, 0.00, 0.00, 0.00, 1.00}};
    return A;
  }

  static std::vector<std::vector<scalar_t>> get_5x5_ut_fixture() {
    const auto KZ                        = KEEP_ZERO<scalar_t>();
    std::vector<std::vector<scalar_t>> A = {{5.00, 1.00, 1.00, 0.00, KZ},
                                            {KZ, 5.00, KZ, 0.00, 1.00},
                                            {0.00, 0.00, 5.00, 1.00, 1.00},
                                            {0.00, 0.00, 0.00, 5.00, 1.00},
                                            {0.00, 0.00, 0.00, 0.00, 5.00}};
    return A;
  }

  static std::vector<std::vector<scalar_t>> get_5x5_lt_fixture() {
    const auto KZ                        = KEEP_ZERO<scalar_t>();
    std::vector<std::vector<scalar_t>> A = {{5.00, KZ, 0.00, 0.00, 0.00},
                                            {2.00, 5.00, 0.00, 0.00, 0.00},
                                            {1.00, KZ, 5.00, 0.00, 0.00},
                                            {0.00, 0.00, 1.00, 5.00, 0.00},
                                            {KZ, 1.00, 1.00, 1.00, 5.00}};
    return A;
  }

  static std::vector<std::vector<scalar_t>> get_5x5_lt_ones_fixture() {
    std::vector<std::vector<scalar_t>> A = {{1.00, 0.00, 0.00, 0.00, 0.00},
                                            {0.00, 1.00, 0.00, 0.00, 0.00},
                                            {1.00, 0.00, 1.00, 0.00, 0.00},
                                            {0.00, 0.00, 1.00, 1.00, 0.00},
                                            {0.00, 1.00, 1.00, 1.00, 1.00}};
    return A;
  }

  static std::vector<std::vector<scalar_t>> get_6x6_lt_ones_fixture() {
    std::vector<std::vector<scalar_t>> A = {{1.00, 0.00, 0.00, 0.00, 0.00, 0.00}, {1.00, 1.00, 0.00, 0.00, 0.00, 0.00},
                                            {0.00, 0.00, 1.00, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 1.00, 0.00, 0.00},
                                            {0.00, 0.00, 0.00, 1.00, 1.00, 0.00}, {0.00, 1.00, 1.00, 1.00, 1.00, 1.00}};
    return A;
  }

  static bool do_cusparse() {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    return (std::is_same<size_type, int>::value && std::is_same<lno_t, int>::value &&
            std::is_same<typename device::execution_space, Kokkos::Cuda>::value);
#else
    return false;
#endif
  }

  struct ReductionCheck {
    ValuesType lhs;

    ReductionCheck(const ValuesType &lhs_) : lhs(lhs_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(lno_t i, scalar_t &tsum) const { tsum += lhs(i); }
  };

  static std::tuple<Crs, ValuesType, ValuesType> create_crs_lhs_rhs(const std::vector<std::vector<scalar_t>> &fixture) {
    RowMapType row_map;
    EntriesType entries;
    ValuesType values;

    compress_matrix(row_map, entries, values, fixture);
    const auto nrows = row_map.size() - 1;
    const auto nnz   = values.size();

    // Create known_lhs, generate rhs, then solve for lhs to compare to
    // known_lhs
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, scalar_t(1));

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    Crs triMtx("triMtx", nrows, nrows, nnz, values, row_map, entries);
    KokkosSparse::spmv("N", scalar_t(1), triMtx, known_lhs, scalar_t(0), rhs);

    return std::make_tuple(triMtx, lhs, rhs);
  }

  template <typename SpMatrix>
  static void basic_check(const SpMatrix &triMtx, const ValuesType &lhs, const ValuesType &rhs, const bool is_lower,
                          const size_type block_size = 0) {
    // FIXME Issues with some integral type combos for SEQLVLSCHED_TP2,
    // currently unavailable
    std::vector<SPTRSVAlgorithm> algs = {SPTRSVAlgorithm::SEQLVLSCHD_RP, SPTRSVAlgorithm::SEQLVLSCHD_TP1};
    if (block_size == 0) {
      // SEQLVLSCHD_TP1CHAIN and SPTRSV_CUSPARSE are not supported for blocks
      algs.push_back(SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN);
      if (do_cusparse()) {
        algs.push_back(SPTRSVAlgorithm::SPTRSV_CUSPARSE);
      }
    }

    auto row_map = triMtx.graph.row_map;
    auto entries = triMtx.graph.entries;
    auto values  = triMtx.values;

    const size_type nrows = row_map.size() - 1;

    for (auto alg : algs) {
      KernelHandle kh;
      kh.create_sptrsv_handle(alg, nrows, is_lower, block_size);
      if (alg == SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN) {
        auto chain_threshold = 1;
        kh.get_sptrsv_handle()->reset_chain_threshold(chain_threshold);
      }

      sptrsv_symbolic(&kh, row_map, entries, values);
      Kokkos::fence();

      sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
      Kokkos::fence();

      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(range_policy_t(0, lhs.extent(0)), ReductionCheck(lhs), sum);
      EXPECT_EQ(sum, lhs.extent(0));

      Kokkos::deep_copy(lhs, scalar_t(0));

      kh.destroy_sptrsv_handle();
    }
  }

  static void run_test_sptrsv() {
    const size_type nrows = 5;

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    using host_crsmat_t = typename KernelHandle::SPTRSVHandleType::host_crsmat_t;
    using host_graph_t  = typename host_crsmat_t::StaticCrsGraphType;

    using row_map_view_t = typename host_graph_t::row_map_type::non_const_type;
    using cols_view_t    = typename host_graph_t::entries_type::non_const_type;
    using values_view_t  = typename host_crsmat_t::values_type::non_const_type;

    // L & U handle for supernodal SpTrsv
    KernelHandle khL;
    KernelHandle khU;

    // right-hand-side and solution
    ValuesType B("rhs", nrows);
    ValuesType X("sol", nrows);

    // host CRS for L & U
    host_crsmat_t L, U, Ut;
#endif

    // Upper tri
    {
      {
        const auto [triMtx, lhs, rhs] = create_crs_lhs_rhs(get_5x5_ut_ones_fixture());

        basic_check(triMtx, lhs, rhs, false);
      }

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
      const scalar_t FIVE    = scalar_t(5);
      const size_type nnz_sp = 14;
      {
        // U in csr
        auto ut_fixture = get_5x5_ut_fixture();
        row_map_view_t hUrowptr;
        cols_view_t hUcolind;
        values_view_t hUvalues;

        // first row -> first supernode
        // second row -> first supernode
        // third row -> second supernode
        // fourth row -> third supernode
        // fifth row -> fourth supernode

        compress_matrix(hUrowptr, hUcolind, hUvalues, ut_fixture);

        // save U for Supernodal Sptrsv
        host_graph_t static_graph(hUcolind, hUrowptr);
        U = host_crsmat_t("CrsMatrixU", nrows, hUvalues, static_graph);

        // create handle for Supernodal Sptrsv
        bool is_lower_tri = false;
        khU.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, is_lower_tri);

        // X = U*ONES to generate B = A*ONES (on device)
        {
          RowMapType Urowptr("Urowptr", nrows + 1);
          EntriesType Ucolind("Ucolind", nnz_sp);
          ValuesType Uvalues("Uvalues", nnz_sp);

          Kokkos::deep_copy(Urowptr, hUrowptr);
          Kokkos::deep_copy(Ucolind, hUcolind);
          Kokkos::deep_copy(Uvalues, hUvalues);

          Crs mtxU("mtxU", nrows, nrows, nnz_sp, Uvalues, Urowptr, Ucolind);
          Kokkos::deep_copy(B, scalar_t(1));
          KokkosSparse::spmv("N", scalar_t(1), mtxU, B, scalar_t(0), X);
        }
      }

      {
        // U in csc (for inverting off-diag)
        row_map_view_t hUcolptr("hUcolptr", nrows + 1);
        cols_view_t hUrowind("hUrowind", nnz_sp);
        values_view_t hUvalues("hUvalues", nnz_sp);

        // The unsorted ordering seems to matter here, so we cannot use our
        // fixture tools.

        hUcolptr(0) = 0;
        hUcolptr(1) = 2;
        hUcolptr(2) = 4;
        hUcolptr(3) = 7;
        hUcolptr(4) = 9;
        hUcolptr(5) = 14;

        // colind
        // first column (first supernode)
        hUrowind(0) = 0;
        hUrowind(1) = 1;
        // second column (first supernode)
        hUrowind(2) = 0;
        hUrowind(3) = 1;
        // third column (second supernode)
        hUrowind(4) = 2;
        hUrowind(5) = 0;
        hUrowind(6) = 1;
        // fourth column (third supernode)
        hUrowind(7) = 3;
        hUrowind(8) = 2;
        // fifth column (fourth supernode)
        hUrowind(9)  = 4;
        hUrowind(10) = 0;
        hUrowind(11) = 1;
        hUrowind(12) = 2;
        hUrowind(13) = 3;

        // values
        // first column (first supernode)
        hUvalues(0) = FIVE;
        hUvalues(1) = scalar_t(0);
        // second column (first supernode)
        hUvalues(2) = scalar_t(1);
        hUvalues(3) = FIVE;
        // third column (second supernode)
        hUvalues(4) = FIVE;
        hUvalues(5) = scalar_t(1);
        hUvalues(6) = scalar_t(0);
        // fourth column (third supernode)
        hUvalues(7) = FIVE;
        hUvalues(8) = scalar_t(1);
        // fifth column (fourth supernode)
        hUvalues(9)  = FIVE;
        hUvalues(10) = scalar_t(0);
        hUvalues(11) = scalar_t(1);
        hUvalues(12) = scalar_t(1);
        hUvalues(13) = scalar_t(1);

        // store Ut in crsmat
        host_graph_t static_graph(hUrowind, hUcolptr);
        Ut = host_crsmat_t("CrsMatrixUt", nrows, hUvalues, static_graph);
      }
#endif
    }

    // Lower tri
    {
      {
        const auto [triMtx, lhs, rhs] = create_crs_lhs_rhs(get_5x5_lt_ones_fixture());

        basic_check(triMtx, lhs, rhs, true);
      }

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
      {
        // L in csc
        const size_type nnz_sp = 14;

        // first column (first supernode)
        // second column (first supernode)
        // third column (second supernode)
        // fourth column (third supernode)
        // fifth column (fourth supernode)

        auto lt_fixture = get_5x5_lt_fixture();
        row_map_view_t hLcolptr;
        cols_view_t hLrowind;
        values_view_t hLvalues;
        compress_matrix<true>(hLcolptr, hLrowind, hLvalues, lt_fixture);

        // store Lt in crsmat
        host_graph_t static_graph(hLrowind, hLcolptr);
        L = host_crsmat_t("CrsMatrixL", nrows, hLvalues, static_graph);

        bool is_lower_tri = true;
        khL.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, is_lower_tri);

        // generate B = A*ONES = L*(U*ONES), where X = U*ONES (on device)
        {
          RowMapType Lcolptr("Lcolptr", nrows + 1);
          EntriesType Lrowind("Lrowind", nnz_sp);
          ValuesType Lvalues("Lvalues", nnz_sp);

          Kokkos::deep_copy(Lcolptr, hLcolptr);
          Kokkos::deep_copy(Lrowind, hLrowind);
          Kokkos::deep_copy(Lvalues, hLvalues);

          Crs mtxL("mtxL", nrows, nrows, nnz_sp, Lvalues, Lcolptr, Lrowind);
          KokkosSparse::spmv("T", scalar_t(1), mtxL, X, scalar_t(0), B);
        }
      }

      {
        // unit-test for supernode SpTrsv (default)
        // > set up supernodes (block size = one)
        size_type nsupers = 4;
        Kokkos::View<int *, Kokkos::HostSpace> supercols("supercols", 1 + nsupers);
        supercols(0) = 0;
        supercols(1) = 2;     // two columns
        supercols(2) = 3;     // one column
        supercols(3) = 4;     // one column
        supercols(4) = 5;     // one column
        int *etree   = NULL;  // we generate graph internally

        // invert diagonal blocks
        bool invert_diag = true;
        khL.set_sptrsv_invert_diagonal(invert_diag);
        khU.set_sptrsv_invert_diagonal(invert_diag);

        // > symbolic (on host)
        sptrsv_supernodal_symbolic(nsupers, supercols.data(), etree, L.graph, &khL, U.graph, &khU);
        // > numeric (on host)
        sptrsv_compute(&khL, L);
        sptrsv_compute(&khU, U);
        Kokkos::fence();

        // > solve
        ValuesType b("b", nrows);
        Kokkos::deep_copy(b, B);
        Kokkos::deep_copy(X, scalar_t(0));
        sptrsv_solve(&khL, &khU, X, b);
        Kokkos::fence();

        // > check
        scalar_t sum = 0.0;
        Kokkos::parallel_reduce(range_policy_t(0, X.extent(0)), ReductionCheck(X), sum);
        EXPECT_EQ(sum, X.extent(0));

        khL.destroy_sptrsv_handle();
        khU.destroy_sptrsv_handle();
      }

      {
        // unit-test for supernode SpTrsv (running TRMM on device for compute)
        // > set up supernodes
        size_type nsupers = 4;
        Kokkos::View<int *, Kokkos::HostSpace> supercols("supercols", 1 + nsupers);
        supercols(0) = 0;
        supercols(1) = 2;     // two columns
        supercols(2) = 3;     // one column
        supercols(3) = 4;     // one column
        supercols(4) = 5;     // one column
        int *etree   = NULL;  // we generate tree internally

        // > create handles
        KernelHandle khLd;
        KernelHandle khUd;
        khLd.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, true);
        khUd.create_sptrsv_handle(SPTRSVAlgorithm::SUPERNODAL_DAG, nrows, false);

        // > invert diagonal blocks
        bool invert_diag = true;
        khLd.set_sptrsv_invert_diagonal(invert_diag);
        khUd.set_sptrsv_invert_diagonal(invert_diag);

        // > invert off-diagonal blocks
        bool invert_offdiag = true;
        khUd.set_sptrsv_column_major(true);
        khLd.set_sptrsv_invert_offdiagonal(invert_offdiag);
        khUd.set_sptrsv_invert_offdiagonal(invert_offdiag);

        // > forcing sptrsv compute to perform TRMM on device
        khLd.set_sptrsv_diag_supernode_sizes(1, 1);
        khUd.set_sptrsv_diag_supernode_sizes(1, 1);

        // > symbolic (on host)
        sptrsv_supernodal_symbolic(nsupers, supercols.data(), etree, L.graph, &khLd, Ut.graph, &khUd);
        // > numeric (on host)
        sptrsv_compute(&khLd, L);
        sptrsv_compute(&khUd, Ut);
        Kokkos::fence();

        // > solve
        ValuesType b("b", nrows);
        Kokkos::deep_copy(b, B);
        Kokkos::deep_copy(X, scalar_t(0));
        sptrsv_solve(&khLd, &khUd, X, b);
        Kokkos::fence();

        // > check
        scalar_t sum = 0.0;
        Kokkos::parallel_reduce(range_policy_t(0, X.extent(0)), ReductionCheck(X), sum);
        EXPECT_EQ(sum, X.extent(0));

        khLd.destroy_sptrsv_handle();
        khUd.destroy_sptrsv_handle();
      }
#endif
    }
  }

  static void run_test_sptrsv_blocks_impl(const bool is_lower, const size_type block_size) {
    auto fixture                      = is_lower ? get_6x6_lt_ones_fixture() : get_6x6_ut_ones_fixture();
    const auto [triMtx_crs, lhs, rhs] = create_crs_lhs_rhs(fixture);

    Bsr triMtx(triMtx_crs, block_size);
    basic_check(triMtx, lhs, rhs, is_lower, block_size);
  }

  static void run_test_sptrsv_blocks() {
    for (size_type block_size : {1, 2, 3}) {
      run_test_sptrsv_blocks_impl(true, block_size);
      run_test_sptrsv_blocks_impl(false, block_size);
    }
  }

  static void run_test_sptrsv_streams(SPTRSVAlgorithm test_algo, int nstreams, const bool is_lower) {
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

    const size_type nrows = 5;
    const size_type nnz   = 10;

    auto instances = Kokkos::Experimental::partition_space(execution_space(), std::vector<int>(nstreams, 1));

    std::vector<KernelHandle> kh_v(nstreams);
    std::vector<KernelHandle *> kh_ptr_v(nstreams);
    std::vector<RowMapType> row_map_v(nstreams);
    std::vector<EntriesType> entries_v(nstreams);
    std::vector<ValuesType> values_v(nstreams);
    std::vector<ValuesType> rhs_v(nstreams);
    std::vector<ValuesType> lhs_v(nstreams);

    auto fixture                  = is_lower ? get_5x5_lt_ones_fixture() : get_5x5_ut_ones_fixture();
    const auto [triMtx, lhs, rhs] = create_crs_lhs_rhs(fixture);

    auto row_map = triMtx.graph.row_map;
    auto entries = triMtx.graph.entries;
    auto values  = triMtx.values;

    for (int i = 0; i < nstreams; i++) {
      // Allocate
      row_map_v[i] = RowMapType("row_map", nrows + 1);
      entries_v[i] = EntriesType("entries", nnz);
      values_v[i]  = ValuesType("values", nnz);

      // Copy
      Kokkos::deep_copy(row_map_v[i], row_map);
      Kokkos::deep_copy(entries_v[i], entries);
      Kokkos::deep_copy(values_v[i], values);

      // Create known_lhs, generate rhs, then solve for lhs to compare to
      // known_lhs
      ValuesType known_lhs("known_lhs", nrows);
      // Create known solution lhs set to all 1's
      Kokkos::deep_copy(known_lhs, scalar_t(1));

      // Solution to find
      lhs_v[i] = ValuesType("lhs", nrows);

      // A*known_lhs generates rhs: rhs is dense, use spmv
      rhs_v[i] = ValuesType("rhs", nrows);

      KokkosSparse::spmv("N", scalar_t(1), triMtx, known_lhs, scalar_t(0), rhs_v[i]);
      Kokkos::fence();

      // Create handle
      kh_v[i] = KernelHandle();
      kh_v[i].create_sptrsv_handle(test_algo, nrows, is_lower);
      kh_ptr_v[i] = &kh_v[i];

      // Symbolic phase
      sptrsv_symbolic(kh_ptr_v[i], row_map_v[i], entries_v[i], values_v[i]);
      Kokkos::fence();
    }  // Done handle creation and sptrsv_symbolic on all streams

    // Solve phase
    sptrsv_solve_streams(instances, kh_ptr_v, row_map_v, entries_v, values_v, rhs_v, lhs_v);

    for (int i = 0; i < nstreams; i++) instances[i].fence();

    // Checking
    for (int i = 0; i < nstreams; i++) {
      scalar_t sum = 0.0;
      Kokkos::parallel_reduce(range_policy_t(0, lhs_v[i].extent(0)), ReductionCheck(lhs_v[i]), sum);
      EXPECT_EQ(sum, lhs_v[i].extent(0));
      kh_v[i].destroy_sptrsv_handle();
    }
  }
};

}  // namespace Test

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_sptrsv() {
  using TestStruct = Test::SptrsvTest<scalar_t, lno_t, size_type, device>;
  TestStruct::run_test_sptrsv();
  TestStruct::run_test_sptrsv_blocks();
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_sptrsv_streams() {
  using TestStruct                  = Test::SptrsvTest<scalar_t, lno_t, size_type, device>;
  std::vector<SPTRSVAlgorithm> algs = {SPTRSVAlgorithm::SEQLVLSCHD_RP, SPTRSVAlgorithm::SEQLVLSCHD_TP1};
  if (TestStruct::do_cusparse()) {
    algs.push_back(SPTRSVAlgorithm::SPTRSV_CUSPARSE);
  }

  for (auto alg : algs) {
    for (int nstreams = 1; nstreams <= 4; ++nstreams) {
      TestStruct::run_test_sptrsv_streams(alg, nstreams, true);
      TestStruct::run_test_sptrsv_streams(alg, nstreams, false);
    }
  }
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                      \
  TEST_F(TestCategory, sparse##_##sptrsv##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_sptrsv<SCALAR, ORDINAL, OFFSET, DEVICE>();                                      \
    test_sptrsv_streams<SCALAR, ORDINAL, OFFSET, DEVICE>();                              \
  }

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
