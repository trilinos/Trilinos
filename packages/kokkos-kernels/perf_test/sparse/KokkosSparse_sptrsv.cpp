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

#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>
#include <unordered_map>

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include <cusparse.h>
#endif

#include <Kokkos_Core.hpp>

#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_sptrsv.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_default_types.hpp"
#include <KokkosKernels_IOUtils.hpp>
#include "KokkosSparse_IOUtils.hpp"

// #define INTERNAL_CUSPARSE

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA) && (!defined(KOKKOS_ENABLE_CUDA) || (8000 <= CUDA_VERSION))
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

// #define PRINTVIEWSSPTRSVPERF
// #define PRINT_HLEVEL_FREQ_PLOT
// #define PRINT_LEVEL_LIST

enum {
  DEFAULT,
  CUSPARSE,
  LVLSCHED_RP,
  LVLSCHED_TP1,
  /*LVLSCHED_TP2,*/ LVLSCHED_TP1CHAIN,
  CUSPARSE_K
};

#ifdef PRINTVIEWSSPTRSVPERF
template <class ViewType>
void print_view1d(const ViewType dv) {
  auto v = Kokkos::create_mirror_view(dv);
  Kokkos::deep_copy(v, dv);
  std::cout << "Output for view " << v.label() << std::endl;
  for (size_t i = 0; i < v.extent(0); ++i) {
    std::cout << "v(" << i << ") = " << v(i) << " , ";
  }
  std::cout << std::endl;
}
#else
template <class ViewType>
void print_view1d(const ViewType /*dv*/) {}
#endif

template <class RowMapType, class EntriesType>
void check_entries_sorted(const RowMapType drow_map, const EntriesType dentries) {
  auto row_map = Kokkos::create_mirror_view(drow_map);
  Kokkos::deep_copy(row_map, drow_map);
  auto entries = Kokkos::create_mirror_view(dentries);
  Kokkos::deep_copy(entries, dentries);

  for (size_t row = 0; row < row_map.extent(0) - 1; ++row) {
    size_t start = row_map(row);
    size_t end   = row_map(row + 1);
    for (size_t offset = start; offset < end - 1; ++offset) {
      size_t pcol = entries(offset);
      size_t ncol = entries(offset + 1);
      if (pcol > ncol) {
        std::cout << "  UNSORTED!!" << std::endl;
      }
    }
  }
}

int test_sptrsv_perf(std::vector<int> tests, const std::string &lfilename, const std::string &ufilename,
                     const int team_size, const int vector_length, const int /*idx_offset*/, const int loop,
                     const int chain_threshold = 0, const float /*dense_row_percent*/ = -1.0) {
  typedef default_scalar scalar_t;
  typedef default_lno_t lno_t;
  typedef default_size_type size_type;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef typename execution_space::memory_space memory_space;

  typedef KokkosSparse::CrsMatrix<scalar_t, lno_t, execution_space, void, size_type> crsmat_t;
  typedef typename crsmat_t::StaticCrsGraphType graph_t;

  typedef Kokkos::View<scalar_t *, memory_space> ValuesType;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, execution_space, memory_space,
                                                           memory_space>
      KernelHandle;

  const scalar_t ZERO = scalar_t(0);
  const scalar_t ONE  = scalar_t(1);

  // Read lmtx
  // Run all requested algorithms
  // Read umtx
  // Run all requested algorithms

  // LOWERTRI
  std::cout << "\n\n" << std::endl;
  if (!lfilename.empty()) {
    std::cout << "Lower Tri Begin: Read matrix filename " << lfilename << std::endl;
    crsmat_t triMtx       = KokkosSparse::Impl::read_kokkos_crst_matrix<crsmat_t>(lfilename.c_str());  // in_matrix
    graph_t graph         = triMtx.graph;                                                              // in_graph
    const size_type nrows = graph.numRows();

    // Create the rhs and lhs_known solution
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    std::cout << "SPMV" << std::endl;
    KokkosSparse::spmv("N", ONE, triMtx, known_lhs, ZERO, rhs);

    auto row_map = graph.row_map;
    auto entries = graph.entries;
    auto values  = triMtx.values;

    std::cout << "Lower Perf: row_map.extent(0) = " << row_map.extent(0) << std::endl;
    std::cout << "Lower Perf: entries.extent(0) = " << entries.extent(0) << std::endl;
    std::cout << "Lower Perf: values.extent(0) = " << values.extent(0) << std::endl;

    std::cout << "Lower Perf: lhs.extent(0) = " << lhs.extent(0) << std::endl;
    std::cout << "Lower Perf: rhs.extent(0) = " << rhs.extent(0) << std::endl;

    check_entries_sorted(row_map, entries);

#ifdef PRINTVIEWSSPTRSVPERF
    print_view1d(row_map);
    print_view1d(entries);
    print_view1d(values);
    print_view1d(known_lhs);
    print_view1d(rhs);
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
    // std::cout << "  cusparse: create handle" << std::endl;
    cusparseStatus_t status;
    cusparseHandle_t handle = 0;
    status                  = cusparseCreate(&handle);
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "handle create status error name " << (status) << std::endl;
    cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
    cusparseMatDescr_t descr = 0;
    csrsv2Info_t info        = 0;
    int pBufferSize;
    void *pBuffer = 0;
    int structural_zero;
    int numerical_zero;
    const double alpha                 = 1.;
    const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
    const cusparseOperation_t trans    = CUSPARSE_OPERATION_NON_TRANSPOSE;

    // step 1: create a descriptor which contains
    // - matrix L is lower triangular
    //   (L may not have all diagonal elements.)
    status = cusparseCreateMatDescr(&descr);
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "matdescr create status error name " << (status) << std::endl;
    // cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ONE);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_LOWER);
    cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    // cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_UNIT);

    // step 2: create a empty info structure
    // std::cout << "  cusparse: create csrsv2info" << std::endl;
    status = cusparseCreateCsrsv2Info(&info);
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "csrsv2info create status error name " << (status) << std::endl;

    // step 3: query how much memory used in csrsv2, and allocate the buffer
    int nnz = triMtx.nnz();
    cusparseDcsrsv2_bufferSize(handle, trans, nrows, nnz, descr, values.data(), row_map.data(), entries.data(), info,
                               &pBufferSize);
    // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
    cudaMalloc((void **)&pBuffer, pBufferSize);
#endif

    for (auto test : tests) {
      std::cout << "\ntest = " << test << std::endl;

      KernelHandle kh;
      bool is_lower_tri = true;

      std::cout << "Create handle (lower)" << std::endl;
      switch (test) {
        case LVLSCHED_RP:
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_RP, nrows, is_lower_tri);
          kh.get_sptrsv_handle()->print_algorithm();
          break;
        case LVLSCHED_TP1:
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
          std::cout << "TP1 set team_size = " << team_size << std::endl;
          if (team_size != -1) kh.get_sptrsv_handle()->set_team_size(team_size);
          kh.get_sptrsv_handle()->print_algorithm();
          break;
        case LVLSCHED_TP1CHAIN:
          printf("TP1 with CHAIN\n");
          printf("chain_threshold %d\n", chain_threshold);
          printf("team_size %d\n", team_size);
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN, nrows, is_lower_tri);
          kh.get_sptrsv_handle()->reset_chain_threshold(chain_threshold);
          if (team_size != -1) kh.get_sptrsv_handle()->set_team_size(team_size);
          if (vector_length != -1) kh.get_sptrsv_handle()->set_vector_size(vector_length);
          kh.get_sptrsv_handle()->print_algorithm();
          break;
          /*
                case LVLSCHED_TP2:
                  kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHED_TP2,
             nrows, is_lower_tri); if (team_size != -1)
             kh.get_sptrsv_handle()->set_team_size(team_size); if (vector_length
             != -1) kh.get_sptrsv_handle()->set_vector_size(vector_length);
                  kh.get_sptrsv_handle()->print_algorithm();
                  break;
          */
        case CUSPARSE_K:
          printf("CUSPARSE WRAPPER\n");
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SPTRSV_CUSPARSE, nrows, is_lower_tri);
          kh.get_sptrsv_handle()->print_algorithm();
          break;
        case CUSPARSE:
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
          std::cout << "CUSPARSE: No kk interface added yet" << std::endl;
          // cusparse_matvec(A, x, y, rows_per_thread, team_size,
          // vector_length);
          break;
#else
          std::cout << "CUSPARSE not enabled: Fall through to defaults" << std::endl;
          [[fallthrough]];
#endif
        default:
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
          if (team_size != -1) kh.get_sptrsv_handle()->set_team_size(team_size);
          kh.get_sptrsv_handle()->print_algorithm();
      }

      // Init run to clear cache etc.
      Kokkos::Timer timer;
      if (test != CUSPARSE) {
        timer.reset();
        if (test == CUSPARSE_K) {
          printf("cusparsek symbolic\n");
          sptrsv_symbolic(&kh, row_map, entries, values);
          printf("  finished cusparsek symbolic\n");
        } else {
          sptrsv_symbolic(&kh, row_map, entries);
        }
        std::cout << "LTRI Symbolic Time: " << timer.seconds() << std::endl;

        // std::cout << "TriSolve Solve" << std::endl;
        timer.reset();
        sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
        Kokkos::fence();
        std::cout << "LTRI Solve Time: " << timer.seconds() << std::endl;

      }
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
      // step 4: perform analysis
      else {
        // int nnz = triMtx.nnz();
        // std::cout << "  cusparse path: analysis" << std::endl;
        // status = cusparseDcsrsv2_analysis(handle, trans, nrows, nnz, descr,
        // (double*)dvalues, (int *)drow_map, (int *)dentries, info, policy,
        // pBuffer);
        timer.reset();
        status = cusparseDcsrsv2_analysis(handle, trans, nrows, triMtx.nnz(), descr, values.data(), row_map.data(),
                                          entries.data(), info, policy, pBuffer);
        std::cout << "LTRI Cusparse Symbolic Time: " << timer.seconds() << std::endl;
        if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "analysis status error name " << (status) << std::endl;
        // L has unit diagonal, so no structural zero is reported.

        // std::cout << "  cusparse path: analysis" << std::endl;
        status = cusparseXcsrsv2_zeroPivot(handle, info, &structural_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
          printf("L(%d,%d) is missing\n", structural_zero, structural_zero);
        }

        // step 5: solve L*y = x
        // std::cout << "  cusparse path: solve" << std::endl;
        // status = cusparseDcsrsv2_solve(handle, trans, nrows, nnz, &alpha,
        // descr, (double*)dvalues, (int *)drow_map, (int *)dentries, info,
        // (double*)drhs, (double*)dlhs, policy, pBuffer);
        timer.reset();
        status = cusparseDcsrsv2_solve(handle, trans, nrows, triMtx.nnz(), &alpha, descr, values.data(), row_map.data(),
                                       entries.data(), info, rhs.data(), lhs.data(), policy, pBuffer);
        Kokkos::fence();
        std::cout << "LTRI Cusparse Solve Time: " << timer.seconds() << std::endl;
        if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "solve status error name " << (status) << std::endl;
        // L has unit diagonal, so no numerical zero is reported.
        status = cusparseXcsrsv2_zeroPivot(handle, info, &numerical_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
          printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
        }
      }
#endif
      // Error Check
      Kokkos::fence();
      {
        scalar_t sum = 0.0;
        Kokkos::parallel_reduce(
            Kokkos::RangePolicy<execution_space>(0, lhs.extent(0)),
            KOKKOS_LAMBDA(const lno_t i, scalar_t &tsum) { tsum += (known_lhs(i) - lhs(i)) * (known_lhs(i) - lhs(i)); },
            sum);

        scalar_t norm_ssd = sqrt(sum / lhs.extent(0));
        std::cout << "  ssd = " << sum << "  norm_sqrt_ssd = " << norm_ssd << std::endl;

        if (norm_ssd > 1e-8) {
          std::cout << "Lower Tri Solve FAILURE: norm_ssd = " << norm_ssd << std::endl;
          return 1;
        } else {
          std::cout << "\nLower Tri Solve Init Test: SUCCESS!\n" << std::endl;
        }

        /*
      sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<execution_space>(0,
      lhs.extent(0)), KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) { tsum +=
      lhs(i);
        }, sum);

      if ( sum != lhs.extent(0) ) {
        std::cout << "Lower Tri Solve FAILURE: sum = " << sum << std::endl;
        auto hsoln = Kokkos::create_mirror_view(lhs);
        Kokkos::deep_copy(hsoln, lhs);
        for ( size_t i = 0; i < hsoln.extent(0); ++i ) {
          std::cout << "lhs(" << i << ") = " << hsoln(i) << std::endl;
        }
        return 1;
      }
      else {
       std::cout << "\nLower Tri Solve Init Test: SUCCESS!\n" << std::endl;
      }
        */
      }

      // Benchmark
      Kokkos::fence();
      double min_time = 1.0e32;
      double max_time = 0.0;
      double ave_time = 0.0;

      for (int titer = 0; titer < loop; titer++) {
        timer.reset();

        if (test != CUSPARSE) {
#ifdef CHECKALLRUNRESULTS
          Kokkos::deep_copy(lhs, 0, 0);
#endif
          sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
#ifdef CHECKALLRUNRESULTS
          {
            scalar_t sum = 0.0;
            Kokkos::parallel_reduce(
                Kokkos::RangePolicy<execution_space>(0, lhs.extent(0)),
                KOKKOS_LAMBDA(const lno_t i, scalar_t &tsum) {
                  tsum += (known_lhs(i) - lhs(i)) * (known_lhs(i) - lhs(i));
                },
                sum);

            scalar_t norm_ssd = sqrt(sum / lhs.extent(0));
            std::cout << "  ssd = " << sum << "  norm_sqrt_ssd = " << norm_ssd << std::endl;
            if (norm_ssd > 1e-8) {
              std::cout << "Lower Tri Solve FAILURE: norm_ssd = " << norm_ssd << std::endl;
              return 1;
            } else {
              std::cout << "\nLower Tri Solve Init Test: SUCCESS!\n" << std::endl;
            }
          }
#endif
        }
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
        else {
          cusparseDcsrsv2_solve(handle, trans, nrows, triMtx.nnz(), &alpha, descr, values.data(), row_map.data(),
                                entries.data(), info, rhs.data(), lhs.data(), policy, pBuffer);
        }
#endif

        Kokkos::fence();
        double time = timer.seconds();
        ave_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }

      std::cout << "LOOP_AVG_TIME:  " << ave_time / loop << std::endl;
      std::cout << "LOOP_MAX_TIME:  " << max_time << std::endl;
      std::cout << "LOOP_MIN_TIME:  " << min_time << std::endl;

// Output for level frequency plot
#ifdef PRINT_HLEVEL_FREQ_PLOT
      if (test != CUSPARSE) {
        auto hnpl              = kh.get_sptrsv_handle()->get_host_nodes_per_level();
        auto nlevels           = kh.get_sptrsv_handle()->get_num_levels();
        std::string algmstring = kh.get_sptrsv_handle()->return_algorithm_string();
        std::cout << algmstring << std::endl;
        // Create filename
        std::string filename = "lower_nodes_per_level_" + algmstring + ".txt";
        std::cout << filename << std::endl;
        std::cout << "  nlevels = " << nlevels << std::endl;
        std::ofstream outfile;
        outfile.open(filename);
        if (outfile.is_open()) {
          for (int i = 0; i < nlevels; ++i) {
            outfile << hnpl(i) << std::endl;
            // std::cout  << hnpl(i) << std::endl;
          }
          outfile.close();
        } else {
          std::cout << "OUTFILE DID NOT OPEN!!!" << std::endl;
        }

        auto hngpl = kh.get_sptrsv_handle()->get_host_nodes_grouped_by_level();
        filename   = "lower_nodes_groupby_level_" + algmstring + ".txt";
        std::cout << filename << std::endl;
        outfile.open(filename);
        if (outfile.is_open()) {
          for (size_t i = 0; i < hngpl.extent(0); ++i) outfile << hngpl(i) << std::endl;
          outfile.close();
        } else {
          std::cout << "OUTFILE DID NOT OPEN!!!" << std::endl;
        }
      }
#endif

#ifdef PRINT_LEVEL_LIST
      if (test != CUSPARSE) {
        auto level_list  = kh.get_sptrsv_handle()->get_level_list();
        auto hlevel_list = Kokkos::create_mirror_view(level_list);
        Kokkos::deep_copy(hlevel_list, level_list);

        auto nlevels = kh.get_sptrsv_handle()->get_num_levels();

        std::string algmstring = kh.get_sptrsv_handle()->return_algorithm_string();
        std::cout << algmstring << std::endl;
        // Create filename
        std::string filename = "lower_level_list_" + algmstring + ".txt";
        std::cout << filename << std::endl;
        std::cout << "  nlevels = " << nlevels << "  nodes = " << hlevel_list.extent(0) << std::endl;
        std::ofstream outfile;
        outfile.open(filename);
        if (outfile.is_open()) {
          for (size_t i = 0; i < hlevel_list.extent(0); ++i) outfile << hlevel_list(i) << std::endl;
          outfile.close();
        } else {
          std::cout << "OUTFILE DID NOT OPEN!!!" << std::endl;
        }
      }
#endif

      kh.destroy_sptrsv_handle();
    }

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
    // step 6: free resources
    cudaFree(pBuffer);
    cusparseDestroyCsrsv2Info(info);
    cusparseDestroyMatDescr(descr);
    cusparseDestroy(handle);
#endif
  }  // end lowertri
  Kokkos::fence();

  std::cout << "\n\n" << std::endl;
  // UPPERTRI
  if (!ufilename.empty()) {
    std::cout << "Upper Tri Begin: Read matrix filename " << ufilename << std::endl;
    crsmat_t triMtx       = KokkosSparse::Impl::read_kokkos_crst_matrix<crsmat_t>(ufilename.c_str());  // in_matrix
    graph_t graph         = triMtx.graph;                                                              // in_graph
    const size_type nrows = graph.numRows();

    // Create the rhs and lhs_known solution
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    std::cout << "SPMV" << std::endl;
    KokkosSparse::spmv("N", ONE, triMtx, known_lhs, ZERO, rhs);

    auto row_map = graph.row_map;
    auto entries = graph.entries;
    auto values  = triMtx.values;

    std::cout << "Upper Perf: row_map.extent(0) = " << row_map.extent(0) << std::endl;
    std::cout << "Upper Perf: entries.extent(0) = " << entries.extent(0) << std::endl;
    std::cout << "Upper Perf: values.extent(0) = " << values.extent(0) << std::endl;

    std::cout << "Upper Perf: lhs.extent(0) = " << lhs.extent(0) << std::endl;
    std::cout << "Upper Perf: rhs.extent(0) = " << rhs.extent(0) << std::endl;

    check_entries_sorted(row_map, entries);

#ifdef PRINTVIEWSSPTRSVPERF
    print_view1d(row_map);
    print_view1d(entries);
    print_view1d(values);
    print_view1d(known_lhs);
    print_view1d(rhs);
#endif

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
    // std::cout << "  cusparse: create handle" << std::endl;
    cusparseStatus_t status;
    cusparseHandle_t handle = 0;
    status                  = cusparseCreate(&handle);
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "handle create status error name " << (status) << std::endl;
    cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
    cusparseMatDescr_t descr = 0;
    csrsv2Info_t info        = 0;
    int pBufferSize;
    void *pBuffer = 0;
    int structural_zero;
    int numerical_zero;
    const double alpha                 = 1.;
    const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
    const cusparseOperation_t trans    = CUSPARSE_OPERATION_NON_TRANSPOSE;

    // step 1: create a descriptor which contains
    //   (L may not have all diagonal elements.)
    status = cusparseCreateMatDescr(&descr);
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "matdescr create status error name " << (status) << std::endl;
    // cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ONE);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_UPPER);
    cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    // cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_TRIANGULAR);
    // cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_UNIT);

    // step 2: create a empty info structure
    // std::cout << "  cusparse: create csrsv2info" << std::endl;
    status = cusparseCreateCsrsv2Info(&info);
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "csrsv2info create status error name " << (status) << std::endl;

    // step 3: query how much memory used in csrsv2, and allocate the buffer
    int nnz = triMtx.nnz();
    cusparseDcsrsv2_bufferSize(handle, trans, nrows, nnz, descr, values.data(), row_map.data(), entries.data(), info,
                               &pBufferSize);
    // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
    cudaMalloc((void **)&pBuffer, pBufferSize);
#endif

    for (auto test : tests) {
      std::cout << "\ntest = " << test << std::endl;

      KernelHandle kh;
      bool is_lower_tri = false;

      std::cout << "Create handle (upper)" << std::endl;
      switch (test) {
        case LVLSCHED_RP:
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_RP, nrows, is_lower_tri);
          kh.get_sptrsv_handle()->print_algorithm();
          break;
        case LVLSCHED_TP1:
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
          std::cout << "TP1 set team_size = " << team_size << std::endl;
          if (team_size != -1) kh.get_sptrsv_handle()->set_team_size(team_size);
          kh.get_sptrsv_handle()->print_algorithm();
          break;
        case LVLSCHED_TP1CHAIN:
          printf("TP1 with CHAIN\n");
          printf("chain_threshold %d\n", chain_threshold);
          printf("team_size %d\n", team_size);
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN, nrows, is_lower_tri);
          kh.get_sptrsv_handle()->reset_chain_threshold(chain_threshold);
          if (team_size != -1) kh.get_sptrsv_handle()->set_team_size(team_size);
          if (vector_length != -1) kh.get_sptrsv_handle()->set_vector_size(vector_length);
          kh.get_sptrsv_handle()->print_algorithm();
          break;
          /*
                case LVLSCHED_TP2:
                  kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHED_TP2,
             nrows, is_lower_tri); if (team_size != -1)
             kh.get_sptrsv_handle()->set_team_size(team_size); if (vector_length
             != -1) kh.get_sptrsv_handle()->set_vector_size(vector_length);
                  kh.get_sptrsv_handle()->print_algorithm();
                  break;
          */
        case CUSPARSE:
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
          std::cout << "CUSPARSE: No kk interface added yet" << std::endl;
          // cusparse_matvec(A, x, y, rows_per_thread, team_size,
          // vector_length);
          break;
#else
          std::cout << "CUSPARSE not enabled: Fall through to defaults" << std::endl;
          [[fallthrough]];
#endif
        default:
          kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
          if (team_size != -1) kh.get_sptrsv_handle()->set_team_size(team_size);
          kh.get_sptrsv_handle()->print_algorithm();
      }

      // Init run to clear cache etc.
      Kokkos::Timer timer;
      if (test != CUSPARSE) {
        timer.reset();
        sptrsv_symbolic(&kh, row_map, entries);
        std::cout << "UTRI Symbolic Time: " << timer.seconds() << std::endl;

        // std::cout << "TriSolve Solve" << std::endl;
        timer.reset();
        sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
        Kokkos::fence();
        std::cout << "UTRI Solve Time: " << timer.seconds() << std::endl;

      }
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
      // step 4: perform analysis
      else {
        // int nnz = triMtx.nnz();
        // std::cout << "  cusparse path: analysis" << std::endl;
        // status = cusparseDcsrsv2_analysis(handle, trans, nrows, nnz, descr,
        // (double*)dvalues, (int *)drow_map, (int *)dentries, info, policy,
        // pBuffer);
        timer.reset();
        status = cusparseDcsrsv2_analysis(handle, trans, nrows, triMtx.nnz(), descr, values.data(), row_map.data(),
                                          entries.data(), info, policy, pBuffer);
        std::cout << "UTRI Cusparse Symbolic Time: " << timer.seconds() << std::endl;
        if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "analysis status error name " << (status) << std::endl;
        // L has unit diagonal, so no structural zero is reported.

        status = cusparseXcsrsv2_zeroPivot(handle, info, &structural_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
          printf("L(%d,%d) is missing\n", structural_zero, structural_zero);
        }

        // step 5: solve L*y = x
        // std::cout << "  cusparse path: solve" << std::endl;
        // status = cusparseDcsrsv2_solve(handle, trans, nrows, nnz, &alpha,
        // descr, (double*)dvalues, (int *)drow_map, (int *)dentries, info,
        // (double*)drhs, (double*)dlhs, policy, pBuffer);
        timer.reset();
        status = cusparseDcsrsv2_solve(handle, trans, nrows, triMtx.nnz(), &alpha, descr, values.data(), row_map.data(),
                                       entries.data(), info, rhs.data(), lhs.data(), policy, pBuffer);
        Kokkos::fence();
        std::cout << "UTRI Cusparse Solve Time: " << timer.seconds() << std::endl;
        if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "solve status error name " << (status) << std::endl;
        // L has unit diagonal, so no numerical zero is reported.
        status = cusparseXcsrsv2_zeroPivot(handle, info, &numerical_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
          printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
        }
      }
#endif
      // Error Check
      Kokkos::fence();
      {
        scalar_t sum = 0.0;
        Kokkos::parallel_reduce(
            Kokkos::RangePolicy<execution_space>(0, lhs.extent(0)),
            KOKKOS_LAMBDA(const lno_t i, scalar_t &tsum) { tsum += (known_lhs(i) - lhs(i)) * (known_lhs(i) - lhs(i)); },
            sum);

        scalar_t norm_ssd = sqrt(sum / lhs.extent(0));
        std::cout << "  ssd = " << sum << "  norm_sqrt_ssd = " << norm_ssd << std::endl;

        if (norm_ssd > 1e-8) {
          std::cout << "Upper Tri Solve FAILURE: norm_ssd = " << norm_ssd << std::endl;
          return 1;
        } else {
          std::cout << "\nUpper Tri Solve Init Test: SUCCESS!\n" << std::endl;
        }

        /*
      sum = 0.0;
      Kokkos::parallel_reduce( Kokkos::RangePolicy<execution_space>(0,
      lhs.extent(0)), KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) { tsum +=
      lhs(i);
        }, sum);

      if ( sum != scalar_t(lhs.extent(0)) ) {
        std::cout << "Upper Tri Solve FAILURE: sum = " << sum << std::endl;
        auto hsoln = Kokkos::create_mirror_view(lhs);
        Kokkos::deep_copy(hsoln, lhs);
        for ( size_t i = 0; i < hsoln.extent(0); ++i ) {
          std::cout << "lhs(" << i << ") = " << hsoln(i) << std::endl;
        }
        return 1;
      }
      else {
       std::cout << "\nUpper Tri Solve Init Test: SUCCESS!\n" << std::endl;
      }
        */
      }

      // Benchmark
      Kokkos::fence();
      double min_time = 1.0e32;
      double max_time = 0.0;
      double ave_time = 0.0;

      for (int titer = 0; titer < loop; titer++) {
        timer.reset();

        if (test != CUSPARSE) {
#ifdef CHECKALLRUNRESULTS
          Kokkos::deep_copy(lhs, 0, 0);
#endif
          sptrsv_solve(&kh, row_map, entries, values, rhs, lhs);
#ifdef CHECKALLRUNRESULTS
          {
            scalar_t sum = 0.0;
            Kokkos::parallel_reduce(
                Kokkos::RangePolicy<execution_space>(0, lhs.extent(0)),
                KOKKOS_LAMBDA(const lno_t i, scalar_t &tsum) {
                  tsum += (known_lhs(i) - lhs(i)) * (known_lhs(i) - lhs(i));
                },
                sum);

            scalar_t norm_ssd = sqrt(sum / lhs.extent(0));
            std::cout << "  ssd = " << sum << "  norm_sqrt_ssd = " << norm_ssd << std::endl;
            if (norm_ssd > 1e-8) {
              std::cout << "Upper Tri Solve FAILURE: norm_ssd = " << norm_ssd << std::endl;
              return 1;
            } else {
              std::cout << "\nUpper Tri Solve Init Test: SUCCESS!\n" << std::endl;
            }
          }
#endif
        }
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
        else {
          cusparseDcsrsv2_solve(handle, trans, nrows, triMtx.nnz(), &alpha, descr, values.data(), row_map.data(),
                                entries.data(), info, rhs.data(), lhs.data(), policy, pBuffer);
        }
#endif

        Kokkos::fence();
        double time = timer.seconds();
        ave_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }

      std::cout << "LOOP_AVG_TIME:  " << ave_time / loop << std::endl;
      std::cout << "LOOP_MAX_TIME:  " << max_time << std::endl;
      std::cout << "LOOP_MIN_TIME:  " << min_time << std::endl;

// Output for level frequency plot
#ifdef PRINT_HLEVEL_FREQ_PLOT
      if (test != CUSPARSE) {
        auto hnpl              = kh.get_sptrsv_handle()->get_host_nodes_per_level();
        auto nlevels           = kh.get_sptrsv_handle()->get_num_levels();
        std::string algmstring = kh.get_sptrsv_handle()->return_algorithm_string();
        std::cout << algmstring << std::endl;
        // Create filename
        std::string filename = "upper_nodes_per_level_" + algmstring + ".txt";
        std::cout << filename << std::endl;
        std::cout << "  nlevels = " << nlevels << std::endl;
        std::ofstream outfile;
        outfile.open(filename);
        if (outfile.is_open()) {
          for (int i = 0; i < nlevels; ++i) {
            outfile << hnpl(i) << std::endl;
            // std::cout  << hnpl(i) << std::endl;
          }
          outfile.close();
        } else {
          std::cout << "OUTFILE DID NOT OPEN!!!" << std::endl;
        }

        auto hngpl = kh.get_sptrsv_handle()->get_host_nodes_grouped_by_level();
        filename   = "upper_nodes_groupby_level_" + algmstring + ".txt";
        std::cout << filename << std::endl;
        outfile.open(filename);
        if (outfile.is_open()) {
          for (size_t i = 0; i < hngpl.extent(0); ++i) outfile << hngpl(i) << std::endl;
          outfile.close();
        } else {
          std::cout << "OUTFILE DID NOT OPEN!!!" << std::endl;
        }
      }
#endif

#ifdef PRINT_LEVEL_LIST
      if (test != CUSPARSE) {
        auto level_list  = kh.get_sptrsv_handle()->get_level_list();
        auto hlevel_list = Kokkos::create_mirror_view(level_list);
        Kokkos::deep_copy(hlevel_list, level_list);

        auto nlevels = kh.get_sptrsv_handle()->get_num_levels();

        std::string algmstring = kh.get_sptrsv_handle()->return_algorithm_string();
        std::cout << algmstring << std::endl;
        // Create filename
        std::string filename = "upper_level_list_" + algmstring + ".txt";
        std::cout << filename << std::endl;
        std::cout << "  nlevels = " << nlevels << "  nodes = " << hlevel_list.extent(0) << std::endl;
        std::ofstream outfile;
        outfile.open(filename);
        if (outfile.is_open()) {
          for (size_t i = 0; i < hlevel_list.extent(0); ++i) outfile << hlevel_list(i) << std::endl;
          outfile.close();
        } else {
          std::cout << "OUTFILE DID NOT OPEN!!!" << std::endl;
        }
      }
#endif

      kh.destroy_sptrsv_handle();
    }

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && defined(INTERNAL_CUSPARSE)
    // step 6: free resources
    cudaFree(pBuffer);
    cusparseDestroyCsrsv2Info(info);
    cusparseDestroyMatDescr(descr);
    cusparseDestroy(handle);
#endif
  }  // end uppertri
  Kokkos::fence();

  return 0;
}

void print_help_sptrsv() {
  printf("Options:\n");
  printf("  --test [OPTION] : Use different kernel implementations\n");
  printf("                    Options:\n");
  printf(
      "                      lvlrp, lvltp1, lvltp2, lvltp1chain, lvldensetp1, "
      "lvldensetp2\n\n");
  printf("                      cusparse           (Vendor Libraries)\n\n");
  printf(
      "  -lf [file]      : Read in Matrix Market formatted text file "
      "'file'.\n");
  printf(
      "  -uf [file]      : Read in Matrix Market formatted text file "
      "'file'.\n");
  printf("  --offset [O]    : Subtract O from every index.\n");
  printf(
      "                    Useful in case the matrix market file is not 0 "
      "based.\n\n");
  printf("  -ts [T]         : Number of threads per team.\n");
  printf(
      "  -vl [V]         : Vector-length (i.e. how many Cuda threads are a "
      "Kokkos 'thread').\n");
  printf(
      "  -ct [V]         : Chain threshold: Only has effect of lvltp1chain "
      "algorithm.\n");
  printf(
      "  -dr [V]         : Dense row percent (as float): Only has effect of "
      "lvldensetp1 algorithm.\n");
  printf("  --loop [LOOP]   : How many spmv to run to aggregate average time. \n");
  //  printf("  --write-lvl-freq: Write output files with number of nodes per
  //  level for each matrix and algorithm.\n"); printf("  -s [N]          :
  //  generate a semi-random banded (band size 0.01xN) NxN matrix\n"); printf("
  //  with average of 10 entries per row.\n"); printf("  --schedule [SCH]: Set
  //  schedule for kk variant (static,dynamic,auto [ default ]).\n"); printf("
  //  -fb [file]      : Read in binary Matrix files 'file'.\n"); printf("
  //  --write-binary  : In combination with -f, generate binary files.\n");
}

int main(int argc, char **argv) {
  std::vector<int> tests;

  std::string lfilename;
  std::string ufilename;

  int vector_length       = -1;
  int team_size           = -1;
  int idx_offset          = 0;
  int loop                = 1;
  int chain_threshold     = 0;
  float dense_row_percent = -1.0;
  // int schedule=AUTO;

  if (argc == 1) {
    print_help_sptrsv();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "--test") == 0)) {
      i++;
      if ((strcmp(argv[i], "lvlrp") == 0)) {
        tests.push_back(LVLSCHED_RP);
      }
      if ((strcmp(argv[i], "lvltp1") == 0)) {
        tests.push_back(LVLSCHED_TP1);
      }
      if ((strcmp(argv[i], "lvltp1chain") == 0)) {
        tests.push_back(LVLSCHED_TP1CHAIN);
      }
      /*
      if((strcmp(argv[i],"lvltp2")==0)) {
        tests.push_back( LVLSCHED_TP2 );
      }
      */
      if ((strcmp(argv[i], "cusparse") == 0)) {
        tests.push_back(CUSPARSE);
      }
      if ((strcmp(argv[i], "cusparsek") == 0)) {
        tests.push_back(CUSPARSE_K);
      }
      continue;
    }
    if ((strcmp(argv[i], "-lf") == 0)) {
      lfilename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-uf") == 0)) {
      ufilename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-ts") == 0)) {
      team_size = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-vl") == 0)) {
      vector_length = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-ct") == 0)) {
      chain_threshold = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-dr") == 0)) {
      dense_row_percent = atof(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-l") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--offset") == 0)) {
      idx_offset = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--loop") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    /*
      if((strcmp(argv[i],"-lfb")==0)) {lfilename = argv[++i]; binaryfile = true;
      continue;} if((strcmp(argv[i],"-ufb")==0)) {ufilename = argv[++i];
      binaryfile = true; continue;} if((strcmp(argv[i],"--schedule")==0)) { i++;
        if((strcmp(argv[i],"auto")==0))
          schedule = AUTO;
        if((strcmp(argv[i],"dynamic")==0))
          schedule = DYNAMIC;
        if((strcmp(argv[i],"static")==0))
          schedule = STATIC;
        continue;
      }
    */
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help_sptrsv();
      return 0;
    }
  }

  if (tests.size() == 0) {
    tests.push_back(DEFAULT);
  }
  for (size_t i = 0; i < tests.size(); ++i) {
    std::cout << "tests[" << i << "] = " << tests[i] << std::endl;
  }

  Kokkos::initialize(argc, argv);
  {
    int total_errors = test_sptrsv_perf(tests, lfilename, ufilename, team_size, vector_length, idx_offset, loop,
                                        chain_threshold, dense_row_percent);

    if (total_errors == 0)
      printf("Kokkos::SPTRSV Test: Passed\n");
    else
      printf("Kokkos::SPTRSV Test: Failed\n");
  }
  Kokkos::finalize();
  return 0;
}
#else
int main() {
  std::cout << "KokkosSparse_sptrsv: This perf_test will do nothing when Cuda "
               "is enabled without lambda support."
            << std::endl;
  return 0;
}
#endif
