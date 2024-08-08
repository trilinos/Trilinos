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

#include <iostream>
#include "KokkosKernels_config.h"
#include "KokkosKernels_Handle.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_Utils_cusparse.hpp"
#include "KokkosSparse_Utils_mkl.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

#include "KokkosSparse_spadd.hpp"

using perf_test::CommonInputParams;

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include <cusparse.h>
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>
#include <mkl_spblas.h>
#endif

struct LocalParams {
  bool use_mkl      = false;
  bool use_cusparse = false;
  bool sorted       = true;
  std::string amtx;
  std::string bmtx;
  std::string cmtx;
  int m             = 10000;
  int n             = 10000;
  int nnzPerRow     = 30;
  bool bDiag        = false;  // Whether B should be diagonal only (requires A square)
  bool verbose      = false;
  int repeat        = 1;
  int numericRepeat = 1;  // how many times to call numeric per overall run
};

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr << perf_test::list_common_options();

  std::cerr << "\t[Optional] --amtx <path> :: 1st input matrix" << std::endl;
  std::cerr << "\t[Optional] --bmtx <path> :: 2nd input matrix" << std::endl;
  std::cerr << "\t[Optional] --cmtx <path> :: output matrix for C = A+B" << std::endl;
  std::cerr << "\t[Optional] --mkl         :: run SpAdd from MKL" << std::endl;
  std::cerr << "\t[Optional] --cusparse    :: run SpAdd from cuSPARSE " << std::endl;
  std::cerr << "\t[Optional] --sorted      :: sort rows of inputs, and run the "
               "sorted algorithm"
            << std::endl;
  std::cerr << "\t[Optional] --unsorted    :: run the unsorted algorithm" << std::endl;
  std::cerr << "\t[Optional] --repeat      :: how many times to repeat overall "
               "spadd (symbolic + repeated numeric)"
            << std::endl;
  std::cerr << "\t[Optional] --numeric-repeat :: how many times to repeat "
               "numeric per symbolic"
            << std::endl;
  std::cerr << "\t[Optional] --verbose     :: enable verbose output" << std::endl;
  std::cerr << "\nSettings for randomly generated A/B matrices" << std::endl;
  std::cerr << "\t[Optional] --m           :: number of rows to generate" << std::endl;
  std::cerr << "\t[Optional] --n           :: number of cols to generate" << std::endl;
  std::cerr << "\t[Optional] --nnz         :: number of entries per row to generate" << std::endl;
  std::cerr << "\t[Optional] --bdiag       :: generate B as a diagonal matrix" << std::endl;
}

int parse_inputs(LocalParams& params, int argc, char** argv) {
  bool discard;
  for (int i = 1; i < argc; ++i) {
    // if (perf_test::check_arg_str(i, argc, argv, "--amtx", params.amtx)) {
    //  ++i;
    if (perf_test::check_arg_bool(i, argc, argv, "--mkl", params.use_mkl)) {
    } else if (perf_test::check_arg_bool(i, argc, argv, "--cusparse", params.use_cusparse)) {
    } else if (perf_test::check_arg_bool(i, argc, argv, "--sorted", params.sorted)) {
    } else if (perf_test::check_arg_bool(i, argc, argv, "--unsorted", discard)) {
      params.sorted = false;
    } else if (perf_test::check_arg_str(i, argc, argv, "--amtx", params.amtx)) {
      // A at C=AxB
      ++i;
    } else if (perf_test::check_arg_str(i, argc, argv, "--bmtx", params.bmtx)) {
      // B at C=AxB.
      // if not provided, C = AxA will be performed.
      ++i;
    } else if (perf_test::check_arg_str(i, argc, argv, "--cmtx", params.cmtx)) {
      // if provided, C will be written to given file.
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--m", params.m)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--n", params.n)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--nnz", params.nnzPerRow)) {
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--bdiag", params.bDiag)) {
    } else if (perf_test::check_arg_int(i, argc, argv, "--repeat", params.repeat)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--numeric-repeat", params.numericRepeat)) {
      // Reuse the symbolic step this many times.
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--verbose", params.verbose)) {
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}

template <typename exec_space>
void run_experiment(int argc, char** argv, CommonInputParams) {
  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;

  using mem_space = typename exec_space::memory_space;
  using device_t  = typename Kokkos::Device<exec_space, mem_space>;
  using size_type = default_size_type;
  using lno_t     = default_lno_t;
  using scalar_t  = default_scalar;
  using crsMat_t  = KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type>;

  using KernelHandle =
      KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, exec_space, mem_space, mem_space>;

  using graph_t   = typename crsMat_t::StaticCrsGraphType;
  using rowmap_t  = typename graph_t::row_map_type::non_const_type;
  using entries_t = typename graph_t::entries_type::non_const_type;
  using values_t  = typename crsMat_t::values_type::non_const_type;

  LocalParams params;
  if (parse_inputs(params, argc, argv)) return;

    // First, make sure that requested TPL (if any) is actually available
#if !defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
  if (params.use_mkl) throw std::invalid_argument("To run MKL SpAdd, must enable the MKL TPL in cmake");
#endif
#if !defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
  if (params.use_cusparse) throw std::invalid_argument("To run cuSPARSE SpAdd, must enable the cuSPARSE TPL in cmake");
#else
  if (params.use_cusparse && !std::is_same<exec_space, Kokkos::Cuda>::value)
    throw std::invalid_argument("To run cuSPARSE SpAdd, must select the Cuda backend");
#endif

  if (params.cmtx.length() && params.use_mkl) {
    throw std::invalid_argument("If running MKL, can't output the result to file");
  }

  // Check that offset/ordinal types are compatible with any requested TPLs
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  if (params.use_mkl) {
    if constexpr (!std::is_same_v<int, MKL_INT>) {
      throw std::runtime_error("MKL configured with long long int not supported in Kokkos Kernels");
    }
    if constexpr (!std::is_same_v<MKL_INT, lno_t> || !std::is_same_v<MKL_INT, size_type>) {
      throw std::runtime_error(
          "Must enable int as both ordinal and offset type in KokkosKernels to "
          "call MKL SpAdd");
    }
  }
#endif

  if (params.use_cusparse) {
    if constexpr (!std::is_same_v<int, lno_t> || !std::is_same_v<int, size_type>) {
      throw std::runtime_error(
          "Must enable int as both ordinal and offset type in KokkosKernels to "
          "call cuSPARSE SpAdd");
    }
  }

  std::cout << "************************************* \n";
  crsMat_t A;
  crsMat_t B;
  lno_t m = params.m;
  lno_t n = params.n;
  if (params.amtx.length()) {
    std::cout << "Loading A from " << params.amtx << '\n';
    A = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(params.amtx.c_str());
    m = A.numRows();
    n = A.numCols();
  } else {
    std::cout << "Randomly generating A\n";
    size_type nnzUnused = m * params.nnzPerRow;
    A                   = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(m, n, nnzUnused, 0, (n + 3) / 3);
  }
  if (params.bmtx.length()) {
    std::cout << "Loading B from " << params.bmtx << '\n';
    B = KokkosSparse::Impl::read_kokkos_crst_matrix<crsMat_t>(params.bmtx.c_str());
  } else if (params.bDiag) {
    std::cout << "Generating B as diagonal matrix.\n";
    int diagLength = std::min(m, n);
    rowmap_t rowmap(Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowmap_view"), m + 1);
    entries_t entries(Kokkos::view_alloc(Kokkos::WithoutInitializing, "colsmap_view"), diagLength);
    values_t values(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values_view"), diagLength);
    auto rowmapHost  = Kokkos::create_mirror_view(rowmap);
    auto entriesHost = Kokkos::create_mirror_view(entries);
    auto valuesHost  = Kokkos::create_mirror_view(values);
    for (int i = 0; i < diagLength; i++) {
      rowmapHost(i)  = i;
      entriesHost(i) = i;
      valuesHost(i)  = 1.0;
    }
    for (int i = diagLength; i <= m; i++) {
      rowmapHost(i) = diagLength;
    }
    Kokkos::deep_copy(rowmap, rowmapHost);
    Kokkos::deep_copy(entries, entriesHost);
    Kokkos::deep_copy(values, valuesHost);
    B = crsMat_t("B", m, n, diagLength, values, rowmap, entries);
  } else {
    std::cout << "Randomly generating B\n";
    size_type nnzUnused = m * params.nnzPerRow;
    B                   = KokkosSparse::Impl::kk_generate_sparse_matrix<crsMat_t>(m, n, nnzUnused, 0, (n + 3) / 3);
  }
  // Make sure dimensions are compatible
  if (A.numRows() != B.numRows() || A.numCols() != B.numCols()) {
    std::cout << "ERROR: A is " << A.numRows() << 'x' << A.numCols() << ", but B is " << B.numRows() << 'x'
              << B.numCols() << '\n';
    exit(1);
  }
  std::cout << "A and B are " << m << "x" << n << ". A, B have " << A.nnz() << " and " << B.nnz() << " entries.\n";

  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;

  lno_view_t row_mapC;
  // entriesC, valuesC and cusparseBuffer are allocated inside
  // the loop, as part of symbolic
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  KernelHandle kh;

  if (params.sorted) {
    std::cout << "Assuming input matrices are sorted (explicitly sorting just "
                 "in case)\n";
    KokkosSparse::sort_crs_matrix(A);
    KokkosSparse::sort_crs_matrix(B);
  } else
    std::cout << "Assuming input matrices are not sorted.\n";
  kh.create_spadd_handle(params.sorted);
  auto addHandle = kh.get_spadd_handle();

  row_mapC = lno_view_t("non_const_lnow_row", m + 1);

  Kokkos::Timer timer;
  double symbolicTime = 0;
  double numericTime  = 0;

  // Do an untimed warm up symbolic, and preallocate space for C entries/values
  spadd_symbolic(exec_space{}, &kh, A.numRows(), A.numCols(), A.graph.row_map, A.graph.entries, B.graph.row_map,
                 B.graph.entries, row_mapC);

  bool use_kk = !params.use_cusparse && !params.use_mkl;

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  cusparseHandle_t cusparseHandle;
  cusparseMatDescr_t A_cusparse;
  cusparseMatDescr_t B_cusparse;
  cusparseMatDescr_t C_cusparse;
  char* cusparseBuffer;
  const double alphabeta = 1.0;

  if (params.use_cusparse) {
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreate(&cusparseHandle));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetPointerMode(cusparseHandle, CUSPARSE_POINTER_MODE_HOST));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&A_cusparse));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&B_cusparse));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&C_cusparse));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(A_cusparse, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(B_cusparse, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(C_cusparse, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(A_cusparse, CUSPARSE_DIAG_TYPE_NON_UNIT));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(B_cusparse, CUSPARSE_DIAG_TYPE_NON_UNIT));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(C_cusparse, CUSPARSE_DIAG_TYPE_NON_UNIT));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(A_cusparse, CUSPARSE_INDEX_BASE_ZERO));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(B_cusparse, CUSPARSE_INDEX_BASE_ZERO));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(C_cusparse, CUSPARSE_INDEX_BASE_ZERO));
  }
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  sparse_matrix_t Amkl = sparse_matrix_t(), Bmkl = sparse_matrix_t(), Cmkl = sparse_matrix_t();
  if (params.use_mkl) {
    if constexpr (std::is_same_v<lno_t, MKL_INT> && std::is_same_v<size_type, MKL_INT>) {
      KOKKOSKERNELS_MKL_SAFE_CALL(
          mkl_sparse_d_create_csr(&Amkl, SPARSE_INDEX_BASE_ZERO, m, n, (int*)A.graph.row_map.data(),
                                  (int*)A.graph.row_map.data() + 1, A.graph.entries.data(), A.values.data()));
      KOKKOSKERNELS_MKL_SAFE_CALL(
          mkl_sparse_d_create_csr(&Bmkl, SPARSE_INDEX_BASE_ZERO, m, n, (int*)B.graph.row_map.data(),
                                  (int*)B.graph.row_map.data() + 1, B.graph.entries.data(), B.values.data()));
    }
  }
#endif

  int c_nnz = 0;

  for (int sumRep = 0; sumRep < params.repeat; sumRep++) {
    timer.reset();
    if (use_kk) {
      spadd_symbolic(exec_space{}, &kh, A.numRows(), A.numCols(), A.graph.row_map, A.graph.entries, B.graph.row_map,
                     B.graph.entries, row_mapC);
      c_nnz = addHandle->get_c_nnz();
    } else if (params.use_cusparse) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      if constexpr (std::is_same_v<lno_t, int> && std::is_same_v<size_type, int>) {
        // Symbolic phase: compute buffer size, then compute nnz
        size_t bufferSize;
        KOKKOS_CUSPARSE_SAFE_CALL(cusparseDcsrgeam2_bufferSizeExt(
            cusparseHandle, A.numRows(), A.numCols(), &alphabeta, A_cusparse, A.nnz(), A.values.data(),
            A.graph.row_map.data(), A.graph.entries.data(), &alphabeta, B_cusparse, B.nnz(), B.values.data(),
            B.graph.row_map.data(), B.graph.entries.data(), C_cusparse, NULL, row_mapC.data(), NULL, &bufferSize));
        // Allocate work buffer
        KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void**)&cusparseBuffer, bufferSize));
        KOKKOS_CUSPARSE_SAFE_CALL(cusparseXcsrgeam2Nnz(cusparseHandle, m, n, A_cusparse, A.nnz(),
                                                       A.graph.row_map.data(), A.graph.entries.data(), B_cusparse,
                                                       B.nnz(), B.graph.row_map.data(), B.graph.entries.data(),
                                                       C_cusparse, row_mapC.data(), &c_nnz, cusparseBuffer));
      } else {
        throw std::runtime_error(
            "Must enable int as both ordinal and offset type in KokkosKernels "
            "to "
            "call cuSPARSE");
      }
#endif
    }
    if (!params.use_mkl) {
      entriesC = lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC (empty)"), c_nnz);
      valuesC  = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC (empty)"), c_nnz);
    }

    // note: symbolic has a fence at the end
    symbolicTime += timer.seconds();
    timer.reset();
    // Just time all numeric repetitions together
    for (int numericRep = 0; numericRep < params.numericRepeat; numericRep++) {
      if (params.use_cusparse) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
        if constexpr (std::is_same_v<lno_t, int> && std::is_same_v<size_type, int>) {
          KOKKOS_CUSPARSE_SAFE_CALL(cusparseDcsrgeam2(
              cusparseHandle, m, n, &alphabeta, A_cusparse, A.nnz(), A.values.data(), A.graph.row_map.data(),
              A.graph.entries.data(), &alphabeta, B_cusparse, B.nnz(), B.values.data(), B.graph.row_map.data(),
              B.graph.entries.data(), C_cusparse, valuesC.data(), row_mapC.data(), entriesC.data(), cusparseBuffer));
        }
#endif
      } else if (params.use_mkl) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
        if constexpr (std::is_same_v<lno_t, int> && std::is_same_v<size_type, int>) {
          KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, Amkl, 1.0, Bmkl, &Cmkl));
          KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_destroy(Cmkl));
        }
#endif
      } else {
        spadd_numeric(exec_space{}, &kh, A.numRows(), A.numCols(), A.graph.row_map, A.graph.entries, A.values,
                      1.0,  // A, alpha
                      B.graph.row_map, B.graph.entries, B.values,
                      1.0,                           // B, beta
                      row_mapC, entriesC, valuesC);  // C
      }
    }
    numericTime += timer.seconds();
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    if (params.use_cusparse) KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(cusparseBuffer));
#endif
  }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  if (params.use_cusparse) KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroy(cusparseHandle));
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  if (params.use_mkl) {
    KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_destroy(Amkl));
    KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_destroy(Bmkl));
  }
#endif

  int symbolicCalls = params.repeat;
  int numericCalls  = params.repeat * params.numericRepeat;

  std::cout << "Mean total time:    " << (symbolicTime / symbolicCalls) + (numericTime / numericCalls) << '\n'
            << "Mean symbolic time: " << (symbolicTime / symbolicCalls) << '\n'
            << "Mean numeric time:  " << (numericTime / numericCalls) << '\n';

  if (params.verbose) {
    std::cout << "row_mapC:" << row_mapC.extent(0) << std::endl;
    std::cout << "entriesC:" << entriesC.extent(0) << std::endl;
    std::cout << "valuesC:" << valuesC.extent(0) << std::endl;
    KokkosKernels::Impl::print_1Dview(valuesC);
    KokkosKernels::Impl::print_1Dview(entriesC);
    KokkosKernels::Impl::print_1Dview(row_mapC);
  }
  if (params.cmtx.length()) {
    std::cout << "Writing C (" << m << "x" << n << ") to " << params.cmtx << "\n";
    crsMat_t C("C", m, n, c_nnz, valuesC, row_mapC, entriesC);
    KokkosSparse::Impl::write_kokkos_crst_matrix<crsMat_t>(C, params.cmtx.c_str());
  }
}

#define KOKKOSKERNELS_PERF_TEST_NAME run_experiment
#include "KokkosKernels_perf_test_instantiation.hpp"
int main(int argc, char** argv) { return main_instantiation(argc, argv); }  // main
