/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <iostream>
#include "KokkosKernels_config.h"
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_SparseUtils_cusparse.hpp"
#include "KokkosSparse_spadd.hpp"
#include "KokkosKernels_TestUtils.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include <cusparse.h>
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>
#include <mkl_spblas.h>

inline void spadd_mkl_internal_safe_call(sparse_status_t mklStatus,
                                         const char* name,
                                         const char* file = nullptr,
                                         const int line   = 0) {
  if (SPARSE_STATUS_SUCCESS != mklStatus) {
    std::ostringstream oss;
    oss << "MKL call \"" << name << "\" encountered error at " << file << ":"
        << line << '\n';
    Kokkos::abort(oss.str().c_str());
  }
}

#define SPADD_MKL_SAFE_CALL(call) \
  spadd_mkl_internal_safe_call(call, #call, __FILE__, __LINE__)
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) &&     \
    defined(KOKKOSKERNELS_INST_OFFSET_INT) && \
    defined(KOKKOSKERNELS_INST_ORDINAL_INT)

struct Params {
  int use_cuda     = 0;
  int use_openmp   = 0;
  int use_threads  = 0;
  int use_mkl      = 0;
  int use_cusparse = 0;
  bool sorted      = true;
  std::string amtx;
  std::string bmtx;
  std::string cmtx;
  int m         = 10000;
  int n         = 10000;
  int nnzPerRow = 30;
  bool bDiag = false;  // Whether B should be diagonal only (requires A square)
  bool verbose      = false;
  int repeat        = 1;
  int numericRepeat = 1;  // how many times to call numeric per overall run
};

template <typename crsMat_t>
void run_experiment(const Params& params) {
  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;

  using size_type  = typename crsMat_t::size_type;
  using lno_t      = typename crsMat_t::ordinal_type;
  using scalar_t   = typename crsMat_t::value_type;
  using device_t   = typename crsMat_t::device_type;
  using exec_space = typename device_t::execution_space;
  using mem_space  = typename device_t::memory_space;

  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
      size_type, lno_t, scalar_t, exec_space, mem_space, mem_space>;

  using graph_t   = typename crsMat_t::StaticCrsGraphType;
  using rowmap_t  = typename graph_t::row_map_type::non_const_type;
  using entries_t = typename graph_t::entries_type::non_const_type;
  using values_t  = typename crsMat_t::values_type::non_const_type;

  std::cout << "************************************* \n";
  std::cout << "************************************* \n";
  crsMat_t A;
  crsMat_t B;
  lno_t m = params.m;
  lno_t n = params.n;
  if (params.amtx.length()) {
    std::cout << "Loading A from " << params.amtx << '\n';
    A = KokkosKernels::Impl::read_kokkos_crst_matrix<crsMat_t>(
        params.amtx.c_str());
    m = A.numRows();
    n = A.numCols();
  } else {
    std::cout << "Randomly generating A\n";
    size_type nnzUnused = m * params.nnzPerRow;
    A = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
        m, n, nnzUnused, 0, (n + 3) / 3);
  }
  if (params.bmtx.length()) {
    std::cout << "Loading B from " << params.bmtx << '\n';
    B = KokkosKernels::Impl::read_kokkos_crst_matrix<crsMat_t>(
        params.bmtx.c_str());
  } else if (params.bDiag) {
    std::cout << "Generating B as diagonal matrix.\n";
    int diagLength = std::min(m, n);
    rowmap_t rowmap(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowmap_view"), m + 1);
    entries_t entries(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "colsmap_view"),
        diagLength);
    values_t values(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "values_view"),
        diagLength);
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
    B = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(
        m, n, nnzUnused, 0, (n + 3) / 3);
  }
  // Make sure dimensions are compatible
  if (A.numRows() != B.numRows() || A.numCols() != B.numCols()) {
    std::cout << "ERROR: A is " << A.numRows() << 'x' << A.numCols()
              << ", but B is " << B.numRows() << 'x' << B.numCols() << '\n';
    exit(1);
  }
  std::cout << "A and B are " << m << "x" << n << ". A, B have " << A.nnz()
            << " and " << B.nnz() << " entries.\n";

  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type
      lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type
      lno_nnz_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type const_lno_view_t;
  typedef
      typename crsMat_t::StaticCrsGraphType::entries_type const_lno_nnz_view_t;

  lno_view_t row_mapC;
  // entriesC, valuesC and cusparseBuffer are allocated inside
  // the loop, as part of symbolic
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  KernelHandle kh;

  if (params.sorted) {
    std::cout << "Assuming input matrices are sorted (explicitly sorting just "
                 "in case)\n";
    KokkosKernels::sort_crs_matrix(A);
    KokkosKernels::sort_crs_matrix(B);
  } else
    std::cout << "Assuming input matrices are not sorted.\n";
  kh.create_spadd_handle(params.sorted);
  auto addHandle = kh.get_spadd_handle();

  row_mapC = lno_view_t("non_const_lnow_row", m + 1);

  Kokkos::Timer timer;
  double symbolicTime = 0;
  double numericTime  = 0;

  // Do an untimed warm up symbolic, and preallocate space for C entries/values
  spadd_symbolic<KernelHandle, const_lno_view_t, const_lno_nnz_view_t,
                 const_lno_view_t, const_lno_nnz_view_t, lno_view_t,
                 lno_nnz_view_t>(&kh, A.graph.row_map, A.graph.entries,
                                 B.graph.row_map, B.graph.entries, row_mapC);

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
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetPointerMode(cusparseHandle, CUSPARSE_POINTER_MODE_HOST));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&A_cusparse));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&B_cusparse));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&C_cusparse));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatType(A_cusparse, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatType(B_cusparse, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatType(C_cusparse, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatDiagType(A_cusparse, CUSPARSE_DIAG_TYPE_NON_UNIT));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatDiagType(B_cusparse, CUSPARSE_DIAG_TYPE_NON_UNIT));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatDiagType(C_cusparse, CUSPARSE_DIAG_TYPE_NON_UNIT));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatIndexBase(A_cusparse, CUSPARSE_INDEX_BASE_ZERO));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatIndexBase(B_cusparse, CUSPARSE_INDEX_BASE_ZERO));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSetMatIndexBase(C_cusparse, CUSPARSE_INDEX_BASE_ZERO));
  }
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  sparse_matrix_t Amkl, Bmkl, Cmkl;
  if (params.use_mkl) {
    SPADD_MKL_SAFE_CALL(mkl_sparse_d_create_csr(
        &Amkl, SPARSE_INDEX_BASE_ZERO, m, n, (int*)A.graph.row_map.data(),
        (int*)A.graph.row_map.data() + 1, A.graph.entries.data(),
        A.values.data()));
    SPADD_MKL_SAFE_CALL(mkl_sparse_d_create_csr(
        &Bmkl, SPARSE_INDEX_BASE_ZERO, m, n, (int*)B.graph.row_map.data(),
        (int*)B.graph.row_map.data() + 1, B.graph.entries.data(),
        B.values.data()));
  }
#endif

  int c_nnz = 0;

  for (int sumRep = 0; sumRep < params.repeat; sumRep++) {
    timer.reset();
    if (use_kk) {
      spadd_symbolic<KernelHandle, const_lno_view_t, const_lno_nnz_view_t,
                     const_lno_view_t, const_lno_nnz_view_t, lno_view_t,
                     lno_nnz_view_t>(&kh, A.graph.row_map, A.graph.entries,
                                     B.graph.row_map, B.graph.entries,
                                     row_mapC);
      c_nnz = addHandle->get_c_nnz();
    } else if (params.use_cusparse) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      // Symbolic phase: compute buffer size, then compute nnz
      size_t bufferSize;
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDcsrgeam2_bufferSizeExt(
          cusparseHandle, A.numRows(), A.numCols(), &alphabeta, A_cusparse,
          A.nnz(), A.values.data(), A.graph.row_map.data(),
          A.graph.entries.data(), &alphabeta, B_cusparse, B.nnz(),
          B.values.data(), B.graph.row_map.data(), B.graph.entries.data(),
          C_cusparse, NULL, row_mapC.data(), NULL, &bufferSize));
      // Allocate work buffer
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          cudaMalloc((void**)&cusparseBuffer, bufferSize));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseXcsrgeam2Nnz(
          cusparseHandle, m, n, A_cusparse, A.nnz(), A.graph.row_map.data(),
          A.graph.entries.data(), B_cusparse, B.nnz(), B.graph.row_map.data(),
          B.graph.entries.data(), C_cusparse, row_mapC.data(), &c_nnz,
          cusparseBuffer));
#endif
    }
    if (!params.use_mkl) {
      entriesC = lno_nnz_view_t(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC (empty)"),
          c_nnz);
      valuesC = scalar_view_t(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC (empty)"),
          c_nnz);
    }

    // note: symbolic has a fence at the end
    symbolicTime += timer.seconds();
    timer.reset();
    // Just time all numeric repetitions together
    for (int numericRep = 0; numericRep < params.numericRepeat; numericRep++) {
      if (params.use_cusparse) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
        KOKKOS_CUSPARSE_SAFE_CALL(cusparseDcsrgeam2(
            cusparseHandle, m, n, &alphabeta, A_cusparse, A.nnz(),
            A.values.data(), A.graph.row_map.data(), A.graph.entries.data(),
            &alphabeta, B_cusparse, B.nnz(), B.values.data(),
            B.graph.row_map.data(), B.graph.entries.data(), C_cusparse,
            valuesC.data(), row_mapC.data(), entriesC.data(), cusparseBuffer));
#endif
      } else if (params.use_mkl) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
        SPADD_MKL_SAFE_CALL(mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE,
                                             Amkl, 1.0, Bmkl, &Cmkl));
        SPADD_MKL_SAFE_CALL(mkl_sparse_destroy(Cmkl));
#endif
      } else {
        spadd_numeric(
            &kh, A.graph.row_map, A.graph.entries, A.values, 1.0,  // A, alpha
            B.graph.row_map, B.graph.entries, B.values, 1.0,       // B, beta
            row_mapC, entriesC, valuesC);                          // C
      }
    }
    numericTime += timer.seconds();
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    if (params.use_cusparse)
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(cusparseBuffer));
#endif
  }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  if (params.use_cusparse)
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroy(cusparseHandle));
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  if (params.use_mkl) {
    SPADD_MKL_SAFE_CALL(mkl_sparse_destroy(Amkl));
    SPADD_MKL_SAFE_CALL(mkl_sparse_destroy(Bmkl));
  }
#endif

  int symbolicCalls = params.repeat;
  int numericCalls  = params.repeat * params.numericRepeat;

  std::cout << "Mean total time:    "
            << (symbolicTime / symbolicCalls) + (numericTime / numericCalls)
            << '\n'
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
    std::cout << "Writing C (" << m << "x" << n << ") to " << params.cmtx
              << "\n";
    crsMat_t C("C", m, n, c_nnz, valuesC, row_mapC, entriesC);
    KokkosKernels::Impl::write_kokkos_crst_matrix<crsMat_t>(
        C, params.cmtx.c_str());
  }
}

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr
      << "\t[Required] BACKEND: '--threads[numThreads]' | '--openmp "
         "[numThreads]' | '--cuda [cudaDeviceIndex]' | '--hip [hipDeviceIndex]'"
      << std::endl;

  std::cerr << "\t[Optional] --amtx <path> :: 1st input matrix" << std::endl;
  std::cerr << "\t[Optional] --bmtx <path> :: 2nd input matrix" << std::endl;
  std::cerr << "\t[Optional] --cmtx <path> :: output matrix for C = A+B"
            << std::endl;
  std::cerr << "\t[Optional] --mkl         :: run SpAdd from MKL" << std::endl;
  std::cerr << "\t[Optional] --cusparse    :: run SpAdd from cuSPARSE "
            << std::endl;
  std::cerr << "\t[Optional] --sorted      :: sort rows of inputs, and run the "
               "sorted algorithm"
            << std::endl;
  std::cerr << "\t[Optional] --unsorted    :: run the unsorted algorithm"
            << std::endl;
  std::cerr << "\t[Optional] --repeat      :: how many times to repeat overall "
               "spadd (symbolic + repeated numeric)"
            << std::endl;
  std::cerr << "\t[Optional] --numeric-repeat :: how many times to repeat "
               "numeric per symbolic"
            << std::endl;
  std::cerr << "\t[Optional] --verbose     :: enable verbose output"
            << std::endl;
  std::cerr << "\nSettings for randomly generated A/B matrices" << std::endl;
  std::cerr << "\t[Optional] --m           :: number of rows to generate"
            << std::endl;
  std::cerr << "\t[Optional] --n           :: number of cols to generate"
            << std::endl;
  std::cerr
      << "\t[Optional] --nnz         :: number of entries per row to generate"
      << std::endl;
  std::cerr
      << "\t[Optional] --nnz         :: number of entries per row to generate"
      << std::endl;
  std::cerr << "\t[Optional] --bdiag       :: generate B as a diagonal matrix"
            << std::endl;
}

int parse_inputs(Params& params, int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      params.use_threads = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      params.use_openmp = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      params.use_cuda = atoi(argv[++i]) + 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--mkl")) {
      params.use_mkl = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--cusparse")) {
      params.use_cusparse = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--sorted")) {
      params.sorted = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "--unsorted")) {
      params.sorted = false;
    } else if (0 == Test::string_compare_no_case(argv[i], "--amtx")) {
      // A at C=AxB
      params.amtx = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--bmtx")) {
      // B at C=AxB.
      // if not provided, C = AxA will be performed.
      params.bmtx = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--cmtx")) {
      // if provided, C will be written to given file.
      // has to have ".bin", or ".crs" extension.
      params.cmtx = argv[++i];
    } else if (0 == Test::string_compare_no_case(argv[i], "--m")) {
      params.m = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--n")) {
      params.n = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--nnz")) {
      params.nnzPerRow = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--bdiag")) {
      params.bDiag = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "--repeat")) {
      // if provided, C will be written to given file.
      // has to have ".bin", or ".crs" extension.
      params.repeat = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--numeric-repeat")) {
      // Reuse the symbolic step this many times.
      params.numericRepeat = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--verbose")) {
      params.verbose = true;
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": "
                << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}

int main(int argc, char** argv) {
  Params params;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }
  const int num_threads =
      params.use_openmp;  // Assumption is that use_openmp variable is provided
                          // as number of threads
  const int device_id = params.use_cuda - 1;

  Kokkos::initialize(Kokkos::InitArguments(num_threads, -1, device_id));
  // Kokkos::print_configuration(std::cout);

  // First, make sure that requested TPL (if any) is actually available
#if !defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
  if (params.use_mkl)
    throw std::invalid_argument(
        "To run MKL SpAdd, must enable the MKL TPL in cmake");
#endif
#if !defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
  if (params.use_cusparse)
    throw std::invalid_argument(
        "To run cuSPARSE SpAdd, must enable the cuSPARSE TPL in cmake");
#endif

  bool useOMP  = params.use_openmp != 0;
  bool useCUDA = params.use_cuda != 0;

  if (params.use_cusparse && !useCUDA) {
    throw std::invalid_argument(
        "To run cuSPARSE SpAdd, must supply the '--cuda <device id>' flag");
  }

  if (params.cmtx.length() && params.use_mkl) {
    throw std::invalid_argument(
        "If running MKL, can't output the result to file");
  }

  bool useSerial = !useOMP && !useCUDA;

  if (useOMP) {
#if defined(KOKKOS_ENABLE_OPENMP)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::OpenMP, void, int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    return 1;
#endif
  }
  if (useCUDA) {
#if defined(KOKKOS_ENABLE_CUDA)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::Cuda, void, int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: CUDA requested, but not available.\n";
    return 1;
#endif
  }
  if (useSerial) {
#if defined(KOKKOS_ENABLE_SERIAL)
    using crsMat_t =
        KokkosSparse::CrsMatrix<double, int, Kokkos::Serial, void, int>;
    run_experiment<crsMat_t>(params);
#else
    std::cout << "ERROR: Serial device requested, but not available.\n";
    return 1;
#endif
  }
  Kokkos::finalize();
  return 0;
}

#else
int main() {
#if !defined(KOKKOSKERNELS_INST_DOUBLE)
  std::cout << " not defined KOKKOSKERNELS_INST_DOUBLE" << std::endl;
#endif

#if !defined(KOKKOSKERNELS_INST_OFFSET_INT)
  std::cout << " not defined KOKKOSKERNELS_INST_OFFSET_INT" << std::endl;

#endif

#if !defined(KOKKOSKERNELS_INST_ORDINAL_INT)
  std::cout << " not defined KOKKOSKERNELS_INST_ORDINAL_INT" << std::endl;

#endif
}
#endif
