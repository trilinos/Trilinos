// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_CUSPARSE_TRISOLVE_HPP__
#define __TACHO_CUSPARSE_TRISOLVE_HPP__

#if defined(KOKKOS_ENABLE_CUDA)
#include "Tacho_Util.hpp"
#include "cusparse_v2.h"

#include <iostream>
#include <type_traits>

namespace Tacho {

class CuSparseTriSolve {
public:
  typedef double value_type;

  typedef typename UseThisDevice<Kokkos::Cuda>::type device_type;

  typedef typename device_type::execution_space exec_space;
  typedef typename device_type::memory_space exec_memory_space;

  typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;
  typedef typename host_device_type::execution_space host_space;
  typedef typename host_device_type::memory_space host_memory_space;

  typedef Kokkos::View<size_type *, device_type> size_type_array;
  typedef Kokkos::View<ordinal_type *, device_type> ordinal_type_array;
  typedef Kokkos::View<value_type *, device_type> value_type_array;
  typedef Kokkos::View<value_type **, Kokkos::LayoutLeft, device_type> value_type_matrix;

  typedef Kokkos::View<size_type *, host_device_type> size_type_array_host;
  typedef Kokkos::View<ordinal_type *, host_device_type> ordinal_type_array_host;
  typedef Kokkos::View<value_type *, host_device_type> value_type_array_host;
  typedef Kokkos::View<value_type **, Kokkos::LayoutLeft, host_device_type> value_type_matrix_host;

private:
  /// pair of cuSparse objects for Lower and Upper
  std::pair<cusparseHandle_t, cusparseHandle_t> _handle;
  std::pair<cusparseMatDescr_t, cusparseMatDescr_t> _desc;
  std::pair<csrsv2Info_t, csrsv2Info_t> _info;
  std::pair<cusparseSolvePolicy_t, cusparseSolvePolicy_t> _policy;
  std::pair<cusparseOperation_t, cusparseOperation_t> _trans;

  cudaStream_t _stream;
  cudaGraph_t _graph;
  cudaGraphExec_t _instance;
  int _graph_status;
  int _status;

  ordinal_type _m;
  size_type _nnz;

  size_type_array _ap; /// ap_ordinal is used to interface to cuSparse
  ordinal_type_array _ap_ordinal, _aj;

  value_type_array _ax, _buf;

  ordinal_type _verbose;

  void checkStatus(const char *s) {
    if (_status != 0) {
      printf("Error: %s, status %d\n", s, _status);
      throw std::runtime_error("Error: checkStatus returns a non-zero status value");
    }
  }

public:
  CuSparseTriSolve() {
    /// let's use a default stream to interoperate with kokkos without execution space
    /// later we can create a separate stream but for now it is not worth
    _stream = NULL;
    ///_status = cudaStreamCreate(&_stream); checkStatus("cudaStreamCreate");

    _status = cusparseCreate(&_handle.first);
    checkStatus("cusparseCreate::Lower");
    _status = cusparseCreate(&_handle.second);
    checkStatus("cusparseCreate::Upper");
    _status = cusparseCreateMatDescr(&_desc.first);
    checkStatus("cusparseCreateMatDescr::Lower");
    _status = cusparseCreateMatDescr(&_desc.second);
    checkStatus("cusparseCreateMatDescr::Upper");
    _status = cusparseCreateCsrsv2Info(&_info.first);
    checkStatus("cusparseCreateCsrsv2Info::Lower");
    _status = cusparseCreateCsrsv2Info(&_info.second);
    checkStatus("cusparseCreateCsrsv2Info::Upper");

    _status = cusparseSetStream(_handle.first, _stream);
    _status = cusparseSetStream(_handle.second, _stream);

    _graph_status = 0;
  }
  virtual ~CuSparseTriSolve() {
    _status = cusparseDestroyCsrsv2Info(_info.second);
    checkStatus("cusparseCreateCsrsv2Info::Upper");
    _status = cusparseDestroyCsrsv2Info(_info.first);
    checkStatus("cusparseCreateCsrsv2Info::Lower");
    _status = cusparseDestroyMatDescr(_desc.second);
    checkStatus("cusparseDestroyMatDescr::Upper");
    _status = cusparseDestroyMatDescr(_desc.first);
    checkStatus("cusparseDestroyMatDescr::Lower");
    _status = cusparseDestroy(_handle.second);
    checkStatus("cusparseDestroy::Upper");
    _status = cusparseDestroy(_handle.first);
    checkStatus("cusparseDestroy::Lower");
    ///_status = cudaStreamDestroy(_stream);
  }

  void setVerbose(const ordinal_type verbose = 1) { _verbose = verbose; }

  void getStream(cudaStream_t *stream) { *stream = _stream; }

  template <typename arg_size_type_array, typename arg_ordinal_type_array, typename arg_value_type_array>
  int analyze(const ordinal_type m, const arg_size_type_array &ap, const arg_ordinal_type_array &aj,
              const arg_value_type_array &ax) {
    if (_verbose) {
      printf("cuSparse: Analyze\n");
      printf("=================\n");
    }
    Kokkos::Timer timer;

    _m = m;

    timer.reset();

    _ap = Kokkos::create_mirror_view(exec_memory_space(), ap);
    Kokkos::deep_copy(_ap, ap);
    _aj = Kokkos::create_mirror_view(exec_memory_space(), aj);
    Kokkos::deep_copy(_aj, aj);
    _ax = Kokkos::create_mirror_view(exec_memory_space(), ax);
    Kokkos::deep_copy(_ax, ax);
    Kokkos::fence();
#if defined(TACHO_USE_INT_INT)
    _ap_ordinal = ap;
#else
    /// LAMBDA cannot capture this pointer; make all variables local
    {
      ordinal_type_array l_ap_ordinal(do_not_initialize_tag("CuSolver::ap_ordinal"), _ap.extent(0));
      auto l_ap = _ap;
      Kokkos::RangePolicy<exec_space, Kokkos::Schedule<Kokkos::Static>> policy(0, l_ap.extent(0));
      Kokkos::parallel_for(
          policy, KOKKOS_LAMBDA(const int i) { l_ap_ordinal(i) = static_cast<int>(l_ap(i)); });
      _ap_ordinal = l_ap_ordinal;
    }
#endif
    auto last = Kokkos::subview(_ap, _m);
    auto h_last = Kokkos::create_mirror_view(host_memory_space(), last);
    Kokkos::deep_copy(h_last, last);
    _nnz = h_last();

    Kokkos::fence();
    const double t_copy = timer.seconds();

    timer.reset();
    std::pair<int, int> bufSizeInBytes;

    _policy.first = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
    _policy.second = CUSPARSE_SOLVE_POLICY_USE_LEVEL;

    _trans.first = CUSPARSE_OPERATION_TRANSPOSE;
    _trans.second = CUSPARSE_OPERATION_NON_TRANSPOSE;

    ///
    /// query buffersize and analyze the matrix; two separate handles are required for book-keeping
    ///
    _status = cusparseSetMatIndexBase(_desc.first, CUSPARSE_INDEX_BASE_ZERO);
    checkStatus("cusparseSetMatIndexBase::Lower");
    _status = cusparseSetMatFillMode(_desc.first, CUSPARSE_FILL_MODE_UPPER);
    checkStatus("cusparseSetMatFillMode::Lower");
    _status = cusparseSetMatDiagType(_desc.first, CUSPARSE_DIAG_TYPE_NON_UNIT);
    checkStatus("cusparseSetMatDiagType::Lower");
    _status = cusparseDcsrsv2_bufferSize(_handle.first, _trans.first, _m, _nnz, _desc.first, _ax.data(),
                                         _ap_ordinal.data(), _aj.data(), _info.first, &bufSizeInBytes.first);
    checkStatus("cusparseDcsrsv2_bufferSize::Lower");

    _status = cusparseSetMatIndexBase(_desc.second, CUSPARSE_INDEX_BASE_ZERO);
    checkStatus("cusparseSetMatIndexBase::Upper");
    _status = cusparseSetMatFillMode(_desc.second, CUSPARSE_FILL_MODE_UPPER);
    checkStatus("cusparseSetMatFillMode::Upper");
    _status = cusparseSetMatDiagType(_desc.second, CUSPARSE_DIAG_TYPE_NON_UNIT);
    checkStatus("cusparseSetMatDiagType::Upper");
    _status = cusparseDcsrsv2_bufferSize(_handle.second, _trans.second, _m, _nnz, _desc.second, _ax.data(),
                                         _ap_ordinal.data(), _aj.data(), _info.second, &bufSizeInBytes.second);
    checkStatus("cusparseDcsrsv2_bufferSize::Upper");

    const int maxBufSizeInBytes = std::max(bufSizeInBytes.first, bufSizeInBytes.second);
    const int bufsize = maxBufSizeInBytes > 0 ? maxBufSizeInBytes / sizeof(double) : 32;

    if (bufsize > int(_buf.extent(0)))
      _buf = value_type_array(do_not_initialize_tag("buf"), bufsize);

    _status = cusparseDcsrsv2_analysis(_handle.first, _trans.first, _m, _nnz, _desc.first, _ax.data(),
                                       _ap_ordinal.data(), _aj.data(), _info.first, _policy.first, _buf.data());
    checkStatus("cusparseDcsrsv2_analysis::Lower");

    _status = cusparseDcsrsv2_analysis(_handle.second, _trans.second, _m, _nnz, _desc.second, _ax.data(),
                                       _ap_ordinal.data(), _aj.data(), _info.second, _policy.second, _buf.data());
    checkStatus("cusparseDcsrsv2_analysis::Upper");
    Kokkos::fence();
    const double t_analyze = timer.seconds();

    if (_verbose) {
      printf("  Linear system A\n");
      printf("             number of equations:                             %10d\n", _m);
      printf("             number of nonzeros:                              %10d\n", _nnz);
      printf("\n");
      printf("  Time\n");
      printf("             time for copying A into U:                       %10.6f s\n", t_copy);
      printf("             time for analysis:                               %10.6f s\n", t_analyze);
      printf("             total time spent:                                %10.6f s\n", (t_copy + t_analyze));
      printf("  Workspace\n");
      printf("             upper solve workspace in MB:                  %10.3f MB\n",
             double(bufSizeInBytes.second) / 1.e6);
      printf("             lower solve workspace in MB:                  %10.3f MB\n",
             double(bufSizeInBytes.first) / 1.e6);
      printf("             max workspace in MB:                          %10.3f MB\n",
             double(maxBufSizeInBytes) / 1.e6);
      printf("\n");
    }

    return 0;
  }

  int solve(const value_type_matrix &x, const value_type_matrix &b, const value_type_matrix &t,
            const bool apply_fence = true) {
    if (_verbose) {
      printf("cuSolver: Solve\n");
      printf("===============\n");
    }
    Kokkos::Timer timer;

    timer.reset();

    /// solve A x = t
    const value_type one(1);
    const ordinal_type len = x.extent(0), nrhs = x.extent(1);
    for (ordinal_type i = 0; i < nrhs; ++i) {
      _status = cusparseDcsrsv2_solve(_handle.first, _trans.first, _m, _nnz, (const double *)&one, _desc.first,
                                      (const double *)_ax.data(), (const int *)_ap_ordinal.data(),
                                      (const int *)_aj.data(), _info.first, (const double *)b.data() + i * len,
                                      (double *)t.data() + i * len, _policy.first, (void *)_buf.data());
      checkStatus("cusparseDcsrsv2_solve::Lower");

      _status = cusparseDcsrsv2_solve(_handle.second, _trans.second, _m, _nnz, (const double *)&one, _desc.second,
                                      (const double *)_ax.data(), (const int *)_ap_ordinal.data(),
                                      (const int *)_aj.data(), _info.second, (const double *)t.data() + i * len,
                                      (double *)x.data() + i * len, _policy.second, (void *)_buf.data());
      checkStatus("cusparseDcsrsv2_solve::Upper");
    }
    if (apply_fence)
      Kokkos::fence();

    const double t_solve = timer.seconds();
    if (_verbose) {
      printf("  Time\n");
      printf("             time for solve:                                  %10.6f s\n", t_solve);
      printf("\n");
    }

    return 0;
  }

  int solve_capture(const value_type_matrix &x, const value_type_matrix &b, const value_type_matrix &t) {
    _graph_status = 1; /// capture begin
    cudaStreamBeginCapture(_stream, cudaStreamCaptureModeGlobal);
    solve(x, b, t, false); /// do not apply fence when it capture
    cudaStreamEndCapture(_stream, &_graph);
    cudaGraphInstantiate(&_instance, _graph, NULL, NULL, 0);

    return 0;
  }

  int solve_launch() {
    cudaGraphLaunch(_instance, _stream);
    cudaStreamSynchronize(_stream);
    return 0;
  }
};
} // namespace Tacho

#endif
#endif
