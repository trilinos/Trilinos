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
#ifndef __TACHO_CUSOLVER_HPP__
#define __TACHO_CUSOLVER_HPP__

#if defined(KOKKOS_ENABLE_CUDA)
#include "Tacho_Util.hpp"
//#include "cusparse.h"
#include "cusolverSp.h"
#include "cusolverSp_LOWLEVEL_PREVIEW.h"
#include "cusparse_v2.h"

#include <iostream>
#include <type_traits>

namespace Tacho {

class CuSolver {
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
  cusolverSpHandle_t _handle;
  csrcholInfo_t _chol_info;
#if defined(TACHO_HAVE_CUSPARSE)
  cusparseMatDescr_t _desc;
#endif
  int _status;

  ordinal_type _m;
  size_type _nnz;

  size_type_array _ap; /// ap_ordinal is used to interface to cuSolver
  ordinal_type_array _ap_ordinal, _aj;

  value_type_array _buf;

  ordinal_type _verbose;

  void checkStatus(const char *s) {
    if (_status != 0)
      printf("Error: %s, status %d\n", s, _status);
  }

public:
  CuSolver() {
    _status = cusolverSpCreate(&_handle);
    checkStatus("cusolverSpCreate");
    _status = cusolverSpCreateCsrcholInfo(&_chol_info);
    checkStatus("cusolverSpCreateCsrcholInfo");
#if defined(TACHO_HAVE_CUSPARSE)
    _status = cusparseCreateMatDescr(&_desc);
    checkStatus("cusparseCreateMatDescr");
#else
    std::logic_error("CuSparse is not enabled");
#endif
  }
  virtual ~CuSolver() {
#if defined(TACHO_HAVE_CUSPARSE)
    _status = cusparseDestroyMatDescr(_desc);
    checkStatus("cusparseDestroyMatDescr");
#endif
    _status = cusolverSpDestroyCsrcholInfo(_chol_info);
    checkStatus("cusolverSpDestroyCsrcholInfo");
    _status = cusolverSpDestroy(_handle);
    checkStatus("cusolverSpDestroy");
  }

  void setVerbose(const ordinal_type verbose = 1) { _verbose = verbose; }

  template <typename arg_size_type_array, typename arg_ordinal_type_array>
  int analyze(const ordinal_type m, const arg_size_type_array &ap, const arg_ordinal_type_array &aj) {
    if (_verbose) {
      printf("cuSolver: Analyze\n");
      printf("=================\n");
    }
    Kokkos::Timer timer;

    _m = m;

    timer.reset();

    _ap = Kokkos::create_mirror_view(exec_memory_space(), ap);
    Kokkos::deep_copy(_ap, ap);
    _aj = Kokkos::create_mirror_view(exec_memory_space(), aj);
    Kokkos::deep_copy(_aj, aj);
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
    const double t_copy = timer.seconds();

    timer.reset();
#if defined(TACHO_HAVE_CUSPARSE)
    _status = cusolverSpXcsrcholAnalysis(_handle, _m, _nnz, _desc, _ap_ordinal.data(), _aj.data(), _chol_info);
#endif
    checkStatus("cusolverSpXcsrcholAnalysis");
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
      printf("\n");
    }

    return 0;
  }

  int factorize(const value_type_array &ax) {
    if (_verbose) {
      printf("cuSolver: Factorize\n");
      printf("===================\n");
    }
    Kokkos::Timer timer;

    timer.reset();
    size_t internalDataInBytes = 0, workspaceInBytes = 0;
#if defined(TACHO_HAVE_CUSPARSE)
    _status = cusolverSpDcsrcholBufferInfo(_handle, _m, _nnz, _desc, ax.data(), _ap_ordinal.data(), _aj.data(),
                                           _chol_info, &internalDataInBytes, &workspaceInBytes);
#endif
    checkStatus("cusolverSpDcsrcholBufferInfo");

    const size_t bufsize = workspaceInBytes / sizeof(value_type);
    if (bufsize > _buf.extent(0))
      _buf = value_type_array(do_not_initialize_tag("cusolver buf"), bufsize);
    Kokkos::fence();
    const double t_alloc = timer.seconds();

    timer.reset();
#if defined(TACHO_HAVE_CUSPARSE)
    _status = cusolverSpDcsrcholFactor(_handle, _m, _nnz, _desc, ax.data(), _ap_ordinal.data(), _aj.data(), _chol_info,
                                       _buf.data());
#endif
    checkStatus("cusolverSpDcsrcholFactor");
    Kokkos::fence();
    const double t_factor = timer.seconds();
    if (_verbose) {
      printf("  Time\n");
      printf("             time for workspace allocation:                   %10.6f s\n", t_alloc);
      printf("             time for numeric factorization:                  %10.6f s\n", t_factor);
      printf("             total time spent:                                %10.6f s\n", (t_alloc + t_factor));
      printf("\n");
      printf("  Workspace\n");
      printf("             internal data in MB:                          %10.3f MB\n",
             double(internalDataInBytes) / 1.e6);
      printf("             workspace in MB:                              %10.3f MB\n", double(workspaceInBytes) / 1.e6);
      printf("\n");
    }

    return 0;
  }

  int solve(const value_type_matrix &x, const value_type_matrix &b) {
    if (_verbose) {
      printf("cuSolver: Solve\n");
      printf("===============\n");
    }
    Kokkos::Timer timer;

    timer.reset();
    /// solve A x = t
    const ordinal_type len = x.extent(0), nrhs = x.extent(1);
    for (ordinal_type i = 0; i < nrhs; ++i) {
      _status = cusolverSpDcsrcholSolve(_handle, _m, b.data() + i * len, x.data() + i * len, _chol_info, _buf.data());
      checkStatus("cusolverSpDcsrcholSolve");
    }
    Kokkos::fence();
    const double t_solve = timer.seconds();
    if (_verbose) {
      printf("  Time\n");
      printf("             time for solve:                                  %10.6f s\n", t_solve);
      printf("\n");
    }

    return 0;
  }
};
} // namespace Tacho

#endif
#endif
