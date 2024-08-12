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

#ifndef KOKKOSLAPACK_CUSOLVER_HPP_
#define KOKKOSLAPACK_CUSOLVER_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#include <cusolverDn.h>

namespace KokkosLapack {
namespace Impl {

// Declaration of the singleton for cusolver
// this is the only header that needs to be
// included when using cusolverDn.
struct CudaLapackSingleton {
  cusolverDnHandle_t handle;

  CudaLapackSingleton();

  static CudaLapackSingleton& singleton();
};

inline void cusolver_internal_error_throw(cusolverStatus_t cusolverStatus, const char* name, const char* file,
                                          const int line) {
  std::ostringstream out;
  out << name << " error( ";
  switch (cusolverStatus) {
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      out << "CUSOLVER_STATUS_NOT_INITIALIZED): cusolver handle was not "
             "created correctly.";
      break;
    case CUSOLVER_STATUS_ALLOC_FAILED:
      out << "CUSOLVER_STATUS_ALLOC_FAILED): you might tried to allocate too "
             "much memory";
      break;
    case CUSOLVER_STATUS_INVALID_VALUE: out << "CUSOLVER_STATUS_INVALID_VALUE)"; break;
    case CUSOLVER_STATUS_ARCH_MISMATCH: out << "CUSOLVER_STATUS_ARCH_MISMATCH)"; break;
    case CUSOLVER_STATUS_EXECUTION_FAILED: out << "CUSOLVER_STATUS_EXECUTION_FAILED)"; break;
    case CUSOLVER_STATUS_INTERNAL_ERROR: out << "CUSOLVER_STATUS_INTERNAL_ERROR)"; break;
    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED: out << "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED)"; break;
    default: out << "unrecognized error code): this is bad!"; break;
  }
  if (file) {
    out << " " << file << ":" << line;
  }
  throw std::runtime_error(out.str());
}

inline void cusolver_internal_safe_call(cusolverStatus_t cusolverStatus, const char* name, const char* file = nullptr,
                                        const int line = 0) {
  if (CUSOLVER_STATUS_SUCCESS != cusolverStatus) {
    cusolver_internal_error_throw(cusolverStatus, name, file, line);
  }
}

// The macro below defines is the public interface for the safe cusolver calls.
// The functions themselves are protected by impl namespace.
#define KOKKOS_CUSOLVER_SAFE_CALL_IMPL(call) \
  KokkosLapack::Impl::cusolver_internal_safe_call(call, #call, __FILE__, __LINE__)

}  // namespace Impl
}  // namespace KokkosLapack
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSOLVER
#endif  // KOKKOSLAPACK_CUSOLVER_HPP_
