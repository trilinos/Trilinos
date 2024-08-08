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

#ifndef KOKKOSBLAS_TPL_SPEC_HPP_
#define KOKKOSBLAS_TPL_SPEC_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include "cuda_runtime.h"
#include "cublas_v2.h"

namespace KokkosBlas {
namespace Impl {

struct CudaBlasSingleton {
  cublasHandle_t handle;

  CudaBlasSingleton();

  static CudaBlasSingleton& singleton();
};

inline void cublas_internal_error_throw(cublasStatus_t cublasState, const char* name, const char* file,
                                        const int line) {
  std::ostringstream out;
  // out << name << " error( " << cublasGetStatusName(cublasState)
  //     << "): " << cublasGetStatusString(cublasState);
  out << name << " error( ";
  switch (cublasState) {
    case CUBLAS_STATUS_NOT_INITIALIZED:
      out << "CUBLAS_STATUS_NOT_INITIALIZED): the library was not initialized.";
      break;
    case CUBLAS_STATUS_ALLOC_FAILED: out << "CUBLAS_STATUS_ALLOC_FAILED): the resource allocation failed."; break;
    case CUBLAS_STATUS_INVALID_VALUE:
      out << "CUBLAS_STATUS_INVALID_VALUE): an invalid numerical value was "
             "used as an argument.";
      break;
    case CUBLAS_STATUS_ARCH_MISMATCH:
      out << "CUBLAS_STATUS_ARCH_MISMATCH): an absent device architectural "
             "feature is required.";
      break;
    case CUBLAS_STATUS_MAPPING_ERROR:
      out << "CUBLAS_STATUS_MAPPING_ERROR): an access to GPU memory space "
             "failed.";
      break;
    case CUBLAS_STATUS_EXECUTION_FAILED:
      out << "CUBLAS_STATUS_EXECUTION_FAILED): the GPU program failed to "
             "execute.";
      break;
    case CUBLAS_STATUS_INTERNAL_ERROR: out << "CUBLAS_STATUS_INTERNAL_ERROR): an internal operation failed."; break;
    case CUBLAS_STATUS_NOT_SUPPORTED:
      out << "CUBLAS_STATUS_NOT_SUPPORTED): the feature required is not "
             "supported.";
      break;
    default: out << "unrecognized error code): this is bad!"; break;
  }
  if (file) {
    out << " " << file << ":" << line;
  }
  throw std::runtime_error(out.str());
}

inline void cublas_internal_safe_call(cublasStatus_t cublasState, const char* name, const char* file = nullptr,
                                      const int line = 0) {
  if (CUBLAS_STATUS_SUCCESS != cublasState) {
    cublas_internal_error_throw(cublasState, name, file, line);
  }
}

// The macro below defines the interface for the safe cublas calls.
// The functions themselves are protected by impl namespace and this
// is not meant to be used by external application or libraries.
#define KOKKOS_CUBLAS_SAFE_CALL_IMPL(call) KokkosBlas::Impl::cublas_internal_safe_call(call, #call, __FILE__, __LINE__)

/// \brief This function converts KK transpose mode to cuBLAS transpose mode
inline cublasOperation_t trans_mode_kk_to_cublas(const char kkMode[]) {
  cublasOperation_t trans;
  if ((kkMode[0] == 'N') || (kkMode[0] == 'n'))
    trans = CUBLAS_OP_N;
  else if ((kkMode[0] == 'T') || (kkMode[0] == 't'))
    trans = CUBLAS_OP_T;
  else
    trans = CUBLAS_OP_C;
  return trans;
}

}  // namespace Impl
}  // namespace KokkosBlas
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCBLAS
#include <rocblas/rocblas.h>

namespace KokkosBlas {
namespace Impl {

struct RocBlasSingleton {
  rocblas_handle handle;

  RocBlasSingleton();

  static RocBlasSingleton& singleton();
};

inline void rocblas_internal_error_throw(rocblas_status rocblasState, const char* name, const char* file,
                                         const int line) {
  std::ostringstream out;
  out << name << " error( ";
  switch (rocblasState) {
    case rocblas_status_invalid_handle:
      out << "rocblas_status_invalid_handle): handle not initialized, invalid "
             "or null.";
      break;
    case rocblas_status_not_implemented: out << "rocblas_status_not_implemented): function is not implemented."; break;
    case rocblas_status_invalid_pointer: out << "rocblas_status_invalid_pointer): invalid pointer argument."; break;
    case rocblas_status_invalid_size: out << "rocblas_status_invalid_size): invalid size argument."; break;
    case rocblas_status_memory_error:
      out << "rocblas_status_memory_error): failed internal memory allocation, "
             "copy or dealloc.";
      break;
    case rocblas_status_internal_error: out << "rocblas_status_internal_error): other internal library failure."; break;
    case rocblas_status_perf_degraded:
      out << "rocblas_status_perf_degraded): performance degraded due to low "
             "device memory.";
      break;
    case rocblas_status_size_query_mismatch: out << "unmatched start/stop size query): ."; break;
    case rocblas_status_size_increased:
      out << "rocblas_status_size_increased): queried device memory size "
             "increased.";
      break;
    case rocblas_status_size_unchanged:
      out << "rocblas_status_size_unchanged): queried device memory size "
             "unchanged.";
      break;
    case rocblas_status_invalid_value: out << "rocblas_status_invalid_value): passed argument not valid."; break;
    case rocblas_status_continue:
      out << "rocblas_status_continue): nothing preventing function to "
             "proceed.";
      break;
    case rocblas_status_check_numerics_fail:
      out << "rocblas_status_check_numerics_fail): will be set if the "
             "vector/matrix has a NaN or an Infinity.";
      break;
    default: out << "unrecognized error code): this is bad!"; break;
  }
  if (file) {
    out << " " << file << ":" << line;
  }
  throw std::runtime_error(out.str());
}

inline void rocblas_internal_safe_call(rocblas_status rocblasState, const char* name, const char* file = nullptr,
                                       const int line = 0) {
  if (rocblas_status_success != rocblasState) {
    rocblas_internal_error_throw(rocblasState, name, file, line);
  }
}

// The macro below defines the interface for the safe rocblas calls.
// The functions themselves are protected by impl namespace and this
// is not meant to be used by external application or libraries.
#define KOKKOS_ROCBLAS_SAFE_CALL_IMPL(call) \
  KokkosBlas::Impl::rocblas_internal_safe_call(call, #call, __FILE__, __LINE__)

/// \brief This function converts KK transpose mode to rocBLAS transpose mode
inline rocblas_operation trans_mode_kk_to_rocblas(const char kkMode[]) {
  rocblas_operation trans;
  if ((kkMode[0] == 'N') || (kkMode[0] == 'n'))
    trans = rocblas_operation_none;
  else if ((kkMode[0] == 'T') || (kkMode[0] == 't'))
    trans = rocblas_operation_transpose;
  else
    trans = rocblas_operation_conjugate_transpose;
  return trans;
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCBLAS

#endif  // KOKKOSBLAS_TPL_SPEC_HPP_
