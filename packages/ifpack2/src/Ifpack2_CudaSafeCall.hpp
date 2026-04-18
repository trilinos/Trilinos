// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_CUDA_SAFE_CALL_HPP
#define IFPACK2_CUDA_SAFE_CALL_HPP

// Error handling below follows Kokkos::Impl::cuda_internal_safe_call and
// KOKKOS_IMPL_CUDA_SAFE_CALL; Ifpack2 defines its own helpers because
// KOKKOS_IMPL_* macros are not part of Kokkos' supported API.
// Source: https://github.com/kokkos/kokkos/blob/develop/core/src/Cuda/Kokkos_Cuda_Error.hpp

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_CUDA)

#include <Kokkos_Abort.hpp>

#include <cuda_runtime.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace Ifpack2 {
namespace Impl {

[[noreturn]] inline void cuda_internal_error_throw(cudaError e, const char *name,
                                                   const char *file, int line) {
  std::ostringstream out;
  out << name << " error( " << cudaGetErrorName(e) << "): " << cudaGetErrorString(e);
  if (file) {
    out << " " << file << ":" << line;
  }
  throw std::runtime_error(out.str());
}

#ifndef KOKKOS_COMPILER_NVHPC
[[noreturn]]
#endif
inline void
cuda_internal_error_abort(cudaError e, const char *name, const char *file,
                          int line) {
  std::ostringstream out;
  out << name << " error( " << cudaGetErrorName(e) << "): " << cudaGetErrorString(e);
  if (file) {
    out << " " << file << ":" << line;
  }
  Kokkos::abort(out.str().c_str());
}

inline void cuda_internal_safe_call(cudaError e, const char *name, const char *file,
                                    int line) {
  switch (e) {
    case cudaSuccess:
      break;
    case cudaErrorIllegalAddress:
    case cudaErrorAssert:
    case cudaErrorHardwareStackError:
    case cudaErrorIllegalInstruction:
    case cudaErrorMisalignedAddress:
    case cudaErrorInvalidAddressSpace:
    case cudaErrorInvalidPc:
    case cudaErrorLaunchFailure:
      cuda_internal_error_abort(e, name, file, line);
      break;
    default:
      cuda_internal_error_throw(e, name, file, line);
      break;
  }
}

}  // namespace Impl
}  // namespace Ifpack2

#define IFPACK2_IMPL_CUDA_SAFE_CALL(call) \
  Ifpack2::Impl::cuda_internal_safe_call(call, #call, __FILE__, __LINE__)

#endif  // KOKKOS_ENABLE_CUDA

#endif  // IFPACK2_CUDA_SAFE_CALL_HPP
