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

#ifndef KOKKOSKERNELS_ERROR_HPP
#define KOKKOSKERNELS_ERROR_HPP

#include <stdexcept>

namespace KokkosKernels {
namespace Impl {

inline void throw_runtime_exception(const std::string &msg) {
  throw std::runtime_error(msg);
}

#if defined(KOKKOS_ENABLE_HIP)
inline void hip_internal_error_throw(hipError_t e, const char *name,
                                     const char *file, const int line) {
  std::ostringstream out;
  out << name << " error( " << hipGetErrorName(e)
      << "): " << hipGetErrorString(e);
  if (file) {
    out << " " << file << ":" << line;
  }
  throw_runtime_exception(out.str());
}

inline void hip_internal_safe_call(hipError_t e, const char *name,
                                   const char *file = nullptr,
                                   const int line   = 0) {
  if (hipSuccess != e) {
    hip_internal_error_throw(e, name, file, line);
  }
}

#define KOKKOSKERNELS_IMPL_HIP_SAFE_CALL(call) \
  hip_internal_safe_call(call, #call, __FILE__, __LINE__)
#endif

}  // namespace Impl
}  // namespace KokkosKernels

#endif  // KOKKOSKERNELS_ERROR_HPP
