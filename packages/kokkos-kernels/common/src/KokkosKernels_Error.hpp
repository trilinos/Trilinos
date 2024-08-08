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
#include <sstream>

namespace KokkosKernels {
namespace Impl {

inline void throw_runtime_exception(const std::string &msg) { throw std::runtime_error(msg); }

#if defined(KOKKOS_ENABLE_HIP)
inline void hip_internal_error_throw(hipError_t e, const char *name, const char *file, const int line) {
  std::ostringstream out;
  out << name << " error( " << hipGetErrorName(e) << "): " << hipGetErrorString(e);
  if (file) {
    out << " " << file << ":" << line;
  }
  throw_runtime_exception(out.str());
}

inline void hip_internal_safe_call(hipError_t e, const char *name, const char *file = nullptr, const int line = 0) {
  if (hipSuccess != e) {
    hip_internal_error_throw(e, name, file, line);
  }
}

#define KOKKOSKERNELS_IMPL_HIP_SAFE_CALL(call) hip_internal_safe_call(call, #call, __FILE__, __LINE__)
#endif

}  // namespace Impl
}  // namespace KokkosKernels

/*
 * Asserts and error checking macros/functions.
 *
 * KK_KERNEL** are for error checking within kokkos kernels.
 *
 * Any check with "assert" in the name is disabled for release builds
 *
 * For _MSG checks, the msg argument can contain '<<' if not a kernel check.
 *
 * KK_USER_REQUIRE* are for checking user inputs
 *
 * This code is adapted from EKAT/src/ekat/ekat_assert.hpp
 */

// Internal do not call directly
#define IMPL_THROW(condition, msg, exception_type)                      \
  do {                                                                  \
    if (!(condition)) {                                                 \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition; \
      _ss_ << "\n" << msg;                                              \
      throw exception_type(_ss_.str());                                 \
    }                                                                   \
  } while (0)

// SYCL cannot printf like the other backends quite yet
#define IMPL_KERNEL_THROW(condition, msg)                                      \
  do {                                                                         \
    if (!(condition)) {                                                        \
      Kokkos::printf("KERNEL CHECK FAILED:\n   %s\n   %s\n", #condition, msg); \
      Kokkos::abort("");                                                       \
    }                                                                          \
  } while (0)

#ifndef NDEBUG
#define KK_ASSERT(condition) IMPL_THROW(condition, "", std::logic_error)
#define KK_ASSERT_MSG(condition, msg) IMPL_THROW(condition, msg, std::logic_error)
#define KK_KERNEL_ASSERT(condition) IMPL_KERNEL_THROW(condition, "")
#define KK_KERNEL_ASSERT_MSG(condition, msg) IMPL_KERNEL_THROW(condition, msg)
#else
#define KK_ASSERT(condition) ((void)(0))
#define KK_ASSERT_MSG(condition, msg) ((void)(0))
#define KK_KERNEL_ASSERT(condition) ((void)(0))
#define KK_KERNEL_ASSERT_MSG(condition, msg) ((void)(0))
#endif

#define KK_REQUIRE(condition) IMPL_THROW(condition, "", std::logic_error)
#define KK_REQUIRE_MSG(condition, msg) IMPL_THROW(condition, msg, std::logic_error)

#define KK_USER_REQUIRE(condition) IMPL_THROW(condition, "", std::runtime_error)
#define KK_USER_REQUIRE_MSG(condition, msg) IMPL_THROW(condition, msg, std::runtime_error)

#define KK_KERNEL_REQUIRE(condition) IMPL_KERNEL_THROW(condition, "")
#define KK_KERNEL_REQUIRE_MSG(condition, msg) IMPL_KERNEL_THROW(condition, msg)

#define KK_ERROR_MSG(msg) KK_REQUIRE_MSG(false, msg)
#define KK_KERNEL_ERROR_MSG(msg) KK_KERNEL_REQUIRE_MSG(false, msg)

#endif  // KOKKOSKERNELS_ERROR_HPP
