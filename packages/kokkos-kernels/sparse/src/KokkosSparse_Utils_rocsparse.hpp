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

#ifndef _KOKKOSKERNELS_SPARSEUTILS_ROCSPARSE_HPP
#define _KOKKOSKERNELS_SPARSEUTILS_ROCSPARSE_HPP

#include <type_traits>
#include <sstream>

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#if __has_include(<rocm-core/rocm_version.h>)
#include <rocm-core/rocm_version.h>
#else
#include <rocm_version.h>
#endif
#include <rocsparse/rocsparse.h>

namespace KokkosSparse {
namespace Impl {

inline void rocsparse_internal_error_throw(rocsparse_status rocsparseStatus, const char* name, const char* file,
                                           const int line) {
  std::ostringstream out;
  out << name << " error( ";
  switch (rocsparseStatus) {
    case rocsparse_status_invalid_handle:
      out << "rocsparse_status_invalid_handle): handle not initialized, "
             "invalid or null.";
      break;
    case rocsparse_status_not_implemented:
      out << "rocsparse_status_not_implemented): function is not implemented.";
      break;
    case rocsparse_status_invalid_pointer:
      out << "rocsparse_status_invalid_pointer): invalid pointer parameter.";
      break;
    case rocsparse_status_invalid_size: out << "rocsparse_status_invalid_size): invalid size parameter."; break;
    case rocsparse_status_memory_error:
      out << "rocsparse_status_memory_error): failed memory allocation, copy, "
             "dealloc.";
      break;
    case rocsparse_status_internal_error:
      out << "rocsparse_status_internal_error): other internal library "
             "failure.";
      break;
    case rocsparse_status_invalid_value: out << "rocsparse_status_invalid_value): invalid value parameter."; break;
    case rocsparse_status_arch_mismatch: out << "rocsparse_status_arch_mismatch): device arch is not supported."; break;
    case rocsparse_status_zero_pivot: out << "rocsparse_status_zero_pivot): encountered zero pivot."; break;
    case rocsparse_status_not_initialized:
      out << "rocsparse_status_not_initialized): descriptor has not been "
             "initialized.";
      break;
    case rocsparse_status_type_mismatch: out << "rocsparse_status_type_mismatch): index types do not match."; break;
    default: out << "unrecognized error code): this is bad!"; break;
  }
  if (file) {
    out << " " << file << ":" << line;
  }
  throw std::runtime_error(out.str());
}

inline void rocsparse_internal_safe_call(rocsparse_status rocsparseStatus, const char* name, const char* file = nullptr,
                                         const int line = 0) {
  if (rocsparse_status_success != rocsparseStatus) {
    rocsparse_internal_error_throw(rocsparseStatus, name, file, line);
  }
}

// The macro below defines is the public interface for the safe cusparse calls.
// The functions themselves are protected by impl namespace.
#define KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(call) \
  KokkosSparse::Impl::rocsparse_internal_safe_call(call, #call, __FILE__, __LINE__)

inline rocsparse_operation mode_kk_to_rocsparse(const char kk_mode[]) {
  rocsparse_operation myRocsparseOperation;
  switch (toupper(kk_mode[0])) {
    case 'N': myRocsparseOperation = rocsparse_operation_none; break;
    case 'T': myRocsparseOperation = rocsparse_operation_transpose; break;
    case 'H': myRocsparseOperation = rocsparse_operation_conjugate_transpose; break;
    default: {
      std::ostringstream out;
      out << "Mode " << kk_mode[0] << " invalid for rocSPARSE SpMV.\n";
      throw std::invalid_argument(out.str());
    }
  }
  return myRocsparseOperation;
}

template <typename index_type>
inline rocsparse_indextype rocsparse_index_type() {
  if (std::is_same<index_type, uint16_t>::value) {
    return rocsparse_indextype_u16;
  } else if (std::is_same<index_type, int32_t>::value) {
    return rocsparse_indextype_i32;
  } else if (std::is_same<index_type, int64_t>::value) {
    return rocsparse_indextype_i64;
  } else {
    std::ostringstream out;
    out << "Trying to call rocSPARSE SpMV with unsupported index type: " << typeid(index_type).name();
    throw std::logic_error(out.str());
  }
}

template <typename data_type>
inline rocsparse_datatype rocsparse_compute_type() {
  std::ostringstream out;
  out << "Trying to call rocSPARSE SpMV with unsupported compute type: " << typeid(data_type).name();
  throw std::logic_error(out.str());
}

template <>
inline rocsparse_datatype rocsparse_compute_type<float>() {
  return rocsparse_datatype_f32_r;
}

template <>
inline rocsparse_datatype rocsparse_compute_type<double>() {
  return rocsparse_datatype_f64_r;
}

template <>
inline rocsparse_datatype rocsparse_compute_type<Kokkos::complex<float>>() {
  return rocsparse_datatype_f32_c;
}

template <>
inline rocsparse_datatype rocsparse_compute_type<Kokkos::complex<double>>() {
  return rocsparse_datatype_f64_c;
}

template <typename T, typename E = void>
struct kokkos_to_rocsparse_type;

// for floats, rocsparse uses c++ builtin types
template <typename T>
struct kokkos_to_rocsparse_type<T, std::enable_if_t<std::is_floating_point_v<T>>> {
  using type = T;
};

// translate complex float
template <>
struct kokkos_to_rocsparse_type<Kokkos::complex<float>> {
  using type = rocsparse_float_complex;
};

// translate complex double
template <>
struct kokkos_to_rocsparse_type<Kokkos::complex<double>> {
  using type = rocsparse_double_complex;
};

// e.g. 5.4 -> 50400
#define KOKKOSSPARSE_IMPL_ROCM_VERSION ROCM_VERSION_MAJOR * 10000 + ROCM_VERSION_MINOR * 100 + ROCM_VERSION_PATCH

// Set the stream on the given rocSPARSE handle when this object
// is constructed, and reset to the default stream when this object is
// destructed.
struct TemporarySetRocsparseStream {
  TemporarySetRocsparseStream(rocsparse_handle handle_, const Kokkos::HIP& exec_) : handle(handle_) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_set_stream(handle, exec_.hip_stream()));
  }

  ~TemporarySetRocsparseStream() { KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_set_stream(handle, NULL)); }

  rocsparse_handle handle;
};

}  // namespace Impl

}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#endif  // _KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP
