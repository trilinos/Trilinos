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

#ifndef _KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP
#define _KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"

namespace KokkosSparse {
namespace Impl {

inline void cusparse_internal_error_throw(cusparseStatus_t cusparseStatus,
                                          const char* name, const char* file,
                                          const int line) {
  std::ostringstream out;
#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)
  out << name << " error( " << cusparseGetErrorName(cusparseStatus)
      << "): " << cusparseGetErrorString(cusparseStatus);
#else
  out << name << " error( ";
  switch (cusparseStatus) {
    case CUSPARSE_STATUS_NOT_INITIALIZED:
      out << "CUSPARSE_STATUS_NOT_INITIALIZED): cusparse handle was not "
             "created correctly.";
      break;
    case CUSPARSE_STATUS_ALLOC_FAILED:
      out << "CUSPARSE_STATUS_ALLOC_FAILED): you might tried to allocate too "
             "much memory";
      break;
    case CUSPARSE_STATUS_INVALID_VALUE:
      out << "CUSPARSE_STATUS_INVALID_VALUE)";
      break;
    case CUSPARSE_STATUS_ARCH_MISMATCH:
      out << "CUSPARSE_STATUS_ARCH_MISMATCH)";
      break;
    case CUSPARSE_STATUS_MAPPING_ERROR:
      out << "CUSPARSE_STATUS_MAPPING_ERROR)";
      break;
    case CUSPARSE_STATUS_EXECUTION_FAILED:
      out << "CUSPARSE_STATUS_EXECUTION_FAILED)";
      break;
    case CUSPARSE_STATUS_INTERNAL_ERROR:
      out << "CUSPARSE_STATUS_INTERNAL_ERROR)";
      break;
    case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      out << "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED)";
      break;
    case CUSPARSE_STATUS_ZERO_PIVOT:
      out << "CUSPARSE_STATUS_ZERO_PIVOT)";
      break;
    default: out << "unrecognized error code): this is bad!"; break;
  }
#endif  // CUSPARSE_VERSION
  if (file) {
    out << " " << file << ":" << line;
  }
  throw std::runtime_error(out.str());
}

inline void cusparse_internal_safe_call(cusparseStatus_t cusparseStatus,
                                        const char* name,
                                        const char* file = nullptr,
                                        const int line   = 0) {
  if (CUSPARSE_STATUS_SUCCESS != cusparseStatus) {
    cusparse_internal_error_throw(cusparseStatus, name, file, line);
  }
}

// The macro below defines is the public interface for the safe cusparse calls.
// The functions themselves are protected by impl namespace.
#define KOKKOS_CUSPARSE_SAFE_CALL(call)                                  \
  KokkosSparse::Impl::cusparse_internal_safe_call(call, #call, __FILE__, \
                                                  __LINE__)

template <typename T>
cudaDataType cuda_data_type_from() {
  // Note:  compile-time failure is disabled to allow for packages such as
  // Ifpack2 to more easily support scalar types that cuSPARSE may not.

  // compile-time failure with a nice message if called on an unsupported type
  // static_assert(!std::is_same<T, T>::value,
  //               "cuSparse TPL does not support scalar type");
  // static_assert(false, ...) is allowed to error even if the code is not
  // instantiated. obfuscate the predicate Despite this function being
  // uncompilable, the compiler may decide that a return statement is missing,
  // so throw to silence that
  throw std::logic_error("unreachable throw after static_assert");
}

/* If half_t is not float, need to define a conversion for both
   otherwise, conversion for half_t IS conversion for float
*/
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
template <>
inline cudaDataType cuda_data_type_from<Kokkos::Experimental::half_t>() {
  return CUDA_R_16F;  // Kokkos half_t is a half
}
#endif
// half_t is defined to be float, so this works for both half_t and float when
// half_t is float
template <>
inline cudaDataType cuda_data_type_from<float>() {
  return CUDA_R_32F;  // Kokkos half_t is a float
}
template <>
inline cudaDataType cuda_data_type_from<double>() {
  return CUDA_R_64F;
}
template <>
inline cudaDataType cuda_data_type_from<Kokkos::complex<float>>() {
  return CUDA_C_32F;
}
template <>
inline cudaDataType cuda_data_type_from<Kokkos::complex<double>>() {
  return CUDA_C_64F;
}

#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)

template <typename T>
cusparseIndexType_t cusparse_index_type_t_from() {
#define AS_STR_LITERAL_IMPL_(x) #x
#define AS_STR_LITERAL(x) AS_STR_LITERAL_IMPL_(x)
  static_assert(!std::is_same<T, T>::value,
                "cuSparse " AS_STR_LITERAL(
                    CUSPARSE_VERSION) " TPL does not support index type");
  // static_assert(false, ...) is allowed to error even if the code is not
  // instantiated. obfuscate the predicate Despite this function being
  // uncompilable, the compiler may decide that a return statement is missing,
  // so throw to silence that
  throw std::logic_error("unreachable throw after static_assert");
#undef AS_STR_LITERAL_IMPL_
#undef AS_STR_LITERAL
}

template <>
inline cusparseIndexType_t cusparse_index_type_t_from<int>() {
  return CUSPARSE_INDEX_32I;
}
template <>
inline cusparseIndexType_t cusparse_index_type_t_from<int64_t>() {
  return CUSPARSE_INDEX_64I;
}
// Currently no CUSPARSE_INDEX_64U but this will work most of the time
template <>
inline cusparseIndexType_t cusparse_index_type_t_from<size_t>() {
  return CUSPARSE_INDEX_64I;
}
template <>
inline cusparseIndexType_t cusparse_index_type_t_from<unsigned short>() {
  return CUSPARSE_INDEX_16U;
}
#endif

}  // namespace Impl

}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#endif  // _KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP
