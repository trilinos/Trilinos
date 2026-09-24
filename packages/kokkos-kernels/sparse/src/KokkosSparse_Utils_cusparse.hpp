// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP
#define KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cuda.h"
#include "cuda_runtime.h"
#include "cusparse.h"

#include <KokkosKernels_Cuda_Utils.hpp>

namespace KokkosSparse {
namespace Impl {

inline void cusparse_internal_error_throw(cusparseStatus_t cusparseStatus, const char* name, const char* file,
                                          const int line) {
  std::ostringstream out;
  out << name << " error( " << cusparseGetErrorName(cusparseStatus) << "): " << cusparseGetErrorString(cusparseStatus);
  if (file) {
    out << " " << file << ":" << line;
  }
  throw std::runtime_error(out.str());
}

inline void cusparse_internal_safe_call(cusparseStatus_t cusparseStatus, const char* name, const char* file = nullptr,
                                        const int line = 0) {
  if (CUSPARSE_STATUS_SUCCESS != cusparseStatus) {
    cusparse_internal_error_throw(cusparseStatus, name, file, line);
  }
}

#define KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(call) \
  KokkosSparse::Impl::cusparse_internal_safe_call(call, #call, __FILE__, __LINE__)

//  [[deprecated("Please use KokkosKernels::Impl::cuda_data_type_from() in KokkosKernels_Cuda_Utils.hpp")]]
template <typename T>
cudaDataType cuda_data_type_from() {
  return KokkosKernels::Impl::cuda_data_type_from<T>();
}

template <typename T>
cusparseIndexType_t cusparse_index_type_t_from() {
#define AS_STR_LITERAL_IMPL_(x) #x
#define AS_STR_LITERAL(x) AS_STR_LITERAL_IMPL_(x)
  static_assert(!std::is_same<T, T>::value,
                "cuSparse " AS_STR_LITERAL(CUSPARSE_VERSION) " TPL does not support index type");
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

// Set the stream on the given cuSPARSE handle when this object
// is constructed, and reset to the default stream when this object is
// destructed.
struct TemporarySetCusparseStream {
  TemporarySetCusparseStream(cusparseHandle_t handle_, const Kokkos::Cuda& exec_) : handle(handle_) {
    KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseSetStream(handle, exec_.cuda_stream()));
  }

  ~TemporarySetCusparseStream() { KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseSetStream(handle, NULL)); }

  cusparseHandle_t handle;
};

}  // namespace Impl

}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#endif  // KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP
