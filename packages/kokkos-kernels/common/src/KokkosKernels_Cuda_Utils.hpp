// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_CUDA_UTILS_HPP
#define KOKKOSKERNELS_CUDA_UTILS_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Half.hpp>
#include <stdexcept>

#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda.h>

namespace KokkosKernels::Impl {

template <typename T>
cudaDataType cuda_data_type_from() {
  // Note:  compile-time failure is disabled to allow for packages such as
  // Ifpack2 to more easily support scalar types that cuBLAS/cuSPARSE/cuSPARSE may not.

  // compile-time failure with a nice message if called on an unsupported type
  // static_assert(!std::is_same<T, T>::value,
  //               "CUDA TPLs do not support scalar type");
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

}  // namespace KokkosKernels::Impl

#endif  // KOKKOS_ENABLE_CUDA
#endif  // KOKKOSKERNELS_CUDA_UTILS_HPP
