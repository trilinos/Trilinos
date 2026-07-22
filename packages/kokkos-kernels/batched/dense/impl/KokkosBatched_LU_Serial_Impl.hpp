// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_LU_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_LU_SERIAL_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_LU_Serial_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// =========

///
/// SerialLU no piv
///

#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
template <>
template <typename AViewType>
KOKKOS_INLINE_FUNCTION int SerialLU<Algo::LU::CompactMKL>::invoke(
    const AViewType &A, const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny) {
  typedef typename AViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = A.extent(0), n = A.extent(1);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  int r_val = 0;
  if (A.stride(0) == 1) {
    mkl_dgetrfnp_compact(MKL_COL_MAJOR, m, n, (double *)A.data(), A.stride(1), (MKL_INT *)&r_val, format,
                         (MKL_INT)vector_type::vector_length);
  } else if (A.stride(1) == 1) {
    mkl_dgetrfnp_compact(MKL_ROW_MAJOR, m, n, (double *)A.data(), A.stride(0), (MKL_INT *)&r_val, format,
                         (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

template <>
template <typename AViewType>
KOKKOS_INLINE_FUNCTION int SerialLU<Algo::LU::Unblocked>::invoke(
    const AViewType &A, const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny) {
  return SerialLU_Internal<Algo::LU::Unblocked>::invoke(A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1),
                                                        tiny);
}

template <>
template <typename AViewType>
KOKKOS_INLINE_FUNCTION int SerialLU<Algo::LU::Blocked>::invoke(
    const AViewType &A, const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny) {
  return SerialLU_Internal<Algo::LU::Blocked>::invoke(A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1),
                                                      tiny);
}

}  // namespace KokkosBatched

#endif
