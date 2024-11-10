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
#ifndef __KOKKOSBATCHED_INVERSELU_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_INVERSELU_SERIAL_IMPL_HPP__

/// \author Vinh Dang (vqdang@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_SolveLU_Decl.hpp"
#include "KokkosBatched_SolveLU_Serial_Impl.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// =========

///
/// InverseLU no piv
///

#if defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL__) && defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_BATCHED__) && \
    defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_COMPACT_BATCHED__)
template <>
template <typename AViewType, typename WViewType>
KOKKOS_INLINE_FUNCTION int SerialInverseLU<Algo::InverseLU::CompactMKL>::invoke(const AViewType &A,
                                                                                const WViewType &W) {
  typedef typename AViewType::value_type vector_type;
  // typedef typename vector_type::value_type value_type;

  const int m = A.extent(0), n = A.extent(1);

  static_assert(is_vector<vector_type>::value, "value type is not vector type");
  static_assert(vector_type::vector_length == 4 || vector_type::vector_length == 8,
                "AVX, AVX2 and AVX512 is supported");
  static_assert(AViewType::rank == 2, "A should have two dimensions");
  static_assert(WViewType::rank == 1, "W should have one dimension");
  static_assert(std::is_same<typename AViewType::memory_space, typename WViewType::memory_space>::value,
                "A and W should be on the same memory space");
  static_assert(!std::is_same<typename WViewType::array_layout, Kokkos::LayoutStride>::value,
                "W should be an contiguous 1D array");
  assert(A.extent(0) * A.extent(1) * sizeof(typename AViewType::value_type) <=
         W.span() * sizeof(typename WViewType::value_type));
  assert(m == n);

  const MKL_COMPACT_PACK format = vector_type::vector_length == 8 ? MKL_COMPACT_AVX512 : MKL_COMPACT_AVX;

  int r_val = 0;
  if (A.stride(0) == 1) {
    mkl_dgetrinp_compact(MKL_COL_MAJOR, n, (double *)A.data(), A.stride(1), (double *)W.data(),
                         (MKL_INT)(n * n * vector_type::vector_length), (MKL_INT *)&r_val, format,
                         (MKL_INT)vector_type::vector_length);

  } else if (A.stride(1) == 1) {
    mkl_dgetrinp_compact(MKL_ROW_MAJOR, n, (double *)A.data(), A.stride(0), (double *)W.data(),
                         (MKL_INT)(n * n * vector_type::vector_length), (MKL_INT *)&r_val, format,
                         (MKL_INT)vector_type::vector_length);
  } else {
    r_val = -1;
  }
  return r_val;
}
#endif

}  // namespace KokkosBatched

#endif
