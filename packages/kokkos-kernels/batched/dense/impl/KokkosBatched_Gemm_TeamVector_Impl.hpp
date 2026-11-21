// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GEMM_TEAMVECTOR_IMPL_HPP
#define KOKKOSBATCHED_GEMM_TEAMVECTOR_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBlas_util.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemm_Common_Impl.hpp"
#include "KokkosBatched_Gemm_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============

template <typename MemberType, typename ArgTransA, typename ArgTransB, typename ArgAlgo>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorGemm<MemberType, ArgTransA, ArgTransB, ArgAlgo>::invoke(
    const MemberType &member, const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta,
    const CViewType &C) {
  const int k = std::same_as<ArgTransA, Trans::NoTranspose> ? Impl::get_extent_int(A, 1) : Impl::get_extent_int(A, 0);
  const int m = Impl::get_extent_int(C, 0), n = Impl::get_extent_int(C, 1);

  // Quick return if possible
  if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkGemmInput<ArgTransA, ArgTransB>(A, B, C);
  if (info) return info;

  using opA = std::conditional_t<std::same_as<ArgTransA, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                 KokkosBlas::Impl::OpID>;
  using opB = std::conditional_t<std::same_as<ArgTransB, Trans::ConjTranspose>, KokkosBlas::Impl::OpConj,
                                 KokkosBlas::Impl::OpID>;

  const std::size_t A_stride_0 =
      std::same_as<ArgTransA, Trans::NoTranspose> ? Impl::get_stride(A, 0) : Impl::get_stride(A, 1);
  const std::size_t A_stride_1 =
      std::same_as<ArgTransA, Trans::NoTranspose> ? Impl::get_stride(A, 1) : Impl::get_stride(A, 0);
  const std::size_t B_stride_0 =
      std::same_as<ArgTransB, Trans::NoTranspose> ? Impl::get_stride(B, 0) : Impl::get_stride(B, 1);
  const std::size_t B_stride_1 =
      std::same_as<ArgTransB, Trans::NoTranspose> ? Impl::get_stride(B, 1) : Impl::get_stride(B, 0);
  const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

  // C = beta C + alpha op(A) op(B)
  return Impl::TeamVectorGemmInternal<ArgAlgo>::invoke(member, opA(), opB(), m, n, k, alpha, A.data(), A_stride_0,
                                                       A_stride_1, B.data(), B_stride_0, B_stride_1, beta, C.data(),
                                                       C_stride_0, C_stride_1);
}

}  // namespace KokkosBatched

#endif
