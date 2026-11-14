// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SYMM_IMPL_HPP
#define KOKKOSBATCHED_SYMM_IMPL_HPP

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include <concepts>
#include <type_traits>
#include "KokkosBatched_Symm_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <typename ArgSide, typename ArgUplo, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION static int checkSymmInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const BViewType &B,
                                                 [[maybe_unused]] const CViewType &C) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::symm: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<BViewType>, "KokkosBatched::symm: BViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<CViewType>, "KokkosBatched::symm: CViewType is not a Kokkos::View.");
  static_assert(!std::is_const_v<typename CViewType::value_type>,
                "KokkosBatched::symm: CViewType must have non-const value type.");

  static_assert(AViewType::rank == 2, "KokkosBatched::symm: AViewType must have rank 2.");
  static_assert(BViewType::rank == 2, "KokkosBatched::symm: BViewType must have rank 2.");
  static_assert(CViewType::rank == 2, "KokkosBatched::symm: CViewType must have rank 2.");

#ifndef NDEBUG
  const int m = C.extent_int(0), n = C.extent_int(1);
  const int nrowa = std::same_as<ArgSide, Side::Left> ? m : n;
  const int lda = A.extent_int(0), ldb = B.extent_int(0);

  if (lda < Kokkos::max(1, nrowa)) {
    Kokkos::printf(
        "KokkosBatched::symm: leading dimension of A must not be smaller than "
        "max(1, nrowa): "
        "lda = %d, nrowa = %d\n",
        lda, nrowa);
    return 1;
  }
  if (ldb < Kokkos::max(1, m)) {
    Kokkos::printf(
        "KokkosBatched::symm: leading dimension of B must not be smaller than "
        "max(1, m): "
        "ldb = %d, m = %d\n",
        ldb, m);
    return 1;
  }

#endif

  return 0;
}
}  // namespace Impl

///
/// Serial Impl
/// ===========

template <typename ArgSide, typename ArgUplo, bool Conj>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialSymm<ArgSide, ArgUplo, Conj>::invoke(const ScalarType alpha, const AViewType &A,
                                                                      const BViewType &B, const ScalarType beta,
                                                                      const CViewType &C) {
  const int m = C.extent_int(0), n = C.extent_int(1);

  // Quick return if possible
  if (m == 0 || n == 0 || (alpha == ScalarType(0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSymmInput<ArgSide, ArgUplo>(A, B, C);
  if (info) return info;

  using op     = std::conditional_t<Conj, KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<Conj, KokkosBlas::Impl::OpReal, KokkosBlas::Impl::OpID>;

  return Impl::SerialSymmInternal::invoke<ArgSide, ArgUplo>(op(), sym_op(), m, n, alpha, A.data(), A.stride(0),
                                                            A.stride(1), B.data(), B.stride(0), B.stride(1), beta,
                                                            C.data(), C.stride(0), C.stride(1));
}

///
/// Team Impl
/// =========

template <typename MemberType, typename ArgSide, typename ArgUplo, bool Conj>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int TeamSymm<MemberType, ArgSide, ArgUplo, Conj>::invoke(const MemberType &member,
                                                                                const ScalarType alpha,
                                                                                const AViewType &A, const BViewType &B,
                                                                                const ScalarType beta,
                                                                                const CViewType &C) {
  const int m = C.extent_int(0), n = C.extent_int(1);

  // Quick return if possible
  if (m == 0 || n == 0 || (alpha == ScalarType(0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSymmInput<ArgSide, ArgUplo>(A, B, C);
  if (info) return info;

  using op     = std::conditional_t<Conj, KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<Conj, KokkosBlas::Impl::OpReal, KokkosBlas::Impl::OpID>;

  return Impl::TeamSymmInternal::invoke<ArgSide, ArgUplo>(member, op(), sym_op(), m, n, alpha, A.data(), A.stride(0),
                                                          A.stride(1), B.data(), B.stride(0), B.stride(1), beta,
                                                          C.data(), C.stride(0), C.stride(1));
}

///
/// TeamVector Impl
/// ===============

template <typename MemberType, typename ArgSide, typename ArgUplo, bool Conj>
template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorSymm<MemberType, ArgSide, ArgUplo, Conj>::invoke(
    const MemberType &member, const ScalarType alpha, const AViewType &A, const BViewType &B, const ScalarType beta,
    const CViewType &C) {
  const int m = C.extent_int(0), n = C.extent_int(1);

  // Quick return if possible
  if (m == 0 || n == 0 || (alpha == ScalarType(0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSymmInput<ArgSide, ArgUplo>(A, B, C);
  if (info) return info;

  using op     = std::conditional_t<Conj, KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;
  using sym_op = std::conditional_t<Conj, KokkosBlas::Impl::OpReal, KokkosBlas::Impl::OpID>;

  return Impl::TeamVectorSymmInternal::invoke<ArgSide, ArgUplo>(member, op(), sym_op(), m, n, alpha, A.data(),
                                                                A.stride(0), A.stride(1), B.data(), B.stride(0),
                                                                B.stride(1), beta, C.data(), C.stride(0), C.stride(1));
}

}  // namespace KokkosBatched

#endif
