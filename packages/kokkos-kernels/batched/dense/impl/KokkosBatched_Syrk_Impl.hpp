// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SYRK_IMPL_HPP_
#define KOKKOSBATCHED_SYRK_IMPL_HPP_

#include <KokkosBlas_util.hpp>
#include <KokkosBatched_Util.hpp>
#include "KokkosBatched_Syrk_Internal.hpp"

namespace KokkosBatched {
namespace Impl {
template <bool is_trans, typename AViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION static int checkSyrkInput([[maybe_unused]] const AViewType &A,
                                                 [[maybe_unused]] const CViewType &C) {
  static_assert(Kokkos::is_view_v<AViewType>, "KokkosBatched::syrk: AViewType is not a Kokkos::View.");
  static_assert(Kokkos::is_view_v<CViewType>, "KokkosBatched::syrk: CViewType is not a Kokkos::View.");
  static_assert(AViewType::rank == 2, "KokkosBatched::syrk: AViewType must have rank 2.");
  static_assert(CViewType::rank == 2, "KokkosBatched::syrk: CViewType must have rank 2.");
#ifndef NDEBUG
  const int ldc = C.extent_int(0), n = C.extent_int(1);
  const int lda   = A.extent_int(0);
  const int k     = is_trans ? A.extent_int(0) : A.extent_int(1);
  const int nrowa = is_trans ? k : n;

  if (lda < Kokkos::max(1, nrowa)) {
    Kokkos::printf(
        "KokkosBatched::syrk: leading dimension of A must not be smaller than "
        "max(1, nrowa): "
        "lda = %d, nrowa = %d\n",
        lda, nrowa);
    return 1;
  }

  if (ldc < Kokkos::max(1, n)) {
    Kokkos::printf(
        "KokkosBatched::syrk: leading dimension of C must not be smaller than "
        "max(1, n): "
        "ldc = %d, n = %d\n",
        ldc, n);
    return 1;
  }
#endif
  return 0;
}
}  // namespace Impl

/// Serial Batched Syrk:
/// {s,d,c,z}syrk and {c,z}herk interface
/// Performs one of the symmetric rank k operation
///   C := alpha*A*A**T + beta*C
///   C := alpha*A*A**H + beta*C
///   C := alpha*A**T*A + beta*C
///   C := alpha*A**H*A + beta*C
template <typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename AViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int SerialSyrk<ArgUplo, ArgTrans>::invoke(const ScalarType alpha, const AViewType &A,
                                                                 const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  constexpr bool is_trans = std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>;
  const int n             = C.extent_int(1);
  const int k             = is_trans ? A.extent_int(0) : A.extent_int(1);
  if (n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSyrkInput<is_trans>(A, C);
  if (info) return info;

  using op = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                 std::same_as<ArgTrans, Trans::ConjTranspose>),
                                KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;

  using sym_op = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                     std::same_as<ArgTrans, Trans::ConjTranspose>),
                                    KokkosBlas::Impl::OpReal, KokkosBlas::Impl::OpID>;

  using value_type = typename CViewType::non_const_value_type;
  using mag_type   = typename KokkosKernels::ArithTraits<value_type>::mag_type;
  using ReduceType = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                         std::same_as<ArgTrans, Trans::ConjTranspose>),
                                        mag_type, value_type>;

  if constexpr (std::same_as<ArgUplo, Uplo::Lower>) {
    return Impl::SerialSyrkInternalLower::invoke<is_trans, ReduceType>(
        op(), sym_op(), n, k, alpha, A.data(), A.stride(0), A.stride(1), beta, C.data(), C.stride(0), C.stride(1));
  } else {
    return Impl::SerialSyrkInternalUpper::invoke<is_trans, ReduceType>(
        op(), sym_op(), n, k, alpha, A.data(), A.stride(0), A.stride(1), beta, C.data(), C.stride(0), C.stride(1));
  }
}

/// Team Batched Syrk:
/// {s,d,c,z}syrk and {c,z}herk interface
/// Performs one of the symmetric rank k operation
///   C := alpha*A*A**T + beta*C
///   C := alpha*A*A**H + beta*C
///   C := alpha*A**T*A + beta*C
///   C := alpha*A**H*A + beta*C
template <typename MemberType, typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename AViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int TeamSyrk<MemberType, ArgUplo, ArgTrans>::invoke(const MemberType &member,
                                                                           const ScalarType alpha, const AViewType &A,
                                                                           const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  constexpr bool is_trans = std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>;
  const int n             = C.extent_int(1);
  const int k             = is_trans ? A.extent_int(0) : A.extent_int(1);
  if (n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSyrkInput<is_trans>(A, C);
  if (info) return info;

  using op = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                 std::same_as<ArgTrans, Trans::ConjTranspose>),
                                KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;

  using sym_op = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                     std::same_as<ArgTrans, Trans::ConjTranspose>),
                                    KokkosBlas::Impl::OpReal, KokkosBlas::Impl::OpID>;

  using value_type = typename CViewType::non_const_value_type;
  using mag_type   = typename KokkosKernels::ArithTraits<value_type>::mag_type;
  using ReduceType = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                         std::same_as<ArgTrans, Trans::ConjTranspose>),
                                        mag_type, value_type>;

  if constexpr (std::same_as<ArgUplo, Uplo::Lower>) {
    return Impl::TeamSyrkInternalLower::invoke<is_trans, ReduceType>(member, op(), sym_op(), n, k, alpha, A.data(),
                                                                     A.stride(0), A.stride(1), beta, C.data(),
                                                                     C.stride(0), C.stride(1));
  } else {
    return Impl::TeamSyrkInternalUpper::invoke<is_trans, ReduceType>(member, op(), sym_op(), n, k, alpha, A.data(),
                                                                     A.stride(0), A.stride(1), beta, C.data(),
                                                                     C.stride(0), C.stride(1));
  }
}

/// TeamVector Batched Syrk:
/// {s,d,c,z}syrk and {c,z}herk interface
/// Performs one of the symmetric rank k operation
///   C := alpha*A*A**T + beta*C
///   C := alpha*A*A**H + beta*C
///   C := alpha*A**T*A + beta*C
///   C := alpha*A**H*A + beta*C
template <typename MemberType, typename ArgUplo, typename ArgTrans>
template <typename ScalarType, typename AViewType, typename CViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorSyrk<MemberType, ArgUplo, ArgTrans>::invoke(
    const MemberType &member, const ScalarType alpha, const AViewType &A, const ScalarType beta, const CViewType &C) {
  // Quick return if possible
  constexpr bool is_trans = std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>;
  const int n             = C.extent_int(1);
  const int k             = is_trans ? A.extent_int(0) : A.extent_int(1);
  if (n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

  auto info = Impl::checkSyrkInput<is_trans>(A, C);
  if (info) return info;

  using op = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                 std::same_as<ArgTrans, Trans::ConjTranspose>),
                                KokkosBlas::Impl::OpConj, KokkosBlas::Impl::OpID>;

  using sym_op = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                     std::same_as<ArgTrans, Trans::ConjTranspose>),
                                    KokkosBlas::Impl::OpReal, KokkosBlas::Impl::OpID>;

  using value_type = typename CViewType::non_const_value_type;
  using mag_type   = typename KokkosKernels::ArithTraits<value_type>::mag_type;
  using ReduceType = std::conditional_t<(std::same_as<ArgTrans, Trans::ConjNoTranspose> ||
                                         std::same_as<ArgTrans, Trans::ConjTranspose>),
                                        mag_type, value_type>;

  if constexpr (std::same_as<ArgUplo, Uplo::Lower>) {
    return Impl::TeamVectorSyrkInternalLower::invoke<is_trans, ReduceType>(member, op(), sym_op(), n, k, alpha,
                                                                           A.data(), A.stride(0), A.stride(1), beta,
                                                                           C.data(), C.stride(0), C.stride(1));
  } else {
    return Impl::TeamVectorSyrkInternalUpper::invoke<is_trans, ReduceType>(member, op(), sym_op(), n, k, alpha,
                                                                           A.data(), A.stride(0), A.stride(1), beta,
                                                                           C.data(), C.stride(0), C.stride(1));
  }
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_SYRK_IMPL_HPP_
