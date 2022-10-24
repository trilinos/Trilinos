#ifndef __KOKKOSBATCHED_GEMV_TEAM_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMV_TEAM_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBlas1_set_impl.hpp"
#include "KokkosBlas1_team_scal_impl.hpp"
#include "KokkosBatched_InnerMultipleDotProduct_Serial_Impl.hpp"

namespace KokkosBatched {

///
/// Team Internal Impl
/// ====================
template <typename ArgAlgo>
struct TeamGemvInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const int m, const int n,
      const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
      const int as1, const ValueType *KOKKOS_RESTRICT x, const int xs0,
      const ScalarType beta,
      /**/ ValueType *KOKKOS_RESTRICT y, const int ys0);

  template <typename MemberType, typename ScalarType, typename layout,
            typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const int N, const int m, const int n,
      const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
      const int as1, const int as2, const ValueType *KOKKOS_RESTRICT x,
      const int xs0, const int xs1, const ScalarType beta,
      /**/ ValueType *KOKKOS_RESTRICT y, const int ys0, const int ys1);
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamGemvInternal<Algo::Gemv::Unblocked>::invoke(
    const MemberType &member, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    const ValueType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT y, const int ys0) {
  const ScalarType one(1.0), zero(0.0);

  // y = beta y + alpha A x
  // y (m), A(m x n), B(n)

  if (beta == zero)
    KokkosBlas::Impl::TeamSetInternal::invoke(member, m, zero, y, ys0);
  else if (beta != one)
    KokkosBlas::Impl::TeamScaleInternal::invoke(member, m, beta, y, ys0);

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    if (beta != one) member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, m),
                         [&](const int &i) {
                           ValueType t(0);
                           const ValueType *KOKKOS_RESTRICT tA = (A + i * as0);
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
                           for (int j = 0; j < n; ++j)
                             t += tA[j * as1] * x[j * xs0];
                           y[i * ys0] += alpha * t;
                         });
  }
  return 0;
}

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamGemvInternal<Algo::Gemv::Blocked>::invoke(
    const MemberType &member, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    const ValueType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT y, const int ys0) {
  const ScalarType one(1.0), zero(0.0);

  // y = beta y + alpha A x
  // y (m), A(m x n), B(n)

  constexpr int mbAlgo = Algo::Gemv::Blocked::mb();

  if (beta == zero)
    KokkosBlas::Impl::TeamSetInternal::invoke(member, m, zero, y, ys0);
  else if (beta != one)
    KokkosBlas::Impl::TeamScaleInternal::invoke(member, m, beta, y, ys0);

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    if (beta != one) member.team_barrier();

    InnerMultipleDotProduct<mbAlgo> inner(as0, as1, xs0, ys0);
    const int tsize = member.team_size();
    const int mb_a = m / tsize + (m % tsize > 0), mb_b = mbAlgo;
    // Made this non-const in order to WORKAROUND issue #349
    int mb = mb_a < mb_b ? mb_a : mb_b, mp = m % mb;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, (m / mb) + (mp > 0)),
                         [&](const int &ii) {
                           const int i = ii * mb;
                           inner.serial_invoke(alpha, A + i * as0, x,
                                               (i + mb) > m ? (m - i) : mb, n,
                                               y + i * ys0);
                         });
    member.team_barrier();
  }

  return 0;
}

template <>
template <typename MemberType, typename ScalarType, typename layout,
          typename ValueType>
KOKKOS_INLINE_FUNCTION int TeamGemvInternal<Algo::Gemv::Unblocked>::invoke(
    const MemberType &member, const int N, const int m, const int n,
    const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
    const int as1, const int as2, const ValueType *KOKKOS_RESTRICT X,
    const int xs0, const int xs1, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT Y, const int ys0, const int ys1) {
  const ScalarType one(1.0), zero(0.0);

  // y_l = beta y_l + alpha A_l x_l for l in range(0, N)
  // y_l (m), A_l(m x n), B_l(n)

  if (beta == zero)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, N * m),
                         [&](const int &iTemp) {
                           int iRow, iMatrix;
                           getIndices<int, layout>(iTemp, m, N, iRow, iMatrix);
                           Y[ys0 * iMatrix + ys1 * iRow] = zero;
                         });
  else if (beta != one)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, N * m),
                         [&](const int &iTemp) {
                           int iRow, iMatrix;
                           getIndices<int, layout>(iTemp, m, N, iRow, iMatrix);
                           Y[ys0 * iMatrix + ys1 * iRow] *= beta;
                         });

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    if (beta != one) member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, N * m),
                         [&](const int &iTemp) {
                           int iRow, iMatrix;
                           ValueType t(0);
                           getIndices<int, layout>(iTemp, m, N, iRow, iMatrix);
                           for (int i = 0; i < n; ++i)
                             t += A[as0 * iMatrix + as1 * iRow + as2 * i] *
                                  X[xs0 * iMatrix + xs1 * i];
                           Y[ys0 * iMatrix + ys1 * iRow] += alpha * t;
                         });
  }
  return 0;
}
}  // namespace KokkosBatched

#endif
