#ifndef __KOKKOSBATCHED_GEMV_TEAMVECTOR_INTERNAL_HPP__
#define __KOKKOSBATCHED_GEMV_TEAMVECTOR_INTERNAL_HPP__

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
struct TeamVectorGemvInternal {
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType & /*member*/, const int /*m*/, const int /*n*/,
      const ScalarType /*alpha*/, const ValueType *KOKKOS_RESTRICT /*A*/,
      const int /*as0*/, const int /*as1*/,
      const ValueType *KOKKOS_RESTRICT /*x*/, const int /*xs0*/,
      const ScalarType /*beta*/,
      /**/ ValueType *KOKKOS_RESTRICT /*y*/, const int /*ys0*/) {
    assert(false && "Error: encounter dummy impl");
    return 0;
  }
  template <typename MemberType, typename ScalarType, typename layout,
            typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType & /*member*/, const int /*N*/, const int /*m*/,
      const int /*n*/, const ScalarType /*alpha*/,
      const ValueType *KOKKOS_RESTRICT /*A*/, const int /*as0*/,
      const int /*as1*/, const int /*as2*/,
      const ValueType *KOKKOS_RESTRICT /*x*/, const int /*xs0*/,
      const int /*xs1*/, const ScalarType /*beta*/,
      /**/ ValueType *KOKKOS_RESTRICT /*y*/, const int /*ys0*/,
      const int /*ys1*/) {
    assert(false && "Error: encounter dummy impl");
    return 0;
  }
};

template <>
template <typename MemberType, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int
TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(
    const MemberType &member, const int m, const int n, const ScalarType alpha,
    const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    const ValueType *KOKKOS_RESTRICT x, const int xs0, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT y, const int ys0) {
  const ScalarType one(1.0), zero(0.0);

  // y = beta y + alpha A x
  // y (m), A(m x n), B(n)

  if (beta == zero)
    KokkosBlas::Impl::TeamVectorSetInternal::invoke(member, m, zero, y, ys0);
  else if (beta != one)
    KokkosBlas::Impl::TeamVectorScaleInternal::invoke(member, m, beta, y, ys0);

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    if (beta != one) member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, m), [&](const int &i) {
      ValueType t(0);
      const ValueType *KOKKOS_RESTRICT tA = (A + i * as0);
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(member, n),
          [&](const int &j, ValueType &update) {
            update += tA[j * as1] * x[j * xs0];
          },
          t);
      Kokkos::single(Kokkos::PerThread(member),
                     [&]() { y[i * ys0] += alpha * t; });
    });
  }
  return 0;
}

template <>
template <typename MemberType, typename ScalarType, typename layout,
          typename ValueType>
KOKKOS_INLINE_FUNCTION int
TeamVectorGemvInternal<Algo::Gemv::Unblocked>::invoke(
    const MemberType &member, const int N, const int m, const int n,
    const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A, const int as0,
    const int as1, const int as2, const ValueType *KOKKOS_RESTRICT X,
    const int xs0, const int xs1, const ScalarType beta,
    /**/ ValueType *KOKKOS_RESTRICT Y, const int ys0, const int ys1) {
  const ScalarType one(1.0), zero(0.0);

  // y_l = beta y_l + alpha A_l x_l for l in range(0, N)
  // y_l (m), A_l(m x n), B_l(n)

  if (beta == zero)
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, N * m),
                         [&](const int &iTemp) {
                           int iRow, iMatrix;
                           getIndices<int, layout>(iTemp, m, N, iRow, iMatrix);
                           Y[ys0 * iMatrix + ys1 * iRow] = zero;
                         });
  else if (beta != one)
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, N * m),
                         [&](const int &iTemp) {
                           int iRow, iMatrix;
                           getIndices<int, layout>(iTemp, m, N, iRow, iMatrix);
                           Y[ys0 * iMatrix + ys1 * iRow] *= beta;
                         });

  if (alpha != zero) {
    if (m <= 0 || n <= 0) return 0;

    if (beta != one) member.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, N * m),
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
