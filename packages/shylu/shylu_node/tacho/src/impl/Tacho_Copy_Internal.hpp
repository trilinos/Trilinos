// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_COPY_INTERNAL_HPP__
#define __TACHO_COPY_INTERNAL_HPP__

/// \file  Tacho_Copy_Internal.hpp
/// \brief Copy
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <> struct Copy<Algo::Internal> {
  template <typename MemberType, typename ViewTypeA, typename ViewTypeB>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeB &B) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeB::non_const_value_type value_type_b;
    static_assert(std::is_same<value_type, value_type_b>::value, "A and B does not have the value_type.");

    /// this should be for contiguous array
    const ordinal_type sA = A.span(), mA = A.extent(0), nA = A.extent(1), as0 = A.stride(0), as1 = A.stride(1);
    const ordinal_type /*sB = B.span(), */ mB = B.extent(0), nB = B.extent(1), bs0 = B.stride(0), bs1 = B.stride(1);
    if (mA == mB && nA == nB) {
      if (sA == (mA * nA) && as0 == bs0 && as1 == bs1) {
        /// contiguous array
        value_type *ptrA(A.data());
        const value_type *ptrB(B.data());
        KOKKOS_IF_ON_DEVICE((
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, sA),
                             [ptrA, ptrB](const ordinal_type &ij) { ptrA[ij] = ptrB[ij]; });
        ))
        KOKKOS_IF_ON_HOST((memcpy((void *)ptrA, (const void *)ptrB, sA * sizeof(value_type));))
      } else {
        KOKKOS_IF_ON_DEVICE((
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, mA * nA), [A, B, mA](const ordinal_type &ij) {
          const ordinal_type i = ij % mA, j = ij / mA;
          A(i, j) = B(i, j);
        });))
        KOKKOS_IF_ON_HOST((
        for (ordinal_type j = 0; j < nA; ++j)
          for (ordinal_type i = 0; i < mA; ++i)
            A(i, j) = B(i, j);))
      }
    } else {
      Kokkos::printf("Error: Copy<Algo::Internal> A and B dimensions are not same\n");
    }
    return 0;
  }

  template <typename MemberType, typename UploType, typename DiagType, typename ViewTypeA, typename ViewTypeB>
  KOKKOS_INLINE_FUNCTION static int invoke(MemberType &member, const ViewTypeA &A, const UploType uploB,
                                           const DiagType diagB, const ViewTypeB &B) {
    using value_type = typename ViewTypeA::non_const_value_type;
    using value_type_b = typename ViewTypeB::non_const_value_type;
    static_assert(std::is_same<value_type, value_type_b>::value, "A and B does not have the value_type.");

    /// this should be for contiguous array
    // const ordinal_type sA = A.span(), sB = B.span();
    if (A.extent(0) == B.extent(0) && A.extent(0) == B.extent(0) && A.span() > 0) {
      if (uploB.param == 'U') {
        KOKKOS_IF_ON_DEVICE((
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, A.extent(1)), [&](const ordinal_type &j) {
          const ordinal_type tmp = diagB.param == 'U' ? j : j + 1;
          const ordinal_type iend = tmp < A.extent(0) ? tmp : A.extent(0);
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, iend),
                               [&](const ordinal_type &i) { A(i, j) = B(i, j); });
        });))
        KOKKOS_IF_ON_HOST((
          for (ordinal_type j = 0, jend = A.extent(1); j < jend; ++j) {
          const ordinal_type tmp = diagB.param == 'U' ? j : j + 1;
          const ordinal_type iend = tmp < A.extent(0) ? tmp : A.extent(0);
          for (ordinal_type i = 0; i < iend; ++i)
            A(i, j) = B(i, j);
        }))
      } else if (uploB.param == 'L') {
        KOKKOS_IF_ON_DEVICE((
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, A.extent(1)), [&](const ordinal_type &j) {
          const ordinal_type ibeg = diagB.param == 'U' ? j + 1 : j;
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, ibeg, A.extent(0)),
                               [&](const ordinal_type &i) { A(i, j) = B(i, j); });
        });))
        KOKKOS_IF_ON_HOST((
        for (ordinal_type j = 0, jend = A.extent(1); j < jend; ++j) {
          const ordinal_type ibeg = diagB.param == 'U' ? j + 1 : j;
          for (ordinal_type i = ibeg, iend = A.extent(0); i < iend; ++i)
            A(i, j) = B(i, j);
        }))
      }
    } else {
      printf("Error: Copy<Algo::Internal> A and B dimensions are not same\n");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
