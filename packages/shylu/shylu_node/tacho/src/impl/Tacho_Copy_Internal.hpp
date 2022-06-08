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
    // const ordinal_type sA = A.span(), sB = B.span();
    if (A.extent(0) == B.extent(0) && A.extent(0) == B.extent(0)) {
      if (A.span() > 0) {
#if defined(__CUDA_ARCH__)
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, A.extent(1)), [&](const ordinal_type &j) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, A.extent(0)),
                               [&](const ordinal_type &i) { A(i, j) = B(i, j); });
        });
#else
        if (A.span() == (A.extent(0) * A.extent(1)) && B.span() == (B.extent(0) * B.extent(1)))
          memcpy((void *)A.data(), (const void *)B.data(), A.span() * sizeof(value_type));
        else
          for (ordinal_type j = 0, jend = A.extent(1); j < jend; ++j)
            for (ordinal_type i = 0, iend = A.extent(0); i < iend; ++i)
              A(i, j) = B(i, j);
#endif
      }
    } else {
      printf("Error: Copy<Algo::Internal> A and B dimensions are not same\n");
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
#if defined(__CUDA_ARCH__)
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, A.extent(1)), [&](const ordinal_type &j) {
          const ordinal_type tmp = diagB.param == 'U' ? j : j + 1;
          const ordinal_type iend = tmp < A.extent(0) ? tmp : A.extent(0);
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, iend),
                               [&](const ordinal_type &i) { A(i, j) = B(i, j); });
        });
#else
        for (ordinal_type j = 0, jend = A.extent(1); j < jend; ++j) {
          const ordinal_type tmp = diagB.param == 'U' ? j : j + 1;
          const ordinal_type iend = tmp < A.extent(0) ? tmp : A.extent(0);
          for (ordinal_type i = 0; i < iend; ++i)
            A(i, j) = B(i, j);
        }
#endif
      } else if (uploB.param == 'L') {
#if defined(__CUDA_ARCH__)
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, A.extent(1)), [&](const ordinal_type &j) {
          const ordinal_type ibeg = diagB.param == 'U' ? j + 1 : j;
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, ibeg, A.extent(0)),
                               [&](const ordinal_type &i) { A(i, j) = B(i, j); });
        });
#else
        for (ordinal_type j = 0, jend = A.extent(1); j < jend; ++j) {
          const ordinal_type ibeg = diagB.param == 'U' ? j + 1 : j;
          for (ordinal_type i = ibeg, iend = A.extent(0); i < iend; ++i)
            A(i, j) = B(i, j);
        }
#endif
      }
    } else {
      printf("Error: Copy<Algo::Internal> A and B dimensions are not same\n");
    }
    return 0;
  }
};

} // namespace Tacho
#endif
