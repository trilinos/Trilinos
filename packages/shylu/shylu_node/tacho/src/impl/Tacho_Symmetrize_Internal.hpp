#ifndef __TACHO_SYMMETRIZE_INTERNAL_HPP__
#define __TACHO_SYMMETRIZE_INTERNAL_HPP__


/// \file  Tacho_Symmetrize_Internal.hpp
/// \brief Symmetrize a square block matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<>
  struct Symmetrize<Uplo::Upper,Algo::Internal> {
    template<typename MemberType,
             typename ViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(MemberType &member,
           const ViewTypeA &A) {
      typedef typename ViewTypeA::non_const_value_type value_type;        

      if (A.extent(0) == A.extent(1)) {
        if (A.span() > 0) {
#if defined(__CUDA_ARCH__)
          Kokkos::parallel_for
            (Kokkos::TeamThreadRange(member, A.extent(1)),
             [&](const ordinal_type &j) {
              Kokkos::parallel_for
                (Kokkos::ThreadVectorRange(member, A.extent(0)),
                 [&](const ordinal_type &i) {
                  A(i,j) = i > j ? A(j,i) : A(i,j);
                });
            });
#else
          for (ordinal_type j=0,jend=A.extent(1);j<jend;++j)
            for (ordinal_type i=0,iend=A.extent(0);i<iend;++i)
              A(i,j) = i > j ? A(j,i) : A(i,j);
#endif
        }
      } else {
        printf("Error: Symmetrize<Algo::Internal> A is not square\n");
      }
      return 0;
    }
  };

  template<>
  struct Symmetrize<Uplo::Lower,Algo::Internal> {
    template<typename MemberType,
             typename ViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(MemberType &member,
           const ViewTypeA &A) {
      typedef typename ViewTypeA::non_const_value_type value_type;        

      if (A.extent(0) == A.extent(1)) {
        if (A.span() > 0) {
#if defined(__CUDA_ARCH__)
          Kokkos::parallel_for
            (Kokkos::TeamThreadRange(member, A.extent(1)),
             [&](const ordinal_type &j) {
              Kokkos::parallel_for
                (Kokkos::ThreadVectorRange(member, A.extent(0)),
                 [&](const ordinal_type &i) {
                  A(i,j) = i < j ? A(j,i) : A(i,j);
                });
            });
#else
          for (ordinal_type j=0,jend=A.extent(1);j<jend;++j)
            for (ordinal_type i=0,iend=A.extent(0);i<iend;++i)
              A(i,j) = i < j ? A(j,i) : A(i,j);
#endif
        }
      } else {
        printf("Error: Symmetrize<Algo::Internal> A is not square\n");
      }
      return 0;
    }
  };

}
#endif
