#ifndef __TACHO_SET_IDENTITY_INTERNAL_HPP__
#define __TACHO_SET_IDENTITY_INTERNAL_HPP__


/// \file  Tacho_SetIdentity_Internal.hpp
/// \brief Set an identity matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {
  
  template<>
  struct SetIdentity<Algo::Internal> {
    template<typename MemberType,
             typename ViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(MemberType &member,
           const ViewTypeA &A) {
      typedef typename ViewTypeA::non_const_value_type value_type;        
      static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");

      const ordinal_type 
        m = A.extent(0),
        n = A.extent(1);

      if (m > 0 && n > 0) {
        const value_type one(1), zero(0);
        Kokkos::parallel_for
          (Kokkos::TeamThreadRange(member, n),
           [&](const ordinal_type &j) {
            Kokkos::parallel_for
              (Kokkos::ThreadVectorRange(member, m),
               [&](const ordinal_type &i) {
                A(i,j) = i==j ? one : zero;
              });
          });
      }
      return 0;
    }
  };

}
#endif
