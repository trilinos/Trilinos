#ifndef __TACHO_SYMMETRIZE_ON_DEVICE_HPP__
#define __TACHO_SYMMETRIZE_ON_DEVICE_HPP__


/// \file  Tacho_Symmetrize_OnDevice.hpp
/// \brief Symmetrize a matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<>
  struct Symmetrize<Uplo::Upper,Algo::OnDevice> {
    template<typename MemberType,
             typename ViewTypeA>
    inline
    static int
    invoke(MemberType &member,
           const ViewTypeA &A) {
      typedef typename ViewTypeA::non_const_value_type value_type;

      static_assert(ViewTypeA::rank == 2,"A is not rank 2 view.");

      const ordinal_type
        m = A.extent(0),
        n = A.extent(1);

      if (m == n) {
        if (A.span() > 0) {
          using exec_space = MemberType;
          using team_policy_type = Kokkos::TeamPolicy<exec_space>;
          
          const auto exec_instance = member;
          const auto policy = team_policy_type(exec_instance, n, Kokkos::AUTO);
          Kokkos::parallel_for
            (policy, KOKKOS_LAMBDA(const typename team_policy_type::member_type &member) {
              const ordinal_type j = member.league_rank();
              Kokkos::parallel_for
                (Kokkos::TeamVectorRange(member, m),
                 [&, A, j](const ordinal_type &i) { // Value capture is a workaround for cuda + gcc-7.2 compiler bug w/c++14
                  A(i,j) = i < j ? A(j,i) : A(i,j);
                });
            });
        }
      } else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "A is not a square matrix");
      }
      return 0;
    }
  };

}
#endif
